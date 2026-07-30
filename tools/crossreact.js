#!/usr/bin/env node
/* crossreact.js — REPRODUCIBLE in-silico cross-reactivity (specificity) screen
 * for the LAMPrime paper. Uses the ACTUAL LAMPrime algorithm: the pure
 * functions below are copied VERBATIM from ../app.js (line numbers cited per
 * block). No numbers are invented — every result printed comes from running
 * LAMPrime's own design + specificity screen.
 *
 * Pipeline (mirrors the browser app exactly):
 *   1. Load target FASTA, run designLAMP() with the app's default params
 *      (app.js `defaults`, lines 265-298 + readParams defaults, 526-550),
 *      take the top-ranked set sets[0] = "Set #1" (our designed primers).
 *   2. Load a background FASTA, run the specificity screen of Set #1's primers
 *      against it, exactly as the btnScreen handler (app.js 938-952) +
 *      renderSpec (840-865) do: effective seed k = max(8, min(seed, minLen)),
 *      default seed=13, default maxMM=2.
 *   3. Print, per primer (F3,F2,F1c,B1c,B2,B3,LF,LB): flagged?, seed length,
 *      best gap-free match (mismatches, length), and 1-based position+strand.
 *
 * Usage:  node tools/crossreact.js
 */
'use strict';
const fs = require('fs');
const path = require('path');

// ===================================================================
// VERBATIM from app.js — Tm NN (SantaLucia 1998) + salts. Lines 389-413.
// ===================================================================
const COMP = { A:'T', T:'A', G:'C', C:'G', N:'N' };
const comp = b => COMP[b] || 'N';
const revComp = s => { let o=''; for (let i=s.length-1;i>=0;i--) o+=comp(s[i]); return o; };
const gcPct = s => { const n=s.length||1; return 100*((s.match(/[GC]/g)||[]).length)/n; };

const NN_DH = {AA:-7.9,AT:-7.2,TA:-7.2,CA:-8.5,GT:-8.4,CT:-7.8,GA:-8.2,CG:-10.6,GC:-9.8,GG:-8.0,AC:-8.4,AG:-7.8,TC:-8.2,TG:-8.5,TT:-7.9,CC:-8.0};
const NN_DS = {AA:-22.2,AT:-20.4,TA:-21.3,CA:-22.7,GT:-22.4,CT:-21.0,GA:-22.2,CG:-27.2,GC:-24.4,GG:-19.9,AC:-22.4,AG:-21.0,TC:-22.2,TG:-22.7,TT:-22.2,CC:-19.9};
const R_GAS = 1.987;
const TM_OLIGO_NM = 50;

// #6: Mg2+ livre por quelacao 1:1 com dNTP (Ka=3e4 M^-1), unificado com app.js/mgspec. Conc. em MOLAR.
function freeMgM(totMgM, dntpM, Ka=3e4){ const B=Ka*dntpM - Ka*totMgM + 1; return (-(B)+Math.sqrt(B*B+4*Ka*totMgM))/(2*Ka); }

function tmNN(seq, naMM, mgMM, dntpMM, oligoNM) {
  seq = (seq||'').toUpperCase().replace(/[^ATGC]/g,'');
  const N = seq.length; if (N < 2) return NaN;
  let dH=0, dS=0;
  for (let i=0;i<N-1;i++){ const d=seq.substr(i,2); if(NN_DH[d]===undefined) return NaN; dH+=NN_DH[d]; dS+=NN_DS[d]; }
  const init = b => (b==='G'||b==='C') ? [0.1,-2.8] : [2.3,4.1];
  const [h5,s5]=init(seq[0]); const [h3,s3]=init(seq[N-1]); dH+=h5+h3; dS+=s5+s3;
  const mgFree = freeMgM((mgMM||0)/1000, (dntpMM||0)/1000) * 1000; // #6: quelacao (M -> mM)
  const naEq = Math.max(1e-3, (naMM||0) + 120*Math.sqrt(Math.max(0,mgFree))) / 1000;
  const dS_salt = dS + 0.368*(N-1)*Math.log(naEq);
  const C = (oligoNM||TM_OLIGO_NM)*1e-9;
  return (dH*1000) / (dS_salt + R_GAS*Math.log(C/4)) - 273.15;
}

// VERBATIM app.js 416-425
function linguisticComplexity(s) {
  s=(s||'').toUpperCase(); const N=s.length; if (N<3) return 100;
  let prod=1, terms=0;
  for (let k=1;k<=3;k++){
    const max=Math.min(Math.pow(4,k), N-k+1); if (max<=0) continue;
    const set=new Set(); for (let i=0;i+k<=N;i++) set.add(s.substr(i,k));
    prod *= set.size/max; terms++;
  }
  return 100*Math.pow(prod, terms?1/terms:1);
}

// VERBATIM app.js 429-436
function endDeltaG(seq, end) {
  seq=(seq||'').toUpperCase().replace(/[^ATGC]/g,'');
  const w = (end==='5') ? seq.slice(0,6) : seq.slice(-6);
  if (w.length<2) return 0;
  let dH=0,dS=0;
  for (let i=0;i<w.length-1;i++){ const d=w.substr(i,2); if(NN_DH[d]===undefined) continue; dH+=NN_DH[d]; dS+=NN_DS[d]; }
  return dH - 310.15*(dS/1000);
}

// VERBATIM app.js 441-498 (secondary structures)
let REACT_C = 63;
const TREACT_K = () => (REACT_C + 273.15);
function segDG(seg) {
  seg = (seg||'').toUpperCase(); const N = seg.length; if (N < 2) return 0;
  let dH=0, dS=0;
  for (let i=0;i<N-1;i++){ const d=seg.substr(i,2); if(NN_DH[d]===undefined) continue; dH+=NN_DH[d]; dS+=NN_DS[d]; }
  const ini = b => (b==='G'||b==='C') ? [0.1,-2.8] : [2.3,4.1];
  const [h5,s5]=ini(seg[0]); const [h3,s3]=ini(seg[N-1]); dH+=h5+h3; dS+=s5+s3;
  return dH - TREACT_K()*dS/1000;
}
function dimerDG(a, b) {
  a=(a||'').toUpperCase(); const B=(b||'').toUpperCase().split('').reverse().join('');
  let best=0;
  for (let off=-(a.length-1); off<B.length; off++){
    let i=Math.max(0,-off);
    while (i<a.length && i+off<B.length){
      if (i+off>=0 && a[i]===comp(B[i+off])){
        let j=i, seg='';
        while (j<a.length && j+off>=0 && j+off<B.length && a[j]===comp(B[j+off])){ seg+=a[j]; j++; }
        if (seg.length>=3){ const dg=segDG(seg); if (dg<best) best=dg; }
        i=j;
      } else i++;
    }
  }
  return best;
}
function hairpinDG(a) {
  a=(a||'').toUpperCase(); let best=0;
  for (let i=0;i<a.length;i++){
    for (let j=a.length-1;j>i+3;j--){
      if (a[i]!==comp(a[j])) continue;
      let k=0, seg='';
      while (i+k < j-k && (j-k)-(i+k)-1>=3 && a[i+k]===comp(a[j-k])){ seg+=a[i+k]; k++; }
      if (seg.length>=3){ const dg=segDG(seg); if (dg<best) best=dg; }
    }
  }
  return best;
}
function structureScan(primers, hpTh, diTh) {
  let pen=0; const warns=[]; let worstHp=0, worstDi=0;
  for (const p of primers){
    const hp=hairpinDG(p.seq); if (hp<worstHp) worstHp=hp;
    if (hp<hpTh){ pen += (hpTh-hp); warns.push({t:'hairpin', n:p.name, v:hp.toFixed(1)}); }
    const sd=dimerDG(p.seq, p.seq); if (sd<worstDi) worstDi=sd;
    if (sd<diTh){ pen += (diTh-sd); warns.push({t:'selfdimer', n:p.name, v:sd.toFixed(1)}); }
  }
  for (let a=0;a<primers.length;a++) for (let b=a+1;b<primers.length;b++){
    const cd=dimerDG(primers[a].seq, primers[b].seq); if (cd<worstDi) worstDi=cd;
    if (cd<diTh){ pen += (diTh-cd)*0.8; warns.push({t:'dimer', n:`${primers[a].name}/${primers[b].name}`, v:cd.toFixed(1)}); }
  }
  return { penalty: pen*3, warns, worstHairpin: worstHp, worstDimer: worstDi };
}

// VERBATIM app.js 500-524
const inRange = (x, lo, hi) => x>=lo && x<=hi;
function bestByStart(S, start, lenMin, lenMax, tmLo, tmHi, P) {
  const mid=(tmLo+tmHi)/2; let best=null;
  for (let L=lenMin; L<=lenMax; L++){
    if (start+L>S.length) break;
    const seg=S.substr(start,L); const tm=tmNN(seg,P.na,P.mg,P.dntp,TM_OLIGO_NM); if (isNaN(tm)) continue;
    const dev=Math.abs(tm-mid);
    if (!best || dev<best.dev) best={start,end:start+L-1,len:L,seg,tm,gc:gcPct(seg),lc:linguisticComplexity(seg),dev};
  }
  return best;
}
function bestByEnd(S, end, lenMin, lenMax, tmLo, tmHi, P) {
  const mid=(tmLo+tmHi)/2; let best=null;
  for (let L=lenMin; L<=lenMax; L++){
    const start=end-L+1; if (start<0) break;
    const seg=S.substr(start,L); const tm=tmNN(seg,P.na,P.mg,P.dntp,TM_OLIGO_NM); if (isNaN(tm)) continue;
    const dev=Math.abs(tm-mid);
    if (!best || dev<best.dev) best={start,end,len:L,seg,tm,gc:gcPct(seg),lc:linguisticComplexity(seg),dev};
  }
  return best;
}

// VERBATIM app.js 552-675 (designLAMP)
function designLAMP(S, P) {
  const N=S.length;
  if (N<120) return {errorKey:'alert.shortSeq'};
  REACT_C = (P.tReact!=null) ? P.tReact : 63;
  const innerLenMin=18, innerLenMax=26;
  const geom={ f2b2:P.d_F2B2, gap:P.gapF3F2, loopF2F1:P.d_F2F1, f1b1:P.d_F1B1 };
  const blockDev = b => (b.F1?b.F1.dev:b.B1.dev)+(b.F2?b.F2.dev:b.B2.dev)+(b.F3?b.F3.dev:b.B3.dev);

  let Fblocks=[];
  for (let f2s=0; f2s<N; f2s++) {
    const F2=bestByStart(S,f2s,P.f3Min,P.f3Max,P.outerLo,P.outerHi,P);
    if (!F2 || !inRange(F2.tm,P.outerLo,P.outerHi)) continue;
    let F1=null;
    for (let d=P.d_F2F1[0]; d<=P.d_F2F1[1]; d++){
      const c=bestByStart(S,f2s+d,innerLenMin,innerLenMax,P.innerLo,P.innerHi,P);
      if (c && inRange(c.tm,P.innerLo,P.innerHi) && (!F1||c.dev<F1.dev)) F1=c;
    }
    if (!F1) continue;
    let F3=null;
    for (let g=P.gapF3F2[0]; g<=P.gapF3F2[1]; g++){
      const f3e=f2s-1-g; if (f3e<P.f3Min-1) break;
      const c=bestByEnd(S,f3e,P.f3Min,P.f3Max,P.outerLo,P.outerHi,P);
      if (c && inRange(c.tm,P.outerLo,P.outerHi) && (!F3||c.dev<F3.dev)) F3=c;
    }
    if (!F3) continue;
    Fblocks.push({F3,F2,F1});
  }
  Fblocks.sort((a,b)=>blockDev(a)-blockDev(b)); Fblocks=Fblocks.slice(0,250);

  let Bblocks=[];
  for (let b1s=0; b1s<N; b1s++) {
    const B1=bestByStart(S,b1s,innerLenMin,innerLenMax,P.innerLo,P.innerHi,P);
    if (!B1 || !inRange(B1.tm,P.innerLo,P.innerHi)) continue;
    let B2=null;
    for (let d=P.d_B2B1[0]; d<=P.d_B2B1[1]; d++){
      const c=bestByStart(S,b1s+d,P.f3Min,P.f3Max,P.outerLo,P.outerHi,P);
      if (c && inRange(c.tm,P.outerLo,P.outerHi) && (!B2||c.dev<B2.dev)) B2=c;
    }
    if (!B2) continue;
    let B3=null;
    for (let g=P.gapF3F2[0]; g<=P.gapF3F2[1]; g++){
      const c=bestByStart(S,B2.end+1+g,P.f3Min,P.f3Max,P.outerLo,P.outerHi,P);
      if (c && inRange(c.tm,P.outerLo,P.outerHi) && (!B3||c.dev<B3.dev)) B3=c;
    }
    if (!B3) continue;
    Bblocks.push({B1,B2,B3});
  }
  Bblocks.sort((a,b)=>blockDev(a)-blockDev(b)); Bblocks=Bblocks.slice(0,250);

  const sets=[];
  for (const F of Fblocks){
    for (const B of Bblocks){
      if (!(F.F1.end < B.B1.start)) continue;
      const dF1B1 = B.B1.start - F.F1.end - 1;
      if (!inRange(dF1B1, P.d_F1B1[0], P.d_F1B1[1])) continue;
      const dF2B2 = B.B2.end - F.F2.start + 1;
      if (!inRange(dF2B2, P.d_F2B2[0], P.d_F2B2[1])) continue;
      if (!(F.F3.start<F.F2.start && F.F2.end<F.F1.start && B.B1.end<B.B2.start && B.B2.end<B.B3.start)) continue;

      const lfRegion = S.substring(F.F2.end+1, F.F1.start);
      const lbRegion = S.substring(B.B1.end+1, B.B2.start);
      let LF=null, LB=null;
      if (P.loops!=='off'){
        if (lfRegion.length>=P.lfMin){ const c=bestByStart(lfRegion,0,P.lfMin,Math.min(P.lfMax,lfRegion.length),P.innerLo,P.innerHi,P); if(c) LF={seq:revComp(c.seg),tm:c.tm,gc:c.gc,len:c.len}; }
        if (lbRegion.length>=P.lfMin){ const c=bestByStart(lbRegion,0,P.lfMin,Math.min(P.lfMax,lbRegion.length),P.innerLo,P.innerHi,P); if(c) LB={seq:c.seg,tm:c.tm,gc:c.gc,len:c.len}; }
      }

      const F1c=revComp(F.F1.seg), F2=F.F2.seg, B1c=B.B1.seg, B2=revComp(B.B2.seg);
      const FIP=F1c+P.loopLinker+F2, BIP=B1c+P.loopLinker+B2;
      const F3=F.F3.seg, B3=revComp(B.B3.seg);

      const midI=(P.innerLo+P.innerHi)/2, midO=(P.outerLo+P.outerHi)/2;
      const pInner=(Math.abs(F.F1.tm-midI)+Math.abs(B.B1.tm-midI))/Math.max(1e-6,P.limF1cB1c);
      const pF2B2=(Math.abs(F.F2.tm-midO)+Math.abs(B.B2.tm-midO))/Math.max(1e-6,P.limF2B2);
      const pF3B3=(Math.abs(F.F3.tm-midO)+Math.abs(B.B3.tm-midO))/Math.max(1e-6,P.limF3B3);
      const segs=[['F3',F.F3],['F2',F.F2],['F1c',F.F1],['B1c',B.B1],['B2',B.B2],['B3',B.B3]];
      let gcPen=0,lcPen=0,warns=[];
      for (const [name,info] of segs){
        if (!inRange(info.gc,P.gcLo,P.gcHi)){ gcPen+=1; warns.push({t:'gc', n:name, v:info.gc.toFixed(0)}); }
        if (info.lc<P.lcMin){ lcPen+=1; warns.push({t:'lc', n:name, v:info.lc.toFixed(0)}); }
      }
      const fipOk=inRange(FIP.length,P.fipMin,P.fipMax), bipOk=inRange(BIP.length,P.bipMin,P.bipMax);
      if (!fipOk) warns.push({t:'len', n:'FIP', v:FIP.length}); if (!bipOk) warns.push({t:'len', n:'BIP', v:BIP.length});
      let dgPen=0;
      [['F3',F3,'3'],['F2',F2,'3'],['B3',B3,'3'],['B2',B2,'3']].forEach(([n,s,e])=>{ const g=endDeltaG(s,e); if(g>P.endDg){dgPen++; warns.push({t:'dg', n, end:e, v:g.toFixed(1)});} });
      [['F1c',F1c,'5'],['B1c',B1c,'5']].forEach(([n,s,e])=>{ const g=endDeltaG(s,e); if(g>P.endDg){dgPen++; warns.push({t:'dg', n, end:e, v:g.toFixed(1)});} });
      const penalty=(pInner+pF2B2+pF3B3)*100 + gcPen*5 + lcPen*5 + dgPen*4 + (fipOk?0:8) + (bipOk?0:8);

      sets.push({
        penalty,
        ampliconF3B3: B.B3.end - F.F3.start + 1, ampliconF2B2: dF2B2, dF1B1,
        coords:{F3:[F.F3.start+1,F.F3.end+1],F2:[F.F2.start+1,F.F2.end+1],F1:[F.F1.start+1,F.F1.end+1],B1:[B.B1.start+1,B.B1.end+1],B2:[B.B2.start+1,B.B2.end+1],B3:[B.B3.start+1,B.B3.end+1]},
        F3:{seq:F3,tm:F.F3.tm,gc:F.F3.gc}, B3:{seq:B3,tm:B.B3.tm,gc:B.B3.gc},
        F2:{seq:F2,tm:F.F2.tm,gc:F.F2.gc}, B2:{seq:B2,tm:B.B2.tm,gc:B.B2.gc},
        F1c:{seq:F1c,tm:F.F1.tm,gc:F.F1.gc}, B1c:{seq:B1c,tm:B.B1.tm,gc:B.B1.gc},
        FIP:{seq:FIP,len:FIP.length}, BIP:{seq:BIP,len:BIP.length},
        LF, LB, warns
      });
    }
  }
  sets.sort((a,b)=>a.penalty-b.penalty);
  const seen=new Set(), uniq=[];
  for (const s of sets){ const k=s.coords.F2[0]+'|'+s.coords.B2[1]; if(seen.has(k))continue; seen.add(k); uniq.push(s); if(uniq.length>=Math.min(P.maxSets,100))break; }
  for (const s of uniq){
    const prims=[{name:'F3',seq:s.F3.seq},{name:'B3',seq:s.B3.seq},{name:'FIP',seq:s.FIP.seq},{name:'BIP',seq:s.BIP.seq}];
    if (s.LF) prims.push({name:'LF',seq:s.LF.seq});
    if (s.LB) prims.push({name:'LB',seq:s.LB.seq});
    const st=structureScan(prims, P.hpTh, P.dimerTh);
    s.penalty += st.penalty;
    s.struct={ hairpin:st.worstHairpin, dimer:st.worstDimer };
    if (st.warns.length) s.warns=(s.warns||[]).concat(st.warns);
  }
  uniq.sort((a,b)=>a.penalty-b.penalty);
  uniq.forEach((s,i)=>s.rank=i+1);
  return { sets:uniq, totalCandidates:sets.length, geom, params:P, targetLen:N, gc:Math.round(gcPct(S)), tReact:P.tReact };
}

// ===================================================================
// VERBATIM from app.js — specificity screen. Lines 802-839.
// ===================================================================
function sanitizeDNA(seq) { return (seq || '').toUpperCase().replace(/U/g, 'T').replace(/[^ATGC]/g, ''); }
function buildKmerIndex(bg, k) {
  const idx=new Map();
  for (let i=0; i+k<=bg.length; i++){ const key=bg.substr(i,k); let a=idx.get(key); if(!a){a=[];idx.set(key,a);} a.push(i); }
  return idx;
}
function countMM(pat, bg, start, cap) {
  const Lp=pat.length; if (start<0 || start+Lp>bg.length) return cap+1;
  let mm=0; for (let i=0;i<Lp;i++){ if (bg[start+i]!==pat[i]){ if (++mm>cap) return mm; } }
  return mm;
}
function screenPrimerOnBg(primer, bg, idx, k, maxMM) {
  const Lp=primer.length; if (Lp<k) return null;
  let best=null;
  const consider=(ws, pat, strand)=>{ const mm=countMM(pat,bg,ws,maxMM); if (mm<=maxMM && (!best||mm<best.mm)) best={mm, start:ws, end:ws+Lp, strand}; };
  const occB=idx.get(primer.slice(Lp-k));
  if (occB) for (const pos of occB) consider(pos-(Lp-k), primer, '+');
  const rc=revComp(primer);
  const occA=idx.get(rc.slice(0,k));
  if (occA) for (const pos of occA) consider(pos, rc, '−');
  return best;
}
function primersFromSet(set) {
  const ps=[{name:'F3',seq:set.F3.seq},{name:'F2',seq:set.F2.seq},{name:'F1c',seq:set.F1c.seq},
            {name:'B1c',seq:set.B1c.seq},{name:'B2',seq:set.B2.seq},{name:'B3',seq:set.B3.seq}];
  if (set.LF) ps.push({name:'LF',seq:set.LF.seq});
  if (set.LB) ps.push({name:'LB',seq:set.LB.seq});
  return ps;
}

// ===================================================================
// app.js `defaults` (265-298) -> readParams output shape (526-550).
// These are the EXACT app defaults, used unmodified for the design.
// ===================================================================
function defaultParams() {
  return {
    innerLo:64, innerHi:66, outerLo:59, outerHi:61,
    f3Min:18, f3Max:23, lfMin:18, lfMax:23,
    fipMin:42, fipMax:52, bipMin:42, bipMax:52,
    gcLo:40, gcHi:60, loops:'auto', lcMin:75,
    na:50, mg:8, dntp:1.4,
    d_F2F1:[40,60], d_F1B1:[0,60], d_F2B2:[120,180], d_B2B1:[40,60], gapF3F2:[0,60],
    loopLinker:'TTTTTT',
    limF1cB1c:15, limF2B2:20, limF3B3:15,
    maxSets:300, tReact:63, hpTh:-3.0, dimerTh:-3.0, endDg:-4.0
  };
}

// ---- harness helpers ----
function loadFasta(p) {
  const text = fs.readFileSync(p, 'utf-8');
  const seq = text.split(/\r?\n/).filter(l => !l.startsWith('>')).join('').trim();
  const header = (text.split(/\r?\n/).find(l => l.startsWith('>')) || '').slice(1).trim();
  return { header, seq: sanitizeDNA(seq) };
}

// Replicates the btnScreen handler (app.js 938-952) + the effective seed
// computed in renderSpec (848-851): k = max(8, min(seed, minPrimerLen)).
function runScreen(set, bgName, bgSeq, seed, maxMM) {
  const primers = primersFromSet(set);
  const minLen = Math.min.apply(null, primers.map(p => p.seq.length));
  const seedReq = Math.max(8, Math.min(20, seed));          // app clamps seed 8..20 (line 942)
  const k = Math.max(8, Math.min(seedReq, minLen));         // effective seed (renderSpec 849)
  const idx = buildKmerIndex(bgSeq, k);
  const rows = [];
  let hits = 0;
  for (const p of primers) {
    const m = screenPrimerOnBg(p.seq, bgSeq, idx, k, maxMM);
    if (m) { hits++; rows.push({ name:p.name, seq:p.seq, len:p.seq.length, flagged:true, seed:k, mm:m.mm, matchLen:p.seq.length, pos:`${m.start+1}-${m.end}`, strand:m.strand }); }
    else rows.push({ name:p.name, seq:p.seq, len:p.seq.length, flagged:false, seed:k });
  }
  return { bgName, bgLen:bgSeq.length, k, maxMM, hits, total:primers.length, rows };
}

function printSet(label, des) {
  const s = des.sets[0];
  console.log(`\n--- ${label}: LAMPrime Set #1 (penalty ${s.penalty.toFixed(1)}; ${des.sets.length} valid candidates of ${des.totalCandidates}) ---`);
  const ps = primersFromSet(s);
  for (const p of ps) console.log(`    ${p.name.padEnd(4)} 5'-${p.seq}-3'  (${p.seq.length} nt)`);
}

function printScreen(res) {
  console.log(`\n  Background: ${res.bgName} (${res.bgLen} nt) | effective seed k=${res.k} | maxMM=${res.maxMM}`);
  console.log(`  ${res.hits} of ${res.total} primers FLAGGED (possible 3'-anchored anneal)`);
  console.log('  primer | flagged | seed | mm | matchLen | position(1-based) | strand');
  for (const r of res.rows) {
    if (r.flagged) console.log(`    ${r.name.padEnd(4)} |   YES   |  ${r.seed}  | ${r.mm}  |   ${String(r.matchLen).padStart(3)}    | ${r.pos.padStart(13)} |  ${r.strand}`);
    else           console.log(`    ${r.name.padEnd(4)} |   no    |  ${r.seed}  | -  |    -     |       -       |  -`);
  }
}

// ===================================================================
// MAIN — two scenarios + own-target sanity controls
// ===================================================================
const DATA = path.join(__dirname, 'data');
const SEED = 13, MAXMM = 2; // app defaults (app.js 942-943)

function scenario(title, targetFile, bgFile, ownLabel) {
  console.log('\n========================================================================');
  console.log(title);
  console.log('========================================================================');
  const tgt = loadFasta(path.join(DATA, targetFile));
  console.log(`Target FASTA: ${targetFile}`);
  console.log(`  header: ${tgt.header}`);
  console.log(`  length: ${tgt.seq.length} nt, GC ${gcPct(tgt.seq).toFixed(1)}%`);
  const des = designLAMP(tgt.seq, defaultParams());
  if (des.errorKey || !des.sets || !des.sets.length) {
    console.log(`  DESIGN FAILED: ${des.errorKey || 'no sets satisfied constraints'}`);
    return;
  }
  printSet(ownLabel, des);
  const set1 = des.sets[0];

  // sanity control: screen Set #1 against its OWN target (expect 100% flagged)
  const ctrl = runScreen(set1, `${ownLabel} OWN TARGET (${targetFile})`, tgt.seq, SEED, MAXMM);
  console.log('\n  [SANITY CONTROL] Set #1 vs its own target:');
  printScreen(ctrl);

  // cross-reactivity: screen Set #1 against the related-organism background
  const bg = loadFasta(path.join(DATA, bgFile));
  const cross = runScreen(set1, `${bg.header.slice(0,55)} [${bgFile}]`, bg.seq, SEED, MAXMM);
  console.log('\n  [CROSS-REACTIVITY] Set #1 vs related-organism background:');
  printScreen(cross);

  // ENSEMBLE cross-reactivity: the top-ranked designed sets are near-tied in penalty,
  // so the single Set #1 count is fragile (a ~0.1 C model change reorders them). Screen the
  // top-N sets to report a ROBUST RANGE instead of one knife-edge number.
  const ens = ensembleCrossReact(des, bg.seq, Math.min(10, des.sets.length), SEED, MAXMM);
  console.log(`\n  [CROSS-REACTIVITY, ENSEMBLE] Top-${ens.n} near-tied designs vs ${bgFile}:`);
  console.log(`    flagged primers /8 per set (sorted): [${ens.counts.join(', ')}]`);
  console.log(`    penalties of these sets: [${ens.penalties.join(', ')}]`);
  console.log(`    ROBUST RANGE: ${ens.min}–${ens.max} of 8 (median ${ens.median}) across the ${ens.n} best-ranked designs`);
  console.log(`    union of primers ever flagged: {${ens.unionNames.join(', ')}}`);
}

// Screen the top-N designed sets against a background; returns the distribution of
// flagged-primer counts (robust to the near-tie in penalty ranking). See §3.3 of the paper.
function ensembleCrossReact(des, bgSeq, n, seed, maxMM) {
  const counts = [], penalties = [], union = new Set();
  for (let i = 0; i < n; i++) {
    const r = runScreen(des.sets[i], '', bgSeq, seed, maxMM);
    counts.push(r.hits);
    penalties.push(Number(des.sets[i].penalty.toFixed(1)));
    for (const row of r.rows) if (row.flagged) union.add(row.name);
  }
  const sorted = counts.slice().sort((a, b) => a - b);
  return {
    n, counts: sorted, penalties,
    min: sorted[0], max: sorted[sorted.length - 1],
    median: sorted[Math.floor(sorted.length / 2)],
    unionNames: [...union],
  };
}

if (require.main === module) {
  console.log('LAMPrime in-silico cross-reactivity screen — reproducible harness');
  console.log('Algorithm: VERBATIM pure functions from app.js (design + 3\'-seed specificity screen)');
  console.log(`Params: app defaults | seed=${SEED} (effective k clamped to min primer len), maxMM=${MAXMM}`);

  scenario(
    'SCENARIO 1 — A. marginale msp1b design  vs  Anaplasma centrale background',
    'amarginale_msp1b_synthetic.fasta',
    'acentrale_genome_CP001759.fasta',
    'A. marginale msp1b'
  );

  scenario(
    'SCENARIO 2 — B. bovis 18S design  vs  Babesia bigemina 18S background',
    'bbovis_18s_L19077.fasta',
    'bbigemina_18s_KP710228.fasta',
    'B. bovis 18S'
  );

  console.log('\nDONE.');
}

// Exportado para reuso (ex.: desenho de novos alvos como P. falciparum 18S) sem re-executar os cenarios.
module.exports = { designLAMP, defaultParams, loadFasta, primersFromSet, runScreen, ensembleCrossReact, tmNN, freeMgM, gcPct, DATA, SEED, MAXMM };
