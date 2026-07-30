#!/usr/bin/env node
/* mgspec.mjs — offline validation of the LAMPrime Mg2+-vs-specificity Tm path.
 *
 * Mirrors the math added to app.js (kept in sync by hand; app.js is a browser
 * IIFE and is not importable). Exercises:
 *   - freeMgM: 1:1 dNTP chelation equilibrium (Ka=3e4 M^-1)
 *   - owczarzy2008: divalent Mg2+ correction (Owczarzy 2008)
 *   - tmDuplexMg: NN duplex Tm with internal-mismatch params (Allawi&SantaLucia
 *     1997/1998; Peyret 1999), separate from the published tmNN path.
 *
 * Sanity checks (assert): perfect duplex Tm > single-mismatch Tm; raising total
 * Mg2+ raises Tm; free Mg2+ < total Mg2+ in the presence of dNTP.
 *
 * Run:  node tools/mgspec.mjs
 * Refs: SantaLucia 1998; Allawi&SantaLucia 1997; Peyret 1999; Owczarzy 2008.
 */

// ---- Watson-Crick NN (SantaLucia 1998 unified) ----
const NN_DH = {AA:-7.9,AT:-7.2,TA:-7.2,CA:-8.5,GT:-8.4,CT:-7.8,GA:-8.2,CG:-10.6,GC:-9.8,GG:-8.0,AC:-8.4,AG:-7.8,TC:-8.2,TG:-8.5,TT:-7.9,CC:-8.0};
const NN_DS = {AA:-22.2,AT:-20.4,TA:-21.3,CA:-22.7,GT:-22.4,CT:-21.0,GA:-22.2,CG:-27.2,GC:-24.4,GG:-19.9,AC:-22.4,AG:-21.0,TC:-22.2,TG:-22.7,TT:-22.2,CC:-19.9};
const R_GAS = 1.987;
const TM_OLIGO_NM = 50;

// ---- Internal single-mismatch NN (see app.js MM_IMM for citations) ----
const MM_IMM = {
'AG/TT':[1.0,0.9],'AT/TG':[-2.5,-8.3],'CG/GT':[-4.1,-11.7],'CT/GG':[-2.8,-8.0],'GG/CT':[3.3,10.4],'GG/TT':[5.8,16.3],'GT/CG':[-4.4,-12.3],'GT/TG':[4.1,9.5],'TG/AT':[-0.1,-1.7],'TG/GT':[-1.4,-6.2],'TT/AG':[-1.3,-5.3],
'AA/TG':[-0.6,-2.3],'AG/TA':[-0.7,-2.3],'CA/GG':[-0.7,-2.3],'CG/GA':[-4.0,-13.2],'GA/CG':[-0.6,-1.0],'GG/CA':[0.5,3.2],'TA/AG':[0.7,0.7],'TG/AA':[3.0,7.4],
'AC/TT':[0.7,0.2],'AT/TC':[-1.2,-6.2],'CC/GT':[-0.8,-4.5],'CT/GC':[-1.5,-6.1],'GC/CT':[2.3,5.4],'GT/CC':[5.2,13.5],'TC/AT':[1.2,0.7],'TT/AC':[1.0,0.7],
'AA/TC':[2.3,4.6],'AC/TA':[5.3,14.6],'CA/GC':[1.9,3.7],'CC/GA':[0.6,-0.6],'GA/CC':[5.2,14.2],'GC/CA':[-0.7,-3.8],'TA/AC':[3.4,8.0],'TC/AA':[7.6,20.2],
'AA/TA':[1.2,1.7],'CA/GA':[-0.9,-4.2],'GA/CA':[-2.9,-9.8],'TA/AA':[4.7,12.9],
'AC/TC':[0.0,-4.4],'CC/GC':[-1.5,-7.2],'GC/CC':[3.6,8.9],'TC/AC':[6.1,16.4],
'AG/TG':[-3.1,-9.5],'CG/GG':[-4.9,-15.3],'GG/CG':[-6.0,-15.8],'TG/AG':[1.6,3.6],
'AT/TT':[-2.7,-10.8],'CT/GT':[-5.0,-15.8],'GT/CT':[-2.2,-8.4],'TT/AT':[0.2,-1.5]
};

const COMP = { A:'T', T:'A', G:'C', C:'G', N:'N' };
const comp = b => COMP[b] || 'N';
const revComp = s => { let o=''; for (let i=s.length-1;i>=0;i--) o+=comp(s[i]); return o; };
const wcBottom = top => top.toUpperCase().split('').map(comp).join(''); // aligned 3'->5'

function freeMgM(totMgM, dntpM, Ka=3e4){ const B=Ka*dntpM - Ka*totMgM + 1; return (-(B)+Math.sqrt(B*B+4*Ka*totMgM))/(2*Ka); }

function owczarzy2008(tmK1M, freeMg, monM, fGC, Nbp){
  const ln=Math.log; let a=3.92,b=-0.911,c=6.26,d=1.42,e=-48.2,f=52.5,g=8.31;
  const R = Math.sqrt(Math.max(freeMg,1e-9))/Math.max(monM,1e-9);
  if (R < 0.22){ return tmK1M-273.15 + 0; }
  if (R < 6.0){ a=3.92*(0.843-0.352*Math.sqrt(monM)*ln(monM)); d=1.42*(1.279-4.03e-3*ln(monM)-8.03e-3*ln(monM)**2); g=8.31*(0.486-0.258*ln(monM)+5.25e-3*ln(monM)**3); }
  const lm=ln(freeMg);
  const corr=(a + b*lm + fGC*(c + d*lm) + (1/(2*(Nbp-1)))*(e + f*lm + g*lm*lm))*1e-5;
  const tmKc = 1/(1/tmK1M + corr);
  return tmKc - 273.15;
}

function tmDuplexMg(top, bottom, opt){
  top=top.toUpperCase(); bottom=bottom.toUpperCase();
  const N=Math.min(top.length, bottom.length); if (N<2) return {tm:NaN, flag:true};
  const isWC=(x,y)=>(x==='A'&&y==='T')||(x==='T'&&y==='A')||(x==='G'&&y==='C')||(x==='C'&&y==='G');
  let dH=0,dS=0,flag=false,gcMatch=0;
  for (let i=0;i<N;i++){ if (isWC(top[i],bottom[i])) gcMatch += (top[i]==='G'||top[i]==='C')?1:0; }
  for (let i=0;i<N-1;i++){
    const a=top[i],b=top[i+1],w=bottom[i],z=bottom[i+1];
    if (isWC(a,w)&&isWC(b,z)){ const d=a+b; if (NN_DH[d]===undefined){flag=true;continue;} dH+=NN_DH[d]; dS+=NN_DS[d]; }
    else { const key=a+b+'/'+w+z; let p=MM_IMM[key]; if(!p){ p=MM_IMM[z+w+'/'+b+a]; }
      if (p){ dH+=p[0]; dS+=p[1]; } else { dS+=12; flag=true; } }
  }
  const ini=(x,y)=>isWC(x,y)?((x==='G'||x==='C')?[0.1,-2.8]:[2.3,4.1]):[0,0];
  const [h5,s5]=ini(top[0],bottom[0]); const [h3,s3]=ini(top[N-1],bottom[N-1]); dH+=h5+h3; dS+=s5+s3;
  const C=((opt&&opt.oligoNM)||TM_OLIGO_NM)*1e-9;
  const denom=dS+R_GAS*Math.log(C/4); if (!(denom<0)) return {tm:NaN,flag:true};
  const tmK1M=(dH*1000)/denom;
  const totMgM=((opt&&opt.mg!=null?opt.mg:8))/1000, dntpM=((opt&&opt.dntp!=null?opt.dntp:1.4))/1000;
  const monM=Math.max(1e-3,((opt&&opt.mon!=null?opt.mon:50))/1000);
  const tmC=owczarzy2008(tmK1M, freeMgM(totMgM,dntpM), monM, gcMatch/N, N);
  return {tm:tmC, flag};
}

// Best gap-free off-target window: slide `primer` over each background (both
// strands), find the window with FEWEST mismatches; return the target base
// aligned 3'->5' (alignedTargetBottom) for tmDuplexMg. VERBATIM mirror of
// app.js bestOffWindow (lines 619-639). O(primerLen*bgLen) brute force.
function bestOffWindow(primer, bgSeqs){
  primer=(primer||'').toUpperCase(); const Lp=primer.length; if (Lp<2) return null;
  let best=null;
  const scan=(seq, strand, name)=>{
    for (let s=0; s+Lp<=seq.length; s++){
      let mm=0; for (let i=0;i<Lp;i++){ if (seq[s+i]!==primer[i]) mm++; }
      if (!best || mm<best.nMismatch){
        let bottom=''; for (let i=0;i<Lp;i++) bottom+=comp(seq[s+i]);
        best={ alignedTargetBottom:bottom, nMismatch:mm, pos:s, strand, src:name };
        if (mm===0) return;
      }
    }
  };
  (bgSeqs||[]).forEach(bg=>{
    const fwd=bg.seq||bg; const nm=bg.name||'';
    scan(fwd, '+', nm);
    scan(revComp(fwd), '-', nm);
  });
  return best;
}

// ---- Mg sweep + suggestion (verbatim mirror of app.js computeMgSweep) ----
// Per cell we keep S=Tm_on−Tm_off (intrinsic discrimination, ~Mg-invariant) AND
// safety=T_rxn−Tm_off (off-target safety margin, Mg-DEPENDENT). "functional" =
// Tm_on inside the role-based LAMP window. Suggested column = lowest total-Mg
// column where every primer is functional AND safe (safety>0); else the column
// maximizing the count of functional&safe primers, lowest Mg on ties.
const inRange=(x,lo,hi)=>x>=lo&&x<=hi;
function lampTmWindow(name){
  if (/^(F3|B3|F2|B2)$/.test(name)) return [59,61];
  if (/^(F1c|B1c|LF|LB)$/.test(name)) return [64,66];
  return [59,66];
}
function computeMgSweep(primers, bgs, opt){
  const mgMin=opt.mgMin, mgMax=opt.mgMax, mgStep=opt.mgStep>0?opt.mgStep:0.5;
  const cols=[]; for (let mg=mgMin; mg<=mgMax+1e-9; mg+=mgStep) cols.push(Math.round(mg*1000)/1000);
  const offs = primers.map(p => bestOffWindow(p.seq, bgs));
  const rows = primers.map((p,i)=>{
    const win=lampTmWindow(p.name);
    const cells = cols.map(mg=>{
      const o={ mg, dntp:opt.dntp, mon:opt.mon, oligoNM:TM_OLIGO_NM };
      const on = tmDuplexMg(p.seq, wcBottom(p.seq), o);
      let off={tm:NaN, flag:true};
      if (offs[i]) off = tmDuplexMg(p.seq, offs[i].alignedTargetBottom, o);
      const S = (isNaN(on.tm)||isNaN(off.tm)) ? NaN : (on.tm - off.tm);
      const safety = isNaN(off.tm) ? Infinity : (opt.tRxn - off.tm);
      const functional = !isNaN(on.tm) && inRange(on.tm, win[0], win[1]);
      const safe = safety > 0;
      const valid = functional && safe;
      const freeMg = freeMgM(mg/1000, opt.dntp/1000);
      return { mg, freeMg, tmOn:on.tm, tmOff:off.tm, S, safety, functional, safe, valid };
    });
    return { name:p.name, seq:p.seq, off:offs[i], cells };
  });
  const colStats = cols.map((mg,j)=>{
    let minSafety=Infinity, minS=Infinity, allPass=true, nPass=0;
    rows.forEach(r=>{
      const c=r.cells[j];
      if (!isNaN(c.S)) minS=Math.min(minS,c.S);
      if (isFinite(c.safety)) minSafety=Math.min(minSafety,c.safety);
      if (c.functional && c.safe) nPass++; else allPass=false;
    });
    if (!isFinite(minS)) minS=NaN;
    if (!isFinite(minSafety)) minSafety=NaN;
    return { mg, freeMg:freeMgM(mg/1000, opt.dntp/1000), minS, minSafety, feasible:allPass, nPass, nValid:nPass };
  });
  let jStar=-1;
  for (let j=0;j<colStats.length;j++){ if (colStats[j].feasible){ jStar=j; break; } }
  if (jStar<0){ let bestN=-1; colStats.forEach((c,j)=>{ if (c.nPass>bestN){ bestN=c.nPass; jStar=j; } }); }
  return { cols, rows, colStats, jStar, feasible: jStar>=0 && colStats[jStar].feasible, opt };
}

// ---- ESM exports (single reproducible math source; importing does NOT run the self-test) ----
export { freeMgM, owczarzy2008, tmDuplexMg, bestOffWindow, computeMgSweep, lampTmWindow, wcBottom, revComp, comp, NN_DH, NN_DS, MM_IMM, R_GAS, TM_OLIGO_NM };

// Run the self-test only when invoked directly (node tools/mgspec.mjs), not on import.
import { fileURLToPath as _f } from 'node:url';
const _isMain = (typeof process!=='undefined' && process.argv[1] && _f(import.meta.url) === process.argv[1]);
if (!_isMain) { /* imported as a library */ }
else { runSelfTest(); }

function assert(cond, msg){ if (!cond){ console.error('FAIL:', msg); process.exitCode=1; } else console.log('ok  -', msg); }

function runSelfTest(){
// ---- Tests ----
const top = 'GTCAGTTACCGTAGCATTCGAT';                 // 22-mer
const perfect = wcBottom(top);
// single internal mismatch at position 10 (force a non-WC pair on the bottom)
const mmBottom = perfect.split('');
mmBottom[10] = mmBottom[10]==='A' ? 'G' : 'A';
const offBottom = mmBottom.join('');

const opt2 = { mg:2, dntp:1.6, mon:50 };
const opt8 = { mg:8, dntp:1.6, mon:50 };

const tmOn2 = tmDuplexMg(top, perfect, opt2).tm;
const tmOn8 = tmDuplexMg(top, perfect, opt8).tm;
const tmOff8 = tmDuplexMg(top, offBottom, opt8).tm;

console.log(`Tm(perfect, 2mM Mg)  = ${tmOn2.toFixed(2)} C`);
console.log(`Tm(perfect, 8mM Mg)  = ${tmOn8.toFixed(2)} C`);
console.log(`Tm(1-mismatch, 8mM)  = ${tmOff8.toFixed(2)} C`);
console.log(`free Mg @ 8mM tot, 1.6mM dNTP = ${(freeMgM(0.008,0.0016)*1000).toFixed(3)} mM`);

assert(Number.isFinite(tmOn8), 'perfect-match Tm is finite');
assert(tmOn8 > tmOff8 + 1, 'perfect-match Tm exceeds single-mismatch Tm by >1 C');
assert(tmOn8 > tmOn2, 'raising total Mg2+ (2->8 mM) raises Tm');
assert(freeMgM(0.008,0.0016) < 0.008, 'free Mg2+ < total Mg2+ when dNTP present');
assert(Math.abs(freeMgM(0.008,0.0) - 0.008) < 1e-9, 'free Mg2+ == total when no dNTP');

// ---- Mg sweep: safety metric (T_rxn − Tm_off) and suggestion (lowest functional&safe Mg) ----
// One synthetic primer (the on-target) vs a near-identical background (the off-target).
const swTop = 'GTCAGTTACCGTAGCATTCGAT';
const bgFasta = [{ name:'bg', seq: swTop.slice(0,10) + (swTop[10]==='A'?'G':'A') + swTop.slice(11) }]; // 1-mm window present
const swPrimers = [{ name:'F3', seq: swTop }];
const Trxn = 63;
const sw = computeMgSweep(swPrimers, bgFasta, { dntp:1.6, mon:50, tRxn:Trxn, mgMin:2, mgMax:10, mgStep:0.5 });

// safety must equal T_rxn − Tm_off exactly in every cell where Tm_off is finite
let safetyOK=true;
sw.rows[0].cells.forEach(c=>{ if (Number.isFinite(c.tmOff) && Math.abs(c.safety-(Trxn-c.tmOff))>1e-9) safetyOK=false; });
assert(safetyOK, 'safety(i,j) == T_rxn − Tm_off in every cell');

// safety is Mg-DEPENDENT: Tm_off rises with Mg ⇒ safety strictly decreases (first vs last col)
const cFirst=sw.rows[0].cells[0], cLast=sw.rows[0].cells[sw.cols.length-1];
assert(cLast.safety < cFirst.safety - 1e-6, 'safety decreases as total Mg2+ rises (Mg-dependent)');

// S (intrinsic discrimination) is ~Mg-invariant: |S_first − S_last| small vs safety swing
assert(Math.abs(cLast.S - cFirst.S) < Math.abs(cFirst.safety - cLast.safety), 'S is far more Mg-invariant than safety');

// suggested column is the LOWEST Mg that is functional AND safe (if feasible at all)
if (sw.feasible){
  const js=sw.jStar;
  assert(sw.colStats[js].feasible, 'suggested column is functional AND safe');
  let firstFeas=-1; for (let j=0;j<sw.colStats.length;j++){ if (sw.colStats[j].feasible){ firstFeas=j; break; } }
  assert(js===firstFeas, 'suggested column is the LOWEST functional&safe Mg column');
}

if (process.exitCode) console.error('\nVALIDATION FAILED'); else console.log('\nALL CHECKS PASSED');
} // end runSelfTest
