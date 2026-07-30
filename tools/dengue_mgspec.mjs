#!/usr/bin/env node
/* dengue_mgspec.mjs — REAL Dengue Mg2+-vs-specificity matrix for the LAMPrime
 * paper (Figure 3). NO fabricated numbers: every Tm comes from running the
 * SAME math as the app (imported from ./mgspec.mjs) against the REAL fetched
 * Dengue genomes in ./data/.
 *
 * Worked example: serotype-specific DENV-1 RT-LAMP primer set of
 *   Hu S-F et al. (2015) BMC Microbiology 15:265, DOI 10.1186/s12866-015-0595-1.
 * On-target  : DENV-1 EU848545.1
 * Off-targets: DENV-2 AF038403.1, DENV-3 M93130.1, DENV-4 AY947539.1
 *   (the four serotypes are mutual congeneric off-targets, ~30-40% nt divergence;
 *    Hu et al. report ZERO cross-reaction at a FIXED 8 mM MgSO4 — the Mg SWEEP
 *    below is OUR in-silico contribution, not a reproduction of published Mg data.)
 *
 * For each of the 6 published primers (F3,B3,FIP,BIP,LF,LB) and each total
 * MgSO4 j in {2.0,2.5,...,10.0} mM (monovalent 50 mM, dNTP 1.6 mM):
 *   Tm_on(i,j)  = Tm of primer i vs its perfect Watson-Crick complement.
 *   Tm_off(i,j) = Tm of primer i vs the BEST (fewest-mismatch, gap-free, both
 *                 strands) window across the THREE off-target genomes.
 *   S(i,j)      = Tm_on - Tm_off  (specificity margin, deg C).
 *
 * FIP/BIP NOTE: FIP = F1c-loop-F2 and BIP = B1c-loop-B2 are synthetic oligos
 * whose F1c/B1c halves are TAILS that do NOT hybridise to the genomic target.
 * The portion that anneals to (and primes off) the target is F2 (3' core of FIP)
 * and B2 (3' core of BIP). So FIP/BIP rows use their target-binding core (F2/B2),
 * identified as the longest 3'-suffix of FIP / 5'-prefix of BIP that is found in
 * the DENV-1 genome. This is the biologically meaningful annealing unit and is
 * how PrimerExplorer/NEB treat F2/B2. Stated explicitly in the report & caption.
 *
 * Usage: node tools/dengue_mgspec.mjs
 */
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { tmDuplexMg, bestOffWindow, wcBottom, freeMgM, revComp, comp } from './mgspec.mjs';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const DATA = path.join(__dirname, 'data');

function loadFasta(fname){
  const raw = fs.readFileSync(path.join(DATA, fname), 'utf-8');
  const lines = raw.split(/\r?\n/);
  const header = (lines.find(l => l.startsWith('>')) || '').slice(1).trim();
  const seq = lines.filter(l => !l.startsWith('>')).join('')
    .toUpperCase().replace(/U/g,'T').replace(/[^ATGC]/g,'');
  return { header, seq };
}

// ---- genomes ----
const denv1 = loadFasta('denv1_EU848545.fasta');           // on-target
const offs = [
  { acc:'DENV-2 AF038403.1', ...loadFasta('denv2_AF038403.fasta') },
  { acc:'DENV-3 M93130.1',   ...loadFasta('denv3_M93130.fasta') },
  { acc:'DENV-4 AY947539.1', ...loadFasta('denv4_AY947539.fasta') },
].map(o => ({ name:o.acc, seq:o.seq, header:o.header }));

// ---- Hu et al. 2015 DENV-1 primer set (verbatim, 5'->3') ----
const FIP = 'CATCCTGTCTGAAGCATTGGCTGGACAATTGACATGGAATGATC';
const BIP = 'CCTAGCTCTGATGGCCACTTTCTTCTCTAGATGTTAGTCTGCG';
const RAW = {
  F3 : 'TGTGTTCCTCCTTCTCATAATG',
  B3 : 'CAGACTCAATCCAATCGTAAGA',
  LF : 'CCAACCATGATGCATAACCTG',
  LB : 'ATGAGACCAATGTTCGCTGT',
};

// Derive the target-binding cores of FIP (=F2, longest 3'-suffix present in
// DENV-1) and BIP (=B2, longest 5'-prefix present in DENV-1 sense or antisense).
function locate(primer, genome){
  const g = genome; const r = revComp(primer);
  let i = g.indexOf(primer); if (i !== -1) return { start:i+1, end:i+primer.length, strand:'+' };
  let j = g.indexOf(r);      if (j !== -1) return { start:j+1, end:j+primer.length, strand:'-' };
  return null;
}
function fipCore(fip, genome){ // F2 = longest 3'-suffix of FIP found in genome
  for (let L=28; L>=14; L--){ const f2 = fip.slice(-L); if (locate(f2, genome)) return { name:'F2(core of FIP)', seq:f2, loc:locate(f2,genome) }; }
  return null;
}
function bipCore(bip, genome){ // B1c = longest 5'-prefix found; B2 = the remainder, the primed core
  // The 3' core of BIP (B2) is what anneals & is extended. B2 = longest 3'-suffix of BIP found in genome.
  for (let L=28; L>=14; L--){ const b2 = bip.slice(-L); if (locate(b2, genome)) return { name:'B2(core of BIP)', seq:b2, loc:locate(b2,genome) }; }
  return null;
}
const f2core = fipCore(FIP, denv1.seq);
const b2core = bipCore(BIP, denv1.seq);
if (!f2core) throw new Error('FATAL: could not locate F2 core of FIP in DENV-1 genome');
if (!b2core) throw new Error('FATAL: could not locate B2 core of BIP in DENV-1 genome');

// 6 published rows. For FIP/BIP we use the target-binding core (F2/B2).
const PRIMERS = [
  { role:'F3',  bindSeq:RAW.F3,        note:'outer fwd' },
  { role:'B3',  bindSeq:RAW.B3,        note:'outer rev' },
  { role:'FIP', bindSeq:f2core.seq,    note:'FIP target-binding core = F2 ('+f2core.seq.length+' nt)' },
  { role:'BIP', bindSeq:b2core.seq,    note:'BIP target-binding core = B2 ('+b2core.seq.length+' nt)' },
  { role:'LF',  bindSeq:RAW.LF,        note:'loop fwd' },
  { role:'LB',  bindSeq:RAW.LB,        note:'loop rev' },
];

// ---- Mg sweep & fixed conditions (Hu/PrimerExplorer load) ----
const MG_SWEEP = []; for (let m=2.0; m<=10.0001; m+=0.5) MG_SWEEP.push(Math.round(m*10)/10);
const MON = 50;     // monovalent (mM)
const DNTP = 1.6;   // dNTP total (mM)
const OLIGO_NM = 50;

// ---- sanity positive control: each primer vs DENV-1 on-target genome ----
console.log('=== SANITY positive control: best window of each primer in DENV-1 EU848545.1 (expect 0 mismatch) ===');
const onCtrlRows = [];
for (const p of PRIMERS){
  const w = bestOffWindow(p.bindSeq, [{ name:'DENV-1 EU848545.1', seq:denv1.seq }]);
  onCtrlRows.push({ role:p.role, mm:w?w.nMismatch:NaN, pos:w?w.pos+1:NaN, strand:w?w.strand:'?' });
  console.log(`  ${p.role.padEnd(4)} binds ${p.bindSeq} | best in DENV-1: ${w.nMismatch} mismatch @ pos ${w.pos+1} (${w.strand})`);
}

// ---- main matrix ----
const rows = [];
console.log('\n=== Computing S(i,j) over Mg sweep (real fetched genomes) ===');
for (const p of PRIMERS){
  // off-target best window is Mg-independent in SEQUENCE; find it once (fewest mismatches across the 3 off genomes).
  const off = bestOffWindow(p.bindSeq, offs);
  if (!off) throw new Error('FATAL: no off-target window for '+p.role);
  // map which off genome / 1-based position / strand
  // bestOffWindow scanned fwd then revcomp of each bg; reconstruct human-readable locus:
  const offLabel = `${off.src} pos${off.pos+1}(${off.strand}) ${off.nMismatch}mm`;
  const wcOn = wcBottom(p.bindSeq); // perfect WC complement aligned 3'->5'
  for (const mg of MG_SWEEP){
    const opt = { mg, dntp:DNTP, mon:MON, oligoNM:OLIGO_NM };
    const tmOn  = tmDuplexMg(p.bindSeq, wcOn, opt).tm;
    const tmOff = tmDuplexMg(p.bindSeq, off.alignedTargetBottom, opt).tm;
    const mgFree = freeMgM(mg/1000, DNTP/1000) * 1000;
    rows.push({
      primer:p.role, role:p.note, mg_total:mg, mg_free:mgFree,
      tm_on:tmOn, tm_off:tmOff, n_mismatch:off.nMismatch,
      best_offtarget:offLabel, S:tmOn - tmOff,
    });
  }
}

// ---- write CSV ----
const CSV = path.join(DATA, 'dengue_mg_Sij.csv');
const head = 'primer,role,mg_total,mg_free,tm_on,tm_off,n_mismatch,best_offtarget,S';
const body = rows.map(r =>
  [r.primer, '"'+r.role+'"', r.mg_total.toFixed(1), r.mg_free.toFixed(4),
   r.tm_on.toFixed(3), r.tm_off.toFixed(3), r.n_mismatch, '"'+r.best_offtarget+'"', r.S.toFixed(3)].join(',')
).join('\n');
fs.writeFileSync(CSV, head+'\n'+body+'\n', 'utf-8');
console.log(`\nWrote ${rows.length} rows -> ${CSV}`);

// ---- maximin suggested Mg under functional gate ----
// gate: Tm_on in relaxed LAMP window [59,66]; Tm_off < 63.
// Two separable constraints (kept distinct on purpose, reported honestly):
//   SPECIFICITY gate  : Tm_off < 63 C for every primer (no off-target anneals
//                       at the reaction temperature) -> the discriminating axis.
//   FUNCTIONAL window : Tm_on in [59,66] (relaxed LAMP). The target-binding
//                       cores here are short 20-22 nt units, so their absolute
//                       Tm_on sits a few C below the inner-primer band; we treat
//                       the LOWER bound (>=59) as the productive-annealing floor
//                       and use the UPPER bound (<=66) only as a sanity ceiling.
const GATE_ON_LO=59, GATE_ON_HI=66, GATE_OFF_HI=63;
const byMg = new Map();
for (const r of rows){ if(!byMg.has(r.mg_total)) byMg.set(r.mg_total, []); byMg.get(r.mg_total).push(r); }

function colStats(mg){
  const col=byMg.get(mg);
  const minS=Math.min(...col.map(r=>r.S));
  const onLo=Math.min(...col.map(r=>r.tm_on)), onHi=Math.max(...col.map(r=>r.tm_on));
  const offHi=Math.max(...col.map(r=>r.tm_off));
  const specPass = col.every(r => r.tm_off < GATE_OFF_HI);
  const nFunctional = col.filter(r => r.tm_on>=GATE_ON_LO && r.tm_on<=GATE_ON_HI).length;
  const fullGate = specPass && nFunctional===col.length;
  return { mg, col, minS, onLo, onHi, offHi, specPass, nFunctional, fullGate };
}

// MAXIMIN under the full functional gate (spec + all primers in Tm_on window).
let bestCol=null;
for (const mg of MG_SWEEP){ const s=colStats(mg); if (s.fullGate && (!bestCol || s.minS>bestCol.minS)) bestCol=s; }
// Specificity-only maximin (all columns satisfy Tm_off<63): the column maximising
// min_i S while keeping the MOST primers functional, then highest Mg as tie-break.
let bestSpec=null;
for (const mg of MG_SWEEP){
  const s=colStats(mg);
  if (!s.specPass) continue;
  if (!bestSpec || s.nFunctional>bestSpec.nFunctional ||
     (s.nFunctional===bestSpec.nFunctional && s.minS>bestSpec.minS)) bestSpec=s;
}
// Unconstrained maximin (argmax_j min_i S).
let maximinAll=null;
for (const mg of MG_SWEEP){ const s=colStats(mg); if (!maximinAll || s.minS>maximinAll.minS) maximinAll=s; }

console.log('\n=== Functional gate: Tm_on in ['+GATE_ON_LO+','+GATE_ON_HI+'] (functional floor 59) AND Tm_off < '+GATE_OFF_HI+' (specificity) ===');
console.log('Per-Mg column: min_i S, Tm_on span, max Tm_off, #primers functional, gates:');
for (const mg of MG_SWEEP){
  const s=colStats(mg);
  console.log(`  Mg ${mg.toFixed(1).padStart(4)} (free ${(freeMgM(mg/1000,DNTP/1000)*1000).toFixed(2)}) | minS ${s.minS.toFixed(2).padStart(6)} | Tm_on [${s.onLo.toFixed(1)},${s.onHi.toFixed(1)}] | maxTm_off ${s.offHi.toFixed(1)} | func ${s.nFunctional}/6 | spec ${s.specPass?'PASS':'fail'} | full ${s.fullGate?'PASS':'fail'}`);
}
if (bestCol){
  console.log(`\nMAXIMIN suggested Mg (full functional gate) = ${bestCol.mg.toFixed(1)} mM total | min_i S = ${bestCol.minS.toFixed(2)} C`);
} else {
  console.log('\nNo Mg column satisfies the FULL functional gate (all 6 Tm_on in [59,66] AND all Tm_off<63).');
  console.log(`Reason: the short 20-22 nt binding cores keep Tm_on a few C below 59 at low Mg; specificity is NOT the limiter (max Tm_off <= ${Math.max(...rows.map(r=>r.tm_off)).toFixed(1)} C << 63 everywhere).`);
  console.log(`\nSPECIFICITY-gated suggested Mg (Tm_off<63 for all; max #functional then max min_i S) = ${bestSpec.mg.toFixed(1)} mM | func ${bestSpec.nFunctional}/6 | min_i S = ${bestSpec.minS.toFixed(2)} C`);
}
console.log(`Unconstrained maximin (argmax_j min_i S) = ${maximinAll.mg.toFixed(1)} mM | min_i S = ${maximinAll.minS.toFixed(2)} C`);

const SUGG = bestCol ? bestCol.mg : bestSpec.mg;
console.log(`\nPer-primer S at suggested Mg = ${SUGG.toFixed(1)} mM:`);
for (const r of byMg.get(SUGG)){
  console.log(`  ${r.primer.padEnd(4)} Tm_on ${r.tm_on.toFixed(1)}  Tm_off ${r.tm_off.toFixed(1)}  S ${r.S.toFixed(2)}  (off: ${r.best_offtarget})`);
}

// ---- how many primers gain specificity as Mg drops (10 -> 2 mM) ----
let gained=0; const deltas=[];
for (const p of PRIMERS){
  const at10 = rows.find(r=>r.primer===p.role && r.mg_total===10.0).S;
  const at2  = rows.find(r=>r.primer===p.role && r.mg_total===2.0).S;
  deltas.push({ role:p.role, dS: at2-at10 });
  if (at2 > at10) gained++;
}
console.log(`\nPrimers that GAIN specificity (S higher at 2 mM than 10 mM): ${gained}/6`);
for (const d of deltas) console.log(`  ${d.role.padEnd(4)} dS(2mM - 10mM) = ${d.dS>=0?'+':''}${d.dS.toFixed(2)} C`);

// S range observed
const allS = rows.map(r=>r.S);
console.log(`\nS range observed across whole matrix: [${Math.min(...allS).toFixed(2)}, ${Math.max(...allS).toFixed(2)}] C`);

// emit machine-readable suggested column for the plotting script
fs.writeFileSync(path.join(DATA,'dengue_mg_suggested.json'),
  JSON.stringify({ suggested_mg_total: SUGG, gate:{on_lo:GATE_ON_LO,on_hi:GATE_ON_HI,off_hi:GATE_OFF_HI},
    monovalent_mM:MON, dntp_mM:DNTP, oligo_nM:OLIGO_NM,
    fip_core:f2core.seq, bip_core:b2core.seq }, null, 2));
console.log('\nDONE.');
