#!/usr/bin/env node
/* mgtunable_search.mjs — computational screen for a REAL Mg2+-TUNABLE specificity
 * case for the LAMPrime paper. NO fabricated sequences or numbers: every Tm comes
 * from the SAME math as app.js (imported from ./mgspec.mjs) run against REAL
 * NCBI-fetched genomes in ./data/.
 *
 * A "Mg-tunable" case = a published LAMP/RT-LAMP primer whose BEST off-target
 * window in a CLOSE relative has FEW mismatches (ideally 1-3) such that the
 * off-target Tm_off CROSSES the LAMP reaction temperature (~63 C) within the
 * 2-10 mM total MgSO4 sweep: Tm_off(10) >= ~62 (would prime at high Mg) but
 * Tm_off(2) <= ~63 (specific at low Mg), while Tm_on stays functional.
 *
 * For each primer x each off-target genome:
 *   bestOffWindow -> (#mismatch, aligned window)
 *   tmDuplexMg of that window at Mg in {2,4,6,8,10} mM (mono 50, dNTP 1.6, oligo 50 nM)
 *   Tm_on = primer vs its perfect WC complement at the same conditions.
 * FLAG a primer where Tm_off(10) >= 62 and Tm_off(2) <= 63 (a crossing), or more
 * generally where Tm_off in ~[58,68] at 8 mM with 1-4 mismatches.
 *
 * Usage: node tools/mgtunable_search.mjs
 * Refs (engine): SantaLucia 1998; Allawi&SantaLucia 1997; Peyret 1999; Owczarzy 2008.
 */
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { tmDuplexMg, bestOffWindow, wcBottom, revComp, comp, freeMgM } from './mgspec.mjs';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const DATA = path.join(__dirname, 'data');

// ---- conditions (identical to dengue_mgspec.mjs / app defaults) ----
const MG_SWEEP = [2, 4, 6, 8, 10];
const MON = 50;     // monovalent (mM)
const DNTP = 1.6;   // dNTP total (mM)
const OLIGO_NM = 50;
const TREACT = 63;  // LAMP reaction temperature (C)

function loadFasta(fname){
  const raw = fs.readFileSync(path.join(DATA, fname), 'utf-8');
  const lines = raw.split(/\r?\n/);
  const header = (lines.find(l => l.startsWith('>')) || '').slice(1).trim();
  const seq = lines.filter(l => !l.startsWith('>')).join('')
    .toUpperCase().replace(/U/g,'T').replace(/[^ATGC]/g,'');
  return { header, seq };
}

// FIP/BIP: only the 3'-core (F2/B2) anneals & primes off the genomic target.
// F2 = longest 3'-suffix of FIP found (either strand) in the on-target genome.
// B2 = longest 3'-suffix of BIP found (either strand) in the on-target genome.
function locate(primer, genome){
  if (genome.indexOf(primer) !== -1) return true;
  if (genome.indexOf(revComp(primer)) !== -1) return true;
  return false;
}
function core3(seq, genome, lo=14, hi=28){
  for (let L=hi; L>=lo; L--){ const s = seq.slice(-L); if (locate(s, genome)) return s; }
  return null; // not found exactly -> caller falls back to whole sequence
}

// Expand a single degenerate base R={A,G} to its WORSE case (higher Tm_off).
// We just return both A and G variants; caller screens both and keeps the worse.
function expandR(seq){
  if (!seq.includes('R')) return [seq];
  const i = seq.indexOf('R');
  return [seq.slice(0,i)+'A'+seq.slice(i+1), seq.slice(0,i)+'G'+seq.slice(i+1)];
}

// ---- on-target genomes (for FIP/BIP core derivation & on-target sanity) ----
const G = {
  amarginale: loadFasta('amarginale_msp1b_synthetic.fasta'),   // Giglioti 2019 target
  mtb:        loadFasta('mtb_is6110_X17348.fasta'),            // Bentaleb 2016 target
  sars2spike: loadFasta('sarscov2_spike_NC045512.2_21563-25384.fasta'), // Prakash 2023
  denv1:      loadFasta('denv1_EU848545.fasta'),
  zika:       loadFasta('zika_NC035889.fasta'),
};

// ---- off-target panels (close relatives) ----
const OFF = {
  amarginale: [
    { name:'A. centrale CP001759.1',        ...loadFasta('acentrale_genome_CP001759.fasta') },
    { name:'A. ovis Haibei CP015994.1',     ...loadFasta('aovis_genome_CP015994.fasta') },
    { name:'A. phagocytophilum CM177739.1', ...loadFasta('aphago_genome_CM177739.fasta') },
  ],
  mtb: [
    { name:'M. bovis AF2122/97 NC_002945.4', ...loadFasta('mbovis_genome_NC002945.fasta') },
  ],
  sars2spike: [
    { name:'Bat CoV RaTG13 MN996532.2', ...loadFasta('ratg13_MN996532.fasta') },
    { name:'SARS-CoV Tor2 NC_004718.3', ...loadFasta('sarscov_NC004718.fasta') },
  ],
  denv1: [
    { name:'DENV-2 AF038403.1', ...loadFasta('denv2_AF038403.fasta') },
    { name:'DENV-3 M93130.1',   ...loadFasta('denv3_M93130.fasta') },
    { name:'DENV-4 AY947539.1', ...loadFasta('denv4_AY947539.fasta') },
  ],
  zika: [
    { name:'Spondweni DQ859064.1', ...loadFasta('spondweni_DQ859064.fasta') },
    { name:'DENV-2 AF038403.1',    ...loadFasta('denv2_AF038403.fasta') },
  ],
};

// ---- primer sets (verbatim from tools/concordance.py + tools/dengue_mgspec.mjs + brief) ----
// amarginale_msp1b (Giglioti 2019), mtb_is6110 (Bentaleb 2016), sarscov2 spike (Prakash 2023).
const SETS = [];

SETS.push({
  setName:'A. marginale msp1b (Giglioti et al. 2019, Exp Appl Acarol)',
  onKey:'amarginale',
  raw:{ F3:'GCACTACCGTTCATGGATGA', B3:'TCCCCTGTGATATCTGTGCC',
        FIP:'TGCCTTGCCAAATTCTTGCTCCCACCTGACACTGGTGAGAAG',
        BIP:'AGCAGGCTTCAAGCGTACAGTTCCGCGAGCATGTGCA',
        LF:'TCACCCGCTGGTACTTCAA', LB:'GCCTGGAGATGTTAGACCGA' },
});

SETS.push({
  setName:'M. tuberculosis IS6110 (Bentaleb et al. 2016, BMC Infect Dis 16:517)',
  onKey:'mtb',
  raw:{ F3:'TCTCGTCCAGCGCCGCTT', B3:'GCGGGTCCAGATGGCTTG',
        // explicit F2/B2 cores from concordance.py (preferred for FIP/BIP):
        F2:'CCAGCACCTAACCGGCTG', B2:'TCGCGTCGAGGACCATGG',
        FIP:'ACGTAGGCGAACCCTGCCCCCAGCACCTAACCGGCTG',
        BIP:'GTCACCGACGCCTACGCTCTCGCGTCGAGGACCATGG',
        LF:'TCGACACATAGGTGAGGTC', LB:'TCGCTTCCACGATGGCCA' },
});

SETS.push({
  setName:'SARS-CoV-2 spike (Prakash et al. 2023, MethodsX)',
  onKey:'sars2spike',
  raw:{ F3:'TGGTGATATTGCTGCTAGA', B3:'GCACTATTAAATTGGTGGGC',
        FIP:'AGGTCCAACCAGAAGTGATTCACCTTTGCTCACAGATG',
        BIP:'GCAGGTGCTGCATTACAATCTGTGTAACTCCAATACCA',
        LF:'GCTAACAGTGCAGAAGTGTATT', LB:'GCTATGCAAATGGCTTATAGGT' },
});

// Dengue DENV-1 (Hu S-F et al. 2015, BMC Microbiology 15:265) — serotype-specific.
SETS.push({
  setName:'DENV-1 (Hu et al. 2015, BMC Microbiology 15:265)',
  onKey:'denv1',
  raw:{ F3:'TGTGTTCCTCCTTCTCATAATG', B3:'CAGACTCAATCCAATCGTAAGA',
        FIP:'CATCCTGTCTGAAGCATTGGCTGGACAATTGACATGGAATGATC',
        BIP:'CCTAGCTCTGATGGCCACTTTCTTCTCTAGATGTTAGTCTGCG',
        LF:'CCAACCATGATGCATAACCTG', LB:'ATGAGACCAATGTTCGCTGT' },
});

// Zika (da Silva et al. 2019) — FIP has a degenerate R (expand to A/G, worse case).
SETS.push({
  setName:'Zika (da Silva et al. 2019)',
  onKey:'zika',
  raw:{ F3:'CAGTTCACACGGCCCTTG', B3:'TGTACCTCCACTGTGACTGT',
        FIP:'GGCGACATTTCAAGTGGCCAGAGAGCTCTRGAGGCTGAGA',
        BIP:'AGGGCGTGTCATACTCCTTGTGAGTGTTTCAGCCGGGATCT',
        LF:'CCTTCCCTTTGCACCATCCA', LB:'TACCGCAGCGTTCACATTCA' },
});

// For each set, build the list of TARGET-BINDING units to screen:
//  F3,B3,LF,LB  -> use verbatim.
//  FIP          -> F2 core (explicit if given, else longest 3'-suffix in on-target genome).
//  BIP          -> B2 core (explicit if given, else longest 3'-suffix in on-target genome).
function bindingUnits(set){
  const onSeq = G[set.onKey].seq;
  const units = [];
  const r = set.raw;
  const push = (role, seq, note) => { if (seq) units.push({ role, seq, note }); };
  push('F3', r.F3, 'outer fwd');
  push('B3', r.B3, 'outer rev');
  // FIP core
  let f2 = r.F2;
  if (!f2 && r.FIP) f2 = core3(r.FIP, onSeq) || r.FIP;
  push('FIP(F2 core)', f2, 'FIP target-binding core');
  // BIP core
  let b2 = r.B2;
  if (!b2 && r.BIP) b2 = core3(r.BIP, onSeq) || r.BIP;
  push('BIP(B2 core)', b2, 'BIP target-binding core');
  push('LF', r.LF, 'loop fwd');
  push('LB', r.LB, 'loop rev');
  return units;
}

// Tm of a primer vs a given aligned bottom strand at each Mg in the sweep.
function tmSweep(top, bottom){
  const out = {};
  for (const mg of MG_SWEEP){
    out[mg] = tmDuplexMg(top, bottom, { mg, dntp:DNTP, mon:MON, oligoNM:OLIGO_NM }).tm;
  }
  return out;
}

// ---- main screen ----
const allRows = [];   // one row per (set, primer, best-off-window)
console.log('=== Mg-tunable specificity screen (real fetched genomes; app Tm engine) ===');
console.log(`Conditions: monovalent ${MON} mM, dNTP ${DNTP} mM, oligo ${OLIGO_NM} nM; reaction T = ${TREACT} C.`);
console.log(`Mg sweep (total MgSO4): ${MG_SWEEP.join(', ')} mM\n`);

for (const set of SETS){
  const units = bindingUnits(set);
  const offPanel = OFF[set.onKey];
  console.log(`---- ${set.setName} | off-target panel: ${offPanel.map(o=>o.name).join('; ')} ----`);
  for (const u of units){
    // expand degenerate R to worst case
    const variants = expandR(u.seq);
    let best = null; // worst-case (highest Tm_off@8mM) variant result
    for (const v of variants){
      const off = bestOffWindow(v, offPanel);
      if (!off) continue;
      const tmOff = tmSweep(v, off.alignedTargetBottom);
      if (!best || tmOff[8] > best.tmOff[8]) best = { v, off, tmOff };
    }
    if (!best){ console.log(`  ${u.role.padEnd(12)} NO off-target window`); continue; }
    // on-target Tm (perfect WC complement of the chosen variant)
    const tmOn = tmSweep(best.v, wcBottom(best.v));
    // sanity: best window of this unit in its OWN on-target genome (expect ~0 mm)
    const onCtrl = bestOffWindow(best.v, [{ name:set.onKey, seq:G[set.onKey].seq }]);

    const row = {
      set:set.setName, onKey:set.onKey, role:u.role, primer:best.v, len:best.v.length,
      offSrc:best.off.src, offPos:best.off.pos+1, offStrand:best.off.strand, mm:best.off.nMismatch,
      tmOff2:best.tmOff[2], tmOff4:best.tmOff[4], tmOff6:best.tmOff[6], tmOff8:best.tmOff[8], tmOff10:best.tmOff[10],
      tmOn8:tmOn[8], tmOn2:tmOn[2], tmOn10:tmOn[10],
      onCtrlMm: onCtrl ? onCtrl.nMismatch : NaN,
    };
    // crossing flags
    row.crosses63 = (best.tmOff[2] <= TREACT && best.tmOff[10] >= TREACT);
    row.crossingMg = null;
    if (row.crosses63){
      // linear-interp the Mg at which Tm_off == 63 (monotone increasing in Mg)
      const xs = MG_SWEEP, ys = MG_SWEEP.map(m=>best.tmOff[m]);
      for (let i=0;i<xs.length-1;i++){
        if (ys[i] <= TREACT && ys[i+1] >= TREACT){
          const t = (TREACT - ys[i])/(ys[i+1]-ys[i]);
          row.crossingMg = +(xs[i] + t*(xs[i+1]-xs[i])).toFixed(2);
          break;
        }
      }
    }
    row.nearBoundary = (best.tmOff[8] >= 58 && best.tmOff[8] <= 68 && best.off.nMismatch>=1 && best.off.nMismatch<=4);
    row.strongCross = (best.tmOff[10] >= 62 && best.tmOff[2] <= 63);
    allRows.push(row);

    const flag = row.strongCross ? ' <<< STRONG CROSSING' : (row.crosses63 ? ' << crosses 63' : (row.nearBoundary ? ' < near boundary' : ''));
    console.log(`  ${u.role.padEnd(12)} ${best.v.padEnd(28)} mm=${best.off.nMismatch} (on-ctrl mm=${row.onCtrlMm}) | Tm_off 2/8/10=${best.tmOff[2].toFixed(1)}/${best.tmOff[8].toFixed(1)}/${best.tmOff[10].toFixed(1)} | Tm_on8=${tmOn[8].toFixed(1)} | ${best.off.src} pos${best.off.pos+1}${best.off.strand}${flag}`);
  }
  console.log('');
}

// ---- write full CSV ----
const CSV = path.join(DATA, 'mgtunable_search.csv');
const head = 'set,on_target,role,primer,len,n_mismatch,on_ctrl_mm,off_target_src,off_pos,off_strand,tm_off_2,tm_off_4,tm_off_6,tm_off_8,tm_off_10,tm_on_2,tm_on_8,tm_on_10,crosses63,crossing_mg,strong_cross,near_boundary';
const body = allRows.map(r => [
  '"'+r.set+'"', r.onKey, r.role, r.primer, r.len, r.mm, r.onCtrlMm,
  '"'+r.offSrc+'"', r.offPos, r.offStrand,
  r.tmOff2.toFixed(2), r.tmOff4.toFixed(2), r.tmOff6.toFixed(2), r.tmOff8.toFixed(2), r.tmOff10.toFixed(2),
  r.tmOn2.toFixed(2), r.tmOn8.toFixed(2), r.tmOn10.toFixed(2),
  r.crosses63, r.crossingMg===null?'':r.crossingMg, r.strongCross, r.nearBoundary,
].join(',')).join('\n');
fs.writeFileSync(CSV, head+'\n'+body+'\n', 'utf-8');

// ---- ranking ----
// Score: prefer rows that CROSS 63 cleanly with FEWEST mismatches and a functional
// on-target (Tm_on8 >= 59). Sort key: crosses63 first, then strongCross, then
// (smaller mm), then closeness of Tm_off8 to 63 (smaller |Tm_off8-63|).
function score(r){
  const funcOn = r.tmOn8 >= 59 ? 1 : 0;
  const cross = r.crosses63 ? 1 : 0;
  const strong = r.strongCross ? 1 : 0;
  const near = Math.abs(r.tmOff8 - TREACT);
  return { funcOn, cross, strong, mm:r.mm, near, r };
}
const ranked = allRows.map(score).sort((a,b)=>{
  if (b.cross !== a.cross) return b.cross - a.cross;
  if (b.strong !== a.strong) return b.strong - a.strong;
  if (b.funcOn !== a.funcOn) return b.funcOn - a.funcOn;
  if (a.near !== b.near) return a.near - b.near;       // closest to 63 first
  if (a.mm !== b.mm) return a.mm - b.mm;               // fewest mismatches
  return 0;
});

console.log('================ RANKED CANDIDATES (top 12) ================');
console.log('rank set/primer | mm | Tm_off@2/4/6/8/10 | Tm_on@8 | crosses63? @Mg | offtarget');
ranked.slice(0,12).forEach((s,i)=>{
  const r = s.r;
  const setShort = r.set.split('(')[0].trim();
  const cr = r.crosses63 ? `YES @${r.crossingMg} mM` : (r.nearBoundary?'near':'no');
  console.log(`${String(i+1).padStart(2)}. ${setShort} / ${r.role} | mm=${r.mm} | ${r.tmOff2.toFixed(1)}/${r.tmOff4.toFixed(1)}/${r.tmOff6.toFixed(1)}/${r.tmOff8.toFixed(1)}/${r.tmOff10.toFixed(1)} | on8=${r.tmOn8.toFixed(1)} | ${cr} | ${r.offSrc} pos${r.offPos}${r.offStrand}`);
});

// ---- best candidate verdict ----
const crossers = ranked.filter(s => s.cross && s.funcOn);
console.log('\n================ VERDICT ================');
if (crossers.length){
  const b = crossers[0].r;
  console.log('BEST Mg-tunable candidate (crosses 63 C with functional on-target):');
  console.log(`  Set      : ${b.set}`);
  console.log(`  Primer   : ${b.role} = ${b.primer} (${b.len} nt)`);
  console.log(`  Off-target: ${b.offSrc} pos ${b.offPos} (${b.offStrand} strand), ${b.mm} mismatch(es); on-target ctrl mm=${b.onCtrlMm}`);
  console.log(`  Tm_off   : 2mM=${b.tmOff2.toFixed(2)}  4mM=${b.tmOff4.toFixed(2)}  6mM=${b.tmOff6.toFixed(2)}  8mM=${b.tmOff8.toFixed(2)}  10mM=${b.tmOff10.toFixed(2)}`);
  console.log(`  Tm_on    : 2mM=${b.tmOn2.toFixed(2)}  8mM=${b.tmOn8.toFixed(2)}  10mM=${b.tmOn10.toFixed(2)}`);
  console.log(`  Crossing : Tm_off crosses 63 C at ~${b.crossingMg} mM total MgSO4 (specific below, primes above).`);
} else {
  console.log('NO primer has Tm_off crossing 63 C across the 2-10 mM Mg sweep with a functional on-target.');
  const maxOff = Math.max(...allRows.map(r=>r.tmOff10));
  console.log(`Highest Tm_off anywhere (at 10 mM) = ${maxOff.toFixed(2)} C.`);
  const nb = ranked.filter(s=>s.r.nearBoundary);
  if (nb.length){
    console.log('Closest-to-boundary candidates (Tm_off@8 in [58,68], 1-4 mm) but no full crossing:');
    nb.slice(0,5).forEach(s=>{ const r=s.r; console.log(`  ${r.set.split('(')[0].trim()} / ${r.role} mm=${r.mm} Tm_off@8=${r.tmOff8.toFixed(1)} (${r.offSrc})`); });
  }
  console.log('\nRECOMMENDATION: fall back to Plan-A (Dengue robust-specificity figure); do NOT force a weak case.');
}

console.log(`\nWrote full matrix -> ${CSV}`);
console.log('DONE.');
