#!/usr/bin/env node
/* pfalciparum_design.js — desenho LAMP do alvo AT-rico Plasmodium falciparum 18S (#9),
 * usando o MESMO motor do app (importado de crossreact.js), com os parametros default.
 * Demonstra o comportamento do desenhador no regime de baixo GC (estilo Tabela 1 do manuscrito).
 * Uso:  node tools/pfalciparum_design.js
 */
'use strict';
const path = require('path');
const cr = require('./crossreact.js');

const tgt = cr.loadFasta(path.join(cr.DATA, 'pfalciparum_18s_M19172.fasta'));
if (!tgt.seq) { console.error('FASTA ausente: pfalciparum_18s_M19172.fasta'); process.exit(1); }

const P = cr.defaultParams();
const t0 = Date.now();
const des = cr.designLAMP(tgt.seq, P);
const dt = ((Date.now() - t0) / 1000).toFixed(1);

console.log('P. falciparum 18S (M19172.1) — desenho AT-rico (#9); motor = app.js via crossreact.js');
console.log(`alvo: ${tgt.seq.length} nt | GC ${cr.gcPct(tgt.seq).toFixed(1)}% | parametros default (janelas 64-66 / 59-61 C)`);

if (des.errorKey || !des.sets || !des.sets.length) {
  console.log(`RESULTADO: nenhum conjunto valido sob os parametros default (${des.errorKey || 'sem conjuntos'}).`);
  console.log(`candidatos avaliados: ${des.totalCandidates || 0} | tempo ${dt} s`);
  console.log('Interpretacao: alvo AT-rico extremo estressa as janelas de Tm do LAMP — informativo por si so.');
  process.exit(0);
}

const s = des.sets[0];
console.log(`\nconjuntos candidatos: ${des.totalCandidates} | validos: ${des.sets.length} | penalidade do topo: ${s.penalty.toFixed(1)}`);
console.log(`F2-B2: ${s.ampliconF2B2} nt | F3-B3: ${s.ampliconF3B3} nt | tempo: ${dt} s`);
console.log('\nSet #1 (top-ranked):');
console.log('  prim  len   Tm     GC');
const rows = [['F3', s.F3], ['F2', s.F2], ['F1c', s.F1c], ['B1c', s.B1c], ['B2', s.B2], ['B3', s.B3]];
if (s.LF) rows.push(['LF', s.LF]);
if (s.LB) rows.push(['LB', s.LB]);
for (const [nm, c] of rows) {
  console.log(`  ${nm.padEnd(4)} ${String(c.seq.length).padStart(3)}  ${c.tm.toFixed(1).padStart(5)}  ${c.gc.toFixed(0).padStart(3)}`);
}
console.log(`  FIP ${s.FIP.len} nt | BIP ${s.BIP.len} nt`);
if (s.warns && s.warns.length) {
  console.log('  avisos:', s.warns.map(w => `${w.t}:${w.n}${w.v!==undefined?('='+w.v):''}`).slice(0, 12).join(', '));
}
