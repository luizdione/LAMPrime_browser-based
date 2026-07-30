#!/usr/bin/env node
/* tm_emit.mjs — emit LAMPrime (mgspec.mjs) duplex Tm for a panel of oligos and
 * Mg conditions, as JSON on stdout, so validate_tm.py can compare against
 * BioPython's MeltingTemp (the reference Owczarzy-2008 + Allawi/Peyret impl.).
 *
 * Reads a JSON spec from stdin:
 *   { mon, dntp, oligoNM,
 *     perfect: [{name, seq}],                    // perfect WC duplexes
 *     mismatch:[{name, top, bottom_3to5}],       // explicit duplex (bottom aligned 3'->5')
 *     mg: [2,4,6,8,10] }
 * Emits: { perfect:{name:{mg:tm}}, mismatch:{name:{mg:tm}} }
 *
 * Uses the SAME tmDuplexMg as app.js (imported from mgspec.mjs). No fabrication.
 */
import { tmDuplexMg, wcBottom } from './mgspec.mjs';

function readStdin(){ return new Promise(res=>{ let d=''; process.stdin.on('data',c=>d+=c); process.stdin.on('end',()=>res(d)); }); }

const spec = JSON.parse(await readStdin());
const out = { perfect:{}, mismatch:{} };
for (const o of (spec.perfect||[])){
  out.perfect[o.name] = {};
  const bottom = wcBottom(o.seq);
  for (const mg of spec.mg){
    out.perfect[o.name][mg] = tmDuplexMg(o.seq, bottom, { mg, dntp:spec.dntp, mon:spec.mon, oligoNM:spec.oligoNM }).tm;
  }
}
for (const o of (spec.mismatch||[])){
  out.mismatch[o.name] = {};
  for (const mg of spec.mg){
    out.mismatch[o.name][mg] = tmDuplexMg(o.top, o.bottom_3to5, { mg, dntp:spec.dntp, mon:spec.mon, oligoNM:spec.oligoNM }).tm;
  }
}
process.stdout.write(JSON.stringify(out));
