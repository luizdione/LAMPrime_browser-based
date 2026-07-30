#!/usr/bin/env node
/* _dump_owcz.mjs - evaluate the LAMPrime Owczarzy-2008 correction and the 1:1
 * Mg-dNTP chelation on a grid supplied as JSON on stdin, emitting JSON on
 * stdout. Lets the Python test compare the engine (../mgspec.mjs, mirror of
 * app.js) against an independent re-implementation of Owczarzy 2008.
 *
 * stdin : [{tmK1M, freeMg, monM, fGC, Nbp, totMg, dntp}, ...]
 * stdout: [{...input, tmC, freeMgOut}, ...]
 */
import { owczarzy2008, freeMgM } from '../mgspec.mjs';

function readStdin(){ return new Promise(r=>{ let d=''; process.stdin.on('data',c=>d+=c); process.stdin.on('end',()=>r(d)); }); }

const grid = JSON.parse(await readStdin());
const out = grid.map(p => ({
  ...p,
  tmC: owczarzy2008(p.tmK1M, p.freeMg, p.monM, p.fGC, p.Nbp),
  freeMgOut: freeMgM(p.totMg, p.dntp),
}));
process.stdout.write(JSON.stringify(out));
