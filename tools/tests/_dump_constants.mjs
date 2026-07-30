#!/usr/bin/env node
/* _dump_constants.mjs - print the LAMPrime runtime NN / mismatch constant tables
 * (imported from ../mgspec.mjs, the importable mirror of app.js) as JSON, so the
 * Python tests can check them against the primary literature and BioPython.
 * No fabrication: values come straight from the engine module.
 */
import { NN_DH, NN_DS, MM_IMM } from '../mgspec.mjs';
process.stdout.write(JSON.stringify({ NN_DH, NN_DS, MM_IMM }));
