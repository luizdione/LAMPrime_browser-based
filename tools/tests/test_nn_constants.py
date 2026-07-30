#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_nn_constants.py - the Watson-Crick nearest-neighbor parameters used by the
LAMPrime engine (app.js / mgspec.mjs / concordance.py) are the SantaLucia (1998)
unified set, checked THREE independent ways:
  1. engine table == SantaLucia 1998 golden (typed from PNAS 95:1460, Table 1);
  2. engine table == BioPython DNA_NN3 (an independent transcription);
  3. the Python (concordance.py) and JS (mgspec.mjs) copies are identical.
"""
from Bio.SeqUtils import MeltingTemp as mt
import golden_literature as G

_ENGINE = G.run_node_dump("_dump_constants.mjs")


def _golden(dn):
    return G.SANTALUCIA_1998.get(dn) or G.SANTALUCIA_1998.get(G.revcomp2(dn))


def test_engine_nn_matches_santalucia1998():
    assert len(_ENGINE["NN_DH"]) == 16
    for dn in _ENGINE["NN_DH"]:
        h, s = _ENGINE["NN_DH"][dn], _ENGINE["NN_DS"][dn]
        gold = _golden(dn)
        assert gold is not None, f"{dn} has no SantaLucia-1998 golden"
        assert abs(h - gold[0]) < 1e-9 and abs(s - gold[1]) < 1e-9, \
            f"{dn}: engine {(h, s)} != SantaLucia1998 {gold}"


def test_biopython_dna_nn3_matches_santalucia1998():
    """Confirms our golden transcription is itself correct: BioPython's
    independent DNA_NN3 agrees with the typed SantaLucia-1998 values."""
    for dn in _ENGINE["NN_DH"]:
        gold = _golden(dn)
        bio = G.biopython_nn_lookup(mt.DNA_NN3, dn)
        assert bio is not None and abs(bio[0] - gold[0]) < 1e-9 and abs(bio[1] - gold[1]) < 1e-9, \
            f"{dn}: BioPython {bio} != golden {gold}"


def test_initiation_constants():
    assert mt.DNA_NN3["init_G/C"] == G.SL98_INIT_GC
    assert mt.DNA_NN3["init_A/T"] == G.SL98_INIT_AT


def test_python_js_nn_parity():
    """The two hard-coded mirrors must not drift apart."""
    for var in ("NN_DH", "NN_DS"):
        py = G.parse_num_map(G.TOOLS + "/concordance.py", var)
        js = G.parse_num_map(G.TOOLS + "/mgspec.mjs", var)
        eng = _ENGINE[var]
        assert py == js, f"{var}: concordance.py != mgspec.mjs"
        assert py == eng, f"{var}: source literal != runtime engine value"
