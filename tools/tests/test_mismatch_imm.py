#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_mismatch_imm.py - the internal single-mismatch NN parameters (MM_IMM in
app.js / mgspec.mjs) come from Allawi & SantaLucia (1997, G.T) and Peyret et al.
(1999, A.A/C.C/G.G/T.T). Each engine entry is checked against BioPython's
DNA_IMM1 table, which is an INDEPENDENT transcription of the same primary
literature. Also verifies all four self-mismatch families are represented.
"""
from Bio.SeqUtils import MeltingTemp as mt
import golden_literature as G

_ENGINE = G.run_node_dump("_dump_constants.mjs")
_IMM = mt.DNA_IMM1


def _bio_imm(key):
    """DNA_IMM1 value for an internal-mismatch key 'ab/wz', honouring the
    engine's own symmetric fallback (z w / b a reversed)."""
    if key in _IMM:
        return _IMM[key]
    top, bot = key.split("/")
    return _IMM.get(bot[::-1] + "/" + top[::-1])


def test_engine_mmimm_matches_biopython_imm1():
    assert len(_ENGINE["MM_IMM"]) >= 48, "unexpectedly small mismatch table"
    for key, (h, s) in _ENGINE["MM_IMM"].items():
        bio = _bio_imm(key)
        assert bio is not None, f"{key} absent from BioPython DNA_IMM1"
        assert abs(h - bio[0]) < 1e-9 and abs(s - bio[1]) < 1e-9, \
            f"{key}: engine {(h, s)} != DNA_IMM1 {bio}"


def test_all_mismatch_families_present():
    """Peyret 1999 (A.A,C.C,G.G,T.T) + Allawi 1997 (G.T) families all covered."""
    families = {"AA", "CC", "GG", "TT", "GT"}
    seen = set()
    for key in _ENGINE["MM_IMM"]:
        top, bot = key.split("/")
        # the mismatched pair is (top[1], bot[0]) at the central position
        pair = "".join(sorted((top[1], bot[0])))
        seen.add(pair)
    canon = {"".join(sorted(f)) for f in families}
    assert canon.issubset(seen), f"missing mismatch families: {canon - seen}"
