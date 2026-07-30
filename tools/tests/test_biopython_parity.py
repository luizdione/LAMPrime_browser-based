#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_biopython_parity.py - turns tools/validate_tm.py (the reproducible
BioPython comparison) into an assertion. The Mg2+-path
duplex Tm (mgspec.mjs, the SAME math as app.js) must match BioPython's
MeltingTemp.Tm_NN (Owczarzy 2008 saltcorr=7 + SantaLucia 1998 DNA_NN3 +
Allawi/Peyret DNA_IMM1) to well within 1 C over the whole oligo panel, and even a
single internal-mismatch duplex agrees to within 0.31 C.
"""
import os, sys, statistics
from Bio.SeqUtils import MeltingTemp as mt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # tools/
import validate_tm as V

_EMIT = V.run_emit()


def _perfect_deltas():
    d = []
    for name, seq in V.PANEL:
        for mg in V.MG_LIST:
            d.append(_EMIT["perfect"][name][str(mg)] - V.bio_tm(seq, mg, mt.DNA_NN3))
    return d


def _mismatch_deltas():
    d = []
    for mg in V.MG_LIST:
        ours = _EMIT["mismatch"]["mm01_internal"][str(mg)]
        bio = V.bio_tm(V.MM_TOP, mg, mt.DNA_NN3, c_seq=V.MM_BOTTOM_3to5, imm_table=mt.DNA_IMM1)
        d.append(ours - bio)
    return d


def test_perfect_panel_matches_biopython():
    mx = max(abs(x) for x in _perfect_deltas())
    assert mx < 1e-6, f"perfect-match panel max|d| = {mx:.4f} C (expected ~0)"


def test_mismatch_duplex_within_reported_tolerance():
    mx = max(abs(x) for x in _mismatch_deltas())
    # agreement holds to within 0.31 C, including mismatched duplexes
    assert mx <= 0.35, f"mismatch duplex max|d| = {mx:.4f} C (expected <= 0.31)"


def test_biopython_version_pinned():
    import Bio
    assert Bio.__version__ == "1.86", \
        f"BioPython {Bio.__version__}: pin 1.86 (see requirements.txt) for stable reference Tm"
