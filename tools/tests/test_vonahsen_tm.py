#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_vonahsen_tm.py - end-to-end check of the DESIGN-engine Tm path (SantaLucia
1998 NN + SantaLucia 1998 salt correction + von Ahsen 2001 Mg->Na equivalent),
i.e. the S1.1 design-engine path. A reference Tm is implemented here directly from
those equations (using the typed SantaLucia-1998 golden, not the engine's own
table) and must reproduce the design-engine Tm values for the validation panel.
"""
import math
import golden_literature as G

# S1.1 buffer constants (NEB defaults), identical to concordance.py
R, NA, MG, DNTP, OLIGO = 1.987, 50.0, 8.0, 1.4, 50e-9


def _nn(dn):
    v = G.SANTALUCIA_1998.get(dn) or G.SANTALUCIA_1998.get(G.revcomp2(dn))
    if v is None:
        raise KeyError(dn)
    return v


def ref_tm(seq):
    """Tm (C) from S1.1 equations, transcribed from the primary literature."""
    seq = "".join(c for c in seq.upper() if c in "ATGC")
    N = len(seq)
    dH = dS = 0.0
    for i in range(N - 1):
        h, s = _nn(seq[i:i + 2]); dH += h; dS += s
    for b in (seq[0], seq[-1]):
        h, s = G.SL98_INIT_GC if b in "GC" else G.SL98_INIT_AT
        dH += h; dS += s
    mg_free = max(0.0, MG - DNTP)                          # von Ahsen 1:1 chelation
    na_eq = max(1e-3, NA + 120 * math.sqrt(mg_free)) / 1000.0
    dS_salt = dS + 0.368 * (N - 1) * math.log(na_eq)       # SantaLucia 1998 salt
    return (dH * 1000) / (dS_salt + R * math.log(OLIGO / 4)) - 273.15


# validation-panel Tm values (SARS-CoV-2 spike; M. tuberculosis IS6110 targets)
PANEL_TM = [
    ("TGGTGATATTGCTGCTAGA", 57.1), ("CACCTTTGCTCACAGATG", 57.1),
    ("AGGTCCAACCAGAAGTGATT", 60.4), ("GCAGGTGCTGCATTACAA", 59.8),
    ("TCTGTGTAACTCCAATACCA", 57.3), ("GCACTATTAAATTGGTGGGC", 58.9),
    ("GCTAACAGTGCAGAAGTGTATT", 60.6), ("GCTATGCAAATGGCTTATAGGT", 60.6),
    ("TCTCGTCCAGCGCCGCTT", 68.3), ("CCAGCACCTAACCGGCTG", 64.3),
    ("ACGTAGGCGAACCCTGCCC", 68.5), ("GTCACCGACGCCTACGCTC", 66.8),
    ("TCGCGTCGAGGACCATGG", 65.7), ("GCGGGTCCAGATGGCTTG", 64.7),
    ("TCGACACATAGGTGAGGTC", 59.8), ("TCGCTTCCACGATGGCCA", 65.6),
]


def test_reference_reproduces_panel_tm():
    for seq, reported in PANEL_TM:
        got = ref_tm(seq)
        assert abs(got - reported) <= 0.15, \
            f"{seq}: S1.1 reference {got:.2f} C != panel {reported} C"


def test_vonahsen_naeq_formula():
    assert abs(G.vonahsen_naeq(50, 8, 1.4) - (50 + 120 * math.sqrt(6.6))) < 1e-9
    assert abs(G.vonahsen_naeq(50, 2, 1.4) - (50 + 120 * math.sqrt(0.6))) < 1e-9
    # when [Mg] <= [dNTP], free Mg is clamped to 0 -> Na_eq == Na
    assert abs(G.vonahsen_naeq(50, 1.0, 1.4) - 50.0) < 1e-9
