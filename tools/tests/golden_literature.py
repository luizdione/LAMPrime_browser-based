#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
golden_literature.py - GOLDEN thermodynamic constants transcribed DIRECTLY from
the primary literature, plus small helpers. These values are the independent
reference against which the LAMPrime engine (app.js / tools/mgspec.mjs /
tools/concordance.py) is checked. They are typed from the papers, NOT copied
from our own source code, so the tests are not tautological.

Primary sources
---------------
SantaLucia (1998) PNAS 95:1460-1465, Table 1 (unified NN parameters).
Allawi & SantaLucia (1997) Biochemistry 36:10581-10594 (internal G.T).
Peyret et al. (1999) Biochemistry 38:3468-3477 (internal A.A/C.C/G.G/T.T).
Owczarzy et al. (2008) Biochemistry 47:5336-5353 (Mg2+ correction).
von Ahsen et al. (2001) Clin. Chem. 47:1956-1961 (Mg2+ -> Na+ equivalent).
"""
import os, re, json, subprocess, math

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.dirname(HERE)

# --------------------------------------------------------------------------
# SantaLucia (1998) unified NN, Table 1. Keys are 5'-XY-3' (top strand); the
# Watson-Crick complement 3'-X'Y'-5' is implied. 10 unique stacks. dH kcal/mol,
# dS cal/(mol.K).
# --------------------------------------------------------------------------
SANTALUCIA_1998 = {
    "AA": (-7.9, -22.2),   # = TT
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),   # = TG
    "GT": (-8.4, -22.4),   # = AC
    "CT": (-7.8, -21.0),   # = AG
    "GA": (-8.2, -22.2),   # = TC
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),   # = CC
}
# initiation contributions (per terminal base pair)
SL98_INIT_GC = (0.1, -2.8)
SL98_INIT_AT = (2.3, 4.1)

# --------------------------------------------------------------------------
# Owczarzy (2008) Mg2+ correction coefficients (their Eqs. for 1/Tm), pure
# divalent case. The general form:
#   1/Tm(Mg) = 1/Tm(1M Na) + a + b*ln[Mg] + fGC*(c + d*ln[Mg])
#              + (1/(2*(Nbp-1))) * (e + f*ln[Mg] + g*(ln[Mg])^2)
# with, in the presence of monovalent ions (0.22 <= R < 6, R=sqrt[Mg]/[Mon]),
# a, d, g rescaled per Owczarzy 2008 Eqs. 18-20.
# --------------------------------------------------------------------------
OWCZARZY_2008 = dict(a=3.92e-5, b=-9.11e-6, c=6.26e-5, d=1.42e-5,
                     e=-4.82e-4, f=5.25e-4, g=8.31e-5)


def owczarzy_ref(tmK_1M, freeMg, monM, fGC, Nbp):
    """Reference Owczarzy-2008 correction implemented from the PAPER (not from
    our JS). Returns Tm in Celsius. Mirrors the regimes in Owczarzy 2008."""
    a, b, c, d, e, f, g = (OWCZARZY_2008[k] for k in "abcdefg")
    R = math.sqrt(max(freeMg, 1e-9)) / max(monM, 1e-9)
    if R < 0.22:
        return tmK_1M - 273.15                      # monovalent-dominated: no divalent term
    if R < 6.0:                                     # competition regime: rescale a, d, g
        ln = math.log
        a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(monM) * ln(monM))
        d = 1.42e-5 * (1.279 - 4.03e-3 * ln(monM) - 8.03e-3 * ln(monM) ** 2)
        g = 8.31e-5 * (0.486 - 0.258 * ln(monM) + 5.25e-3 * ln(monM) ** 3)
    lm = math.log(freeMg)
    corr = a + b * lm + fGC * (c + d * lm) + (1.0 / (2 * (Nbp - 1))) * (e + f * lm + g * lm * lm)
    tmKc = 1.0 / (1.0 / tmK_1M + corr)
    return tmKc - 273.15


def vonahsen_naeq(na_mM, mg_mM, dntp_mM):
    """von Ahsen (2001): monovalent-equivalent [Na+] for a given [Mg2+],[dNTP]
    under a 1:1 Mg-dNTP chelation. Na_eq = Na + 120*sqrt(max(0, Mg - dNTP))."""
    mg_free = max(0.0, mg_mM - dntp_mM)
    return na_mM + 120.0 * math.sqrt(mg_free)


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def comp2(dn):
    """Base-complement of a dinucleotide (NOT reversed): comp('AT') -> 'TA'."""
    return _COMP[dn[0]] + _COMP[dn[1]]


def revcomp2(dn):
    """Reverse-complement of a dinucleotide: revcomp('CA') -> 'TG'."""
    return _COMP[dn[1]] + _COMP[dn[0]]


def biopython_nn_lookup(table, dn):
    """Look up dinucleotide `dn` (5'->3') in a BioPython DNA_NN* table, honouring
    NN symmetry. BioPython key is 'XY/comp(X)comp(Y)'. If the direct key is
    absent, the symmetric partner (reverse complement) carries the value."""
    k1 = dn + "/" + comp2(dn)
    if k1 in table:
        return table[k1]
    dn2 = revcomp2(dn)
    return table.get(dn2 + "/" + comp2(dn2))


def run_node_dump(script_name, stdin_text=None):
    """Run a node dumper in tools/tests/ and parse its JSON stdout."""
    p = subprocess.run(["node", os.path.join(HERE, script_name)],
                       input=stdin_text, text=True, capture_output=True, cwd=HERE)
    if p.returncode != 0:
        raise RuntimeError(f"{script_name} failed:\n{p.stderr}")
    return json.loads(p.stdout)


def parse_num_map(path, varname):
    """Parse a JS/Python object literal `varname = { KEY: number, ... }` into a
    dict, tolerant of quotes and whitespace. Used to compare the Python
    (concordance.py) and JS (mgspec.mjs) copies of the NN tables."""
    text = open(path, encoding="utf-8").read()
    i = text.find(varname)
    if i < 0:
        raise KeyError(f"{varname} not found in {path}")
    lo = text.find("{", i)
    hi = text.find("}", lo)
    body = text[lo + 1:hi]
    out = {}
    for k, v in re.findall(r"['\"]?([ACGT]{2})['\"]?\s*:\s*(-?\d+\.?\d*)", body):
        out[k] = float(v)
    return out
