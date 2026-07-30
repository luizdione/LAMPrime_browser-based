#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mg_offset.py - Re-score the published, experimentally validated SARS-CoV-2 spike
LAMP set (Prakash et al., 2023) under LAMPrime's own thermodynamics at two
magnesium conventions: NEB 8 mM (LAMPrime default) and PrimerExplorer 4 mM.

Purpose: show quantitatively that the ~4-6 C inner-primer offset is NOT a
magnesium effect. Within
LAMPrime's SantaLucia (1998) + von Ahsen (2001) model, LOWERING Mg to
PrimerExplorer's 4 mM lowers every Tm by ~2 C, i.e. WIDENS the gap to the design
windows. The offset therefore reflects the difference between PrimerExplorer's
Tm algorithm/target windows and LAMPrime's model (plus primer length), not the
reaction Mg2+ concentration.

Engine identical to tools/concordance.py (constants validated in tools/tests/).
Deterministic, standard library only.  Run:  python tools/mg_offset.py
"""
import math, os, sys
sys.stdout.reconfigure(encoding="utf-8")

NN_DH = {'AA':-7.9,'AT':-7.2,'TA':-7.2,'CA':-8.5,'GT':-8.4,'CT':-7.8,'GA':-8.2,'CG':-10.6,
         'GC':-9.8,'GG':-8.0,'AC':-8.4,'AG':-7.8,'TC':-8.2,'TG':-8.5,'TT':-7.9,'CC':-8.0}
NN_DS = {'AA':-22.2,'AT':-20.4,'TA':-21.3,'CA':-22.7,'GT':-22.4,'CT':-21.0,'GA':-22.2,'CG':-27.2,
         'GC':-24.4,'GG':-19.9,'AC':-22.4,'AG':-21.0,'TC':-22.2,'TG':-22.7,'TT':-22.2,'CC':-19.9}
R, NA, DNTP, OLIGO = 1.987, 50.0, 1.4, 50e-9


def tm(seq, mg):
    seq = ''.join(c for c in seq.upper() if c in 'ATGC'); N = len(seq)
    dH = dS = 0.0
    for i in range(N - 1):
        d = seq[i:i + 2]; dH += NN_DH[d]; dS += NN_DS[d]
    for b in (seq[0], seq[-1]):
        h, s = (0.1, -2.8) if b in 'GC' else (2.3, 4.1); dH += h; dS += s
    mgf = max(0.0, mg - DNTP)
    naeq = max(1e-3, (NA + 120 * math.sqrt(mgf))) / 1000.0
    dS_salt = dS + 0.368 * (N - 1) * math.log(naeq)
    return (dH * 1000) / (dS_salt + R * math.log(OLIGO / 4)) - 273.15


# Prakash et al. (2023) SARS-CoV-2 spike set; standard LAMP primer roles.
SARS = [('F3','TGGTGATATTGCTGCTAGA'), ('F2','CACCTTTGCTCACAGATG'),
        ('F1c','AGGTCCAACCAGAAGTGATT'), ('B1c','GCAGGTGCTGCATTACAA'),
        ('B2','TCTGTGTAACTCCAATACCA'), ('B3','GCACTATTAAATTGGTGGGC'),
        ('LF','GCTAACAGTGCAGAAGTGTATT'), ('LB','GCTATGCAAATGGCTTATAGGT')]
WIN = {'F3':(59,61),'B3':(59,61),'F2':(59,61),'B2':(59,61),
       'F1c':(64,66),'B1c':(64,66),'LF':(64,66),'LB':(64,66)}


def offset(t, lo, hi):
    return (lo - t) if t < lo else (0.0 if t <= hi else t - hi)


def main():
    lines = []
    def out(s=""):
        print(s); lines.append(s)
    out("SARS-CoV-2 spike set (Prakash 2023) re-scored under LAMPrime's model")
    out("at NEB 8 mM (default) vs PrimerExplorer 4 mM Mg2+.  dNTP 1.4 mM, Na 50 mM, oligo 50 nM.")
    out("")
    out(f"{'primer':5} {'window':8} {'Tm@8mM':>7} {'Tm@4mM':>7} {'d(4-8)':>7} {'off@8':>6} {'off@4':>6}")
    d_inner8 = d_inner4 = 0.0; n_inner = 0
    for name, seq in SARS:
        t8, t4 = tm(seq, 8.0), tm(seq, 4.0)
        lo, hi = WIN[name]
        out(f"{name:5} {str((lo,hi)):8} {t8:7.1f} {t4:7.1f} {t4-t8:7.2f} {offset(t8,lo,hi):6.1f} {offset(t4,lo,hi):6.1f}")
        if name in ('F1c', 'B1c'):
            d_inner8 += offset(t8, lo, hi); d_inner4 += offset(t4, lo, hi); n_inner += 1
    out("")
    out(f"Inner-primer (F1c,B1c) mean offset below window:  8 mM = {d_inner8/n_inner:.1f} C   4 mM = {d_inner4/n_inner:.1f} C")
    out("Conclusion: 4 mM Mg2+ LOWERS Tm (~2 C) and WIDENS the offset; the offset is a")
    out("Tm-model/convention difference (PrimerExplorer vs SantaLucia/von Ahsen), not a Mg2+ effect.")
    rep = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mg_offset_report.txt")
    open(rep, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    out("\nWrote " + rep)


if __name__ == "__main__":
    main()
