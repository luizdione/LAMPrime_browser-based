#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_owczarzy_freemg.py - the Mg2+ correction (owczarzy2008) and the 1:1
Mg-dNTP chelation (freeMgM) in the engine reproduce Owczarzy et al. (2008):
  1. engine Tm(Mg) == a from-the-paper re-implementation over a grid;
  2. the low-Mg regime (R < 0.22) applies no divalent term;
  3. free Mg2+ obeys the 1:1 binding equilibrium (Ka = 3e4) and the limits.
"""
import json, math
import golden_literature as G


def _grid():
    g = []
    for tmK1M in (298.0, 315.0, 333.0, 350.0):
        for freeMg in (0.0008, 0.002, 0.004, 0.006, 0.008):
            for fGC in (0.30, 0.45, 0.60, 0.75):
                for N in (18, 20, 22, 25):
                    g.append(dict(tmK1M=tmK1M, freeMg=freeMg, monM=0.05, fGC=fGC,
                                  Nbp=N, totMg=freeMg + 0.0014, dntp=0.0014))
    return g


_RES = G.run_node_dump("_dump_owcz.mjs", json.dumps(_grid()))


def test_owczarzy_matches_literature_reimpl():
    maxd = max(abs(G.owczarzy_ref(r["tmK1M"], r["freeMg"], r["monM"], r["fGC"], r["Nbp"]) - r["tmC"])
               for r in _RES)
    assert maxd < 1e-9, f"engine Owczarzy vs literature reimpl max|d| = {maxd:.2e} C"


def test_low_mg_regime_has_no_divalent_term():
    # R = sqrt(freeMg)/monM < 0.22  ->  correction is zero (Tm == 1 M-Na Tm)
    grid = [dict(tmK1M=330.0, freeMg=1e-5, monM=0.30, fGC=0.5, Nbp=20,
                 totMg=1e-5 + 0.0014, dntp=0.0014)]
    r = G.run_node_dump("_dump_owcz.mjs", json.dumps(grid))[0]
    assert math.sqrt(r["freeMg"]) / r["monM"] < 0.22
    assert abs(r["tmC"] - (r["tmK1M"] - 273.15)) < 1e-9


def test_freemg_equilibrium_and_limits():
    Ka = 3e4
    for tot, dntp in [(0.008, 0.0014), (0.002, 0.0016), (0.010, 0.0), (0.004, 0.004)]:
        r = G.run_node_dump("_dump_owcz.mjs",
                            json.dumps([dict(tmK1M=330.0, freeMg=0.005, monM=0.05,
                                             fGC=0.5, Nbp=20, totMg=tot, dntp=dntp)]))[0]
        free = r["freeMgOut"]
        assert 0 <= free <= tot + 1e-12
        if dntp == 0:
            assert abs(free - tot) < 1e-12, "no dNTP -> free == total"
        else:
            assert free < tot, "dNTP present -> free < total"
            bound = tot - free                       # [MgdNTP]
            free_dntp = dntp - bound
            if free_dntp > 1e-9:                      # check 1:1 binding: Ka = bound/(free*free_dntp)
                assert abs(Ka - bound / (free * free_dntp)) / Ka < 1e-6
