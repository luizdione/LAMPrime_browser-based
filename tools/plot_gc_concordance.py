#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""plot_gc_concordance.py — figura complementar do item #11.

Le tools/data/gc_concordance.csv (saida de gc_concordance.py; numeros REAIS) e
desenha 2 paineis mostrando que a Tm do LAMPrime concorda com uma implementacao
INDEPENDENTE (BioPython Tm_NN, NN4 = SantaLucia & Hicks 2004, sal Owczarzy 2004)
sobre um painel de conjuntos LAMP PUBLICADOS e desenhados por OUTRAS ferramentas
(PrimerExplorer V3/V4/V5, NEB), cobrindo ~30-67% de GC e 12 organismos:

  (A) Tm do LAMPrime vs Tm do BioPython, um ponto por oligo (linha de identidade);
  (B) ΔTm = LAMPrime - BioPython vs %GC do amplicon — a concordancia e plana e
      apertada (|ΔTm| <= ~0,15 C) em toda a faixa de GC.

Reproduzir:  python tools/gc_concordance.py && python tools/plot_gc_concordance.py
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "data", "gc_concordance.csv")
OUTDIR = os.path.join(HERE, "figures")
OUT = os.path.join(OUTDIR, "figS_tm_concordance_gc.png")


def load():
    rows = []
    with open(CSV, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if not r["tm_bio5"] or not r["amplicon_gc"]:
                continue
            rows.append(dict(assay=r["assay"], key=r["key"], gc=float(r["amplicon_gc"]),
                             lam=float(r["tm_lamprime"]), bio=float(r["tm_bio5"]),
                             d=float(r["d"]), oligo=r["oligo"]))
    return rows


def main():
    rows = load()
    if not rows:
        raise SystemExit("CSV vazio ou sem coluna bio5 — rode 'python tools/gc_concordance.py' com BioPython instalado.")
    os.makedirs(OUTDIR, exist_ok=True)

    gc = np.array([r["gc"] for r in rows])
    lam = np.array([r["lam"] for r in rows])
    bio = np.array([r["bio"] for r in rows])
    d = np.array([r["d"] for r in rows])
    nass = len({r["key"] for r in rows})
    r_pear = np.corrcoef(bio, lam)[0, 1]
    maxabs = np.max(np.abs(d))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.6, 5.0), gridspec_kw={"width_ratios": [1.0, 1.15]})
    cmap = plt.cm.viridis
    vmin, vmax = 28, 68

    # --- Panel A: identity ---
    lo, hi = min(bio.min(), lam.min()) - 2, max(bio.max(), lam.max()) + 2
    axA.plot([lo, hi], [lo, hi], "-", color="0.55", lw=1.0, zorder=1, label="identity")
    axA.scatter(bio, lam, c=gc, cmap=cmap, vmin=vmin, vmax=vmax, s=34,
                edgecolor="white", linewidth=0.4, zorder=3)
    axA.set_xlim(lo, hi); axA.set_ylim(lo, hi); axA.set_aspect("equal")
    axA.set_xlabel("BioPython $T_m$ (°C)   [NN4, Owczarzy 2004]", fontsize=9)
    axA.set_ylabel("LAMPrime $T_m$ (°C)", fontsize=9)
    axA.set_title("A  Per-oligo $T_m$ agreement", fontsize=10, loc="left")
    axA.text(0.04, 0.96, "n = %d oligos, %d assays\nPearson r = %.5f\nmax |Δ$T_m$| = %.2f °C"
             % (len(rows), nass, r_pear, maxabs),
             transform=axA.transAxes, va="top", ha="left", fontsize=8.2,
             bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.7", alpha=0.9))
    axA.tick_params(labelsize=8)

    # --- Panel B: delta vs GC ---
    axB.axhspan(-0.1, 0.1, color="0.90", zorder=0)
    axB.axhline(0, color="0.55", lw=1.0, zorder=1)
    sc = axB.scatter(gc, d, c=gc, cmap=cmap, vmin=vmin, vmax=vmax, s=34,
                     edgecolor="white", linewidth=0.4, zorder=3)
    # per-assay mean marker
    keys = {}
    for r in rows:
        keys.setdefault(r["key"], []).append(r)
    for k, rs in keys.items():
        gk = np.mean([x["gc"] for x in rs]); dk = np.mean([x["d"] for x in rs])
        axB.plot(gk, dk, marker="D", ms=7, mfc="none", mec="0.25", mew=1.1, zorder=4)
    axB.set_xlabel("amplicon GC content (%)", fontsize=9)
    axB.set_ylabel("Δ$T_m$ = LAMPrime − BioPython (°C)", fontsize=9)
    axB.set_title("B  Agreement is flat across GC (~30–67%)", fontsize=10, loc="left")
    ymax = max(0.2, np.max(np.abs(d)) * 1.25)
    axB.set_ylim(-ymax, ymax)
    axB.text(0.5, 0.965, "grey band = ±0.1 °C  ·  ◇ = per-assay mean",
             transform=axB.transAxes, va="top", ha="center", fontsize=8.0, color="0.35")
    axB.tick_params(labelsize=8)

    cbar = fig.colorbar(sc, ax=axB, pad=0.015, fraction=0.05)
    cbar.set_label("amplicon GC (%)", fontsize=8.5); cbar.ax.tick_params(labelsize=7.5)

    fig.suptitle("LAMPrime $T_m$ engine vs independent implementation across GC content — "
                 "13 published LAMP sets (PrimerExplorer / NEB), 12 organisms",
                 fontsize=10.5, y=1.005)
    plt.tight_layout()
    plt.savefig(OUT, dpi=300, bbox_inches="tight")
    print("saved", OUT)
    print("  n=%d oligos, %d assays, GC %.0f–%.0f%%, Pearson r=%.5f, max|Δ|=%.3f °C"
          % (len(rows), nass, gc.min(), gc.max(), r_pear, maxabs))


if __name__ == "__main__":
    main()
