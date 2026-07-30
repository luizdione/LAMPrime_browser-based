#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_fig3.py - Figure 3 do LAMPrime: especificidade x Mg2+ (exemplo Dengue).
Le tools/data/dengue_mg_Sij.csv (saida de dengue_mgspec.mjs; numeros REAIS) e
desenha 2 paineis:
  (A) heatmap da margem de seguranca do off-target = T_reacao - Tm_off (Mg-dependente);
  (B) Tm_off de cada primer vs MgSO4 total, mostrando que ficam abaixo da T de reacao.
Mensagem honesta: o conjunto DENV-1 (Hu 2015) e robustamente especifico contra DENV-2/3/4
em toda a faixa 2-10 mM; o Mg desloca o Tm absoluto, nao a discriminacao.
Reproduzir: python tools/plot_fig3.py
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, 'data', 'dengue_mg_Sij.csv')
OUT = os.path.normpath(os.path.join(HERE, '..', 'PAPER', 'figuras', 'fig3_mg_specificity.png'))
TRXN = 63.0  # temperatura de reacao LAMP (C)

rows = list(csv.DictReader(open(CSV)))
order, seen = [], set()
for r in rows:
    p = r['primer'].strip()
    if p not in seen:
        seen.add(p); order.append(p)
mgs = sorted({float(r['mg_total']) for r in rows})
freemap = {float(r['mg_total']): float(r['mg_free']) for r in rows}
mm = {}
for r in rows:
    mm.setdefault(r['primer'].strip(), r.get('n_mismatch', '').strip())

def grid(col):
    G = np.full((len(order), len(mgs)), np.nan)
    for r in rows:
        i = order.index(r['primer'].strip()); j = mgs.index(float(r['mg_total']))
        G[i, j] = float(r[col])
    return G

TmOff = grid('tm_off')
TmOn  = grid('tm_on')
Safety = TRXN - TmOff

fig, (axA, axB) = plt.subplots(1, 2, figsize=(15, 4.8), gridspec_kw={'width_ratios': [1.45, 1.0]})

# ---- Painel A: heatmap da margem de seguranca ----
norm = TwoSlopeNorm(vmin=-5, vcenter=0, vmax=30)
im = axA.imshow(Safety, aspect='auto', cmap='RdYlGn', norm=norm)
axA.set_xticks(range(len(mgs)))
axA.set_xticklabels(['%.1f\n(%.1f)' % (m, freemap[m]) for m in mgs], fontsize=6.5)
axA.set_yticks(range(len(order)))
axA.set_yticklabels(['%s (%s mm)' % (p, mm.get(p, '?')) for p in order], fontsize=9)
axA.set_xlabel('total MgSO$_4$ (mM)  [free Mg$^{2+}$, mM]', fontsize=9)
axA.set_ylabel('DENV-1 primer (mismatches to nearest off-target)', fontsize=9)
axA.set_title('A  off-target safety margin  $T_{rxn}-T_m^{off}$ ($^{\\circ}$C)', fontsize=10, loc='left')
for i in range(len(order)):
    for j in range(len(mgs)):
        axA.text(j, i, '%.0f' % TmOff[i, j], ha='center', va='center', fontsize=6, color='black')
cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.02)
cb.set_label('$T_{rxn}-T_m^{off}$ ($^{\\circ}$C)  ( >0 safe )', fontsize=8)
axA.text(0.0, -0.32, 'cell = off-target $T_m$ ($^{\\circ}$C); all far below the %g $^{\\circ}$C reaction temp $\\Rightarrow$ robustly specific'
         % TRXN, transform=axA.transAxes, fontsize=7.5, color='#444')

# ---- Painel B: Tm_off vs Mg (todos os primers) + linha da T de reacao ----
colors = plt.cm.viridis(np.linspace(0, 0.9, len(order)))
for i, p in enumerate(order):
    axB.plot(mgs, TmOff[i, :], '-o', ms=3, lw=1.4, color=colors[i], label='%s' % p)
top = TRXN + 8
axB.axhline(TRXN, color='crimson', ls='--', lw=1.4)
axB.axhspan(TRXN, top, color='crimson', alpha=0.06)
axB.text(mgs[-1], TRXN + 0.5, 'reaction temp %g $^{\\circ}$C  (off-target priming risk above)' % TRXN,
         ha='right', va='bottom', fontsize=7.5, color='crimson')
axB.set_xlabel('total MgSO$_4$ (mM)', fontsize=9)
axB.set_ylabel('off-target $T_m$ ($^{\\circ}$C)', fontsize=9)
axB.set_title('B  off-target $T_m$ rises with Mg$^{2+}$ but stays below $T_{rxn}$', fontsize=10, loc='left')
axB.set_ylim(min(np.nanmin(TmOff) - 4, 10), top)
axB.legend(fontsize=7, ncol=2, loc='lower right', framealpha=0.9)
axB.grid(alpha=0.25)

plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches='tight')
print('saved', OUT)
print('safety margin min/max: %.1f / %.1f C  (>0 = safe)' % (np.nanmin(Safety), np.nanmax(Safety)))
print('off-target Tm min/max: %.1f / %.1f C  (reaction temp %g)' % (np.nanmin(TmOff), np.nanmax(TmOff), TRXN))
