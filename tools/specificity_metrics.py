#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
specificity_metrics.py - Metricas quantitativas da triagem de especificidade (#13).

Transforma a triagem de reatividade cruzada do LAMPrime (flag / nao-flag) em
metricas: sensibilidade, especificidade, matriz de confusao e curva ROC (com AUC),
sobre uma verdade-terreno baseada em SEQUENCIA:
  - POSITIVO: primer contra o SEU PROPRIO alvo (deve parear/amplificar);
  - NEGATIVO: primer contra sequencias de OUTROS organismos (nao deve parear).

Reutiliza o motor e os conjuntos de primers publicados de concordance.py (fonte
unica) e reproduz a triagem do app: semente 3' exata de k nt + extensao gap-free,
melhor janela por numero de mismatches, nas duas fitas (mirror de tools/crossreact.js
e app.js screenPrimerOnBg). Offline/deterministico (stdlib).

LIMITE (ROC PARCIAL): a verdade-terreno aqui e por IDENTIDADE DE SEQUENCIA, um proxy
do pareamento — NAO os `.eds` do autor (saida PrimerExplorer = referencia experimental).
A ROC completa exige os `.eds`; ao fornece-los, ver ROC_REFERENCE abaixo.

Refs: Notomi 2000; especificidade LAMPrime (README); crossreact.js.
Uso:  python tools/specificity_metrics.py [--seed 13] [--maxmm 2]
"""
import os, sys, math, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import concordance as C

sys.stdout.reconfigure(encoding='utf-8')

# --- painel: conjuntos de primers (com seu alvo) + fundos de outros organismos ---
SETS = {
    'A. marginale msp1b': dict(own='amarginale_msp1b_synthetic.fasta', primers=C.amarginale_msp1b, exclude=set()),
    'SARS-CoV-2 S':       dict(own='sarscov2_spike_NC045512.2_21563-25384.fasta', primers=C.sarscov2, exclude=set()),
    # M. bovis (complexo MTB) carrega IS6110 -> deteccao ESPERADA, excluida dos negativos.
    'M. tuberculosis IS6110': dict(own='mtb_is6110_X17348.fasta', primers=C.mtb_is6110,
                                   exclude={'mbovis_genome_NC002945.fasta'}),
}
BACKGROUNDS = [
    'acentrale_genome_CP001759.fasta', 'aphago_genome_CM177739.fasta', 'aovis_genome_CP015994.fasta',
    'mbovis_genome_NC002945.fasta', 'bbovis_18s_L19077.fasta', 'bbigemina_18s_KP710228.fasta',
    'ratg13_MN996532.fasta', 'sarscov_NC004718.fasta', 'pfalciparum_18s_M19172.fasta',
    'zika_NC035889.fasta', 'denv2_AF038403.fasta',
]

# preenchido pelo autor com os .eds de referencia (PrimerExplorer) para a ROC completa.
ROC_REFERENCE = None  # ex.: {'A. marginale msp1b': {'F3': True/False, ...}, ...}


def components(P, target):
    """8 primers de ligacao (F3,F2,F1c,B1c,B2,B3,LF,LB), como em concordance.report:
    usa componentes explicitos; senao divide FIP/BIP contra o alvo proprio."""
    f1c, f2 = (P['F1c'], P['F2']) if P.get('F1c') and P.get('F2') else C.split_fip(target, P['FIP'])
    b1c, b2 = (P['B1c'], P['B2']) if P.get('B1c') and P.get('B2') else C.split_bip(target, P['BIP'])
    out = [('F3', P['F3']), ('F2', f2), ('F1c', f1c), ('B1c', b1c), ('B2', b2), ('B3', P['B3'])]
    if P.get('LF'): out.append(('LF', P['LF']))
    if P.get('LB'): out.append(('LB', P['LB']))
    return [(n, s) for n, s in out if s]


def count_mm(pat, bg, start):
    """mismatches de pat (comprimento total) alinhado gap-free em bg a partir de start; None se fora."""
    Lp = len(pat)
    if start < 0 or start + Lp > len(bg):
        return None
    mm = 0
    for i in range(Lp):
        if bg[start + i] != pat[i]:
            mm += 1
    return mm


def best_mm(primer, bg, k):
    """Menor n de mismatches do primer (comprimento total) ancorado por semente 3' de k nt
    (extensao gap-free), nas duas fitas. None se nao ha semente. Mirror de screenPrimerOnBg."""
    Lp = len(primer)
    if Lp < k or len(bg) < k:
        return None
    best = None
    # fita '+': semente = ultimos k do primer; janela recua (Lp-k)
    seed = primer[Lp - k:]
    pos = bg.find(seed)
    while pos != -1:
        mm = count_mm(primer, bg, pos - (Lp - k))
        if mm is not None and (best is None or mm < best):
            best = mm
        pos = bg.find(seed, pos + 1)
    # fita '-': semente = primeiros k do revcomp
    rcp = C.rc(primer); rseed = rcp[:k]
    pos = bg.find(rseed)
    while pos != -1:
        mm = count_mm(rcp, bg, pos)
        if mm is not None and (best is None or mm < best):
            best = mm
        pos = bg.find(rseed, pos + 1)
    return best


def load(fname):
    return C.load_fasta(fname)


def auc_trapezoid(points):
    """AUC por trapezio; points = lista de (fpr, tpr) ja ordenavel por fpr."""
    pts = sorted(set(points))
    if len(pts) < 2:
        return float('nan')
    a = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        a += (x1 - x0) * (y0 + y1) / 2.0
    return a


def main():
    ap = argparse.ArgumentParser(description='Metricas de especificidade (sens/spec/ROC) do screen LAMPrime.')
    ap.add_argument('--seed', type=int, default=13, help='semente 3\' do ponto de operacao (nt); padrao 13 (clamp 8..20)')
    ap.add_argument('--maxmm', type=int, default=2, help='mismatches maximos no ponto de operacao; padrao 2')
    ap.add_argument('--roc-seed', type=int, default=8, help='semente 3\' permissiva p/ a ROC (nt); padrao 8 (piso do app)')
    args = ap.parse_args()
    k = max(8, min(20, args.seed)); kr = max(8, min(20, args.roc_seed))

    print('specificity_metrics.py — sens/spec/ROC da triagem de especificidade (#13)')
    print('  verdade-terreno por IDENTIDADE DE SEQUENCIA (proxy). ROC completa exige os .eds do autor (ROC_REFERENCE).')
    print('  screen: semente 3\' + extensao gap-free, duas fitas (mirror de crossreact.js/app.js).')
    print('  ponto de operacao: semente k=%d, maxMM=%d | ROC: semente permissiva k=%d, varrendo maxMM.' % (k, args.maxmm, kr))

    seqs = {}
    for f in set(BACKGROUNDS) | {v['own'] for v in SETS.values()}:
        s = load(f)
        if s:
            seqs[f] = s
        else:
            print('  AVISO: FASTA ausente/vazio: %s' % f)

    # instancias: (setname, primer, seqname, label_pos, bm_op, bm_roc)
    instances = []
    variants = []  # own-target primers sem semente k (variantes primer-vs-referencia; idem concordance NAO ACHADO)
    print('\n== varredura por conjunto ==')
    for setname, cfg in SETS.items():
        own = cfg['own']
        if own not in seqs:
            print('  %s: alvo proprio ausente (%s) — pulado' % (setname, own)); continue
        prims = components(cfg['primers'], seqs[own])
        neg_files = [f for f in seqs if f != own and f not in cfg['exclude']]
        pos_flag = 0
        for nm, p in prims:
            bo = best_mm(p, seqs[own], k); br = best_mm(p, seqs[own], kr)
            instances.append((setname, nm, own, True, bo, br))
            if bo is not None and bo <= args.maxmm:
                pos_flag += 1
            else:
                variants.append((setname, nm))
        neg_hits = []
        for nm, p in prims:
            for f in neg_files:
                bo = best_mm(p, seqs[f], k); br = best_mm(p, seqs[f], kr)
                instances.append((setname, nm, f, False, bo, br))
                if bo is not None and bo <= args.maxmm:
                    neg_hits.append((nm, f, bo))
        print('  %-26s primers=%d | POS sinalizados %d/%d | fundos=%d | cross-hits=%d %s'
              % (setname, len(prims), pos_flag, len(prims), len(neg_files), len(neg_hits),
                 ('[%s]' % ', '.join('%s@%s(mm%d)' % (n, f.split('_')[0], m) for n, f, m in neg_hits[:6])) if neg_hits else ''))

    # ---- ponto de operacao (semente=k, maxmm) ----
    def flagged(bm):
        return bm is not None and bm <= args.maxmm
    TP = sum(1 for _, _, _, pos, bo, _ in instances if pos and flagged(bo))
    FN = sum(1 for _, _, _, pos, bo, _ in instances if pos and not flagged(bo))
    FP = sum(1 for _, _, _, pos, bo, _ in instances if not pos and flagged(bo))
    TN = sum(1 for _, _, _, pos, bo, _ in instances if not pos and not flagged(bo))
    P, Ncnt = TP + FN, TN + FP
    sens = TP / P if P else float('nan')
    spec = TN / Ncnt if Ncnt else float('nan')
    sens_screen = TP / (P - len(variants)) if (P - len(variants)) else float('nan')
    print('\n== ponto de operacao (semente=%d, maxMM=%d) ==' % (k, args.maxmm))
    print('  positivos (primer vs alvo proprio) = %d | negativos (vs outro organismo) = %d' % (P, Ncnt))
    print('  TP=%d FN=%d FP=%d TN=%d' % (TP, FN, FP, TN))
    print('  sensibilidade (a referencia) = %.3f | especificidade = %.3f' % (sens, spec))
    if variants:
        print('  * %d FN sao VARIANTES primer-vs-referencia (idem concordance NAO ACHADO): %s'
              % (len(variants), ', '.join('%s/%s' % (s.split()[0], n) for s, n in variants)))
        print('    sensibilidade do SCREEN (primers presentes na referencia) = %.3f' % sens_screen)

    # ---- ROC: semente permissiva kr, varre maxMM 0..N -> traca FPR de 0 a ~1 ----
    thresholds = list(range(-1, 16))
    roc = []
    for tau in thresholds:
        tp = sum(1 for _, _, _, pos, _, br in instances if pos and br is not None and br <= tau)
        fp = sum(1 for _, _, _, pos, _, br in instances if not pos and br is not None and br <= tau)
        roc.append((fp / Ncnt if Ncnt else 0.0, tp / P if P else 0.0, tau))
    pts = [(f, t) for f, t, _ in roc]
    pauc = auc_trapezoid(pts)
    fpr_max = max((f for f, _ in pts), default=0.0)
    pauc_norm = pauc / fpr_max if fpr_max > 0 else float('nan')  # TPR medio sobre a faixa de FPR acessivel
    print('\n== ROC (semente permissiva=%d, varrendo maxMM 0..%d) ==' % (kr, max(thresholds)))
    print('  %5s %7s %7s' % ('maxMM', 'FPR', 'TPR'))
    for fpr, tpr, tau in roc:
        if tau < 0:
            continue
        print('  %5d %7.3f %7.3f' % (tau, fpr, tpr))
    print('  pAUC (parcial) = %.4f sobre FPR em [0, %.3f]; TPR medio nessa faixa = %.3f' % (pauc, fpr_max, pauc_norm))
    print('  (a semente 3\' ancora o screen como FILTRO: FPR satura < 1 — nao ha AUC[0,1] natural.')
    print('   O resumo do desempenho e o PONTO DE OPERACAO acima; a ROC definitiva exige os .eds do autor.)')

    print('\n  interpretacao: no ponto de operacao do app (semente 13, maxMM 2) o screen opera no canto de')
    print('  FPR baixissimo (alta especificidade); a ROC (semente permissiva) mostra a curva completa.')
    print('  Cross-hits reais e esperados: primers do S de SARS-CoV-2 batem em RaTG13 (~96%% identico).')
    print('  M. bovis excluido dos negativos de IS6110 (deteccao esperada — complexo MTB).')
    print('  ROC definitiva: preencher ROC_REFERENCE com os .eds do autor e comparar.')


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        try:
            sys.stdout.close()
        except Exception:
            pass
