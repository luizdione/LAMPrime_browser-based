#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
benchmark_independent.py - Verificacao INDEPENDENTE do motor de Tm do LAMPrime.

Cruza o Tm do LAMPrime (modelo NN de SantaLucia 1998 + Mg2+ livre por quelacao,
via concordance.tm()) contra uma implementacao de terceiros — Biopython
(Bio.SeqUtils.MeltingTemp.Tm_NN) — sobre:
  - os conjuntos de primers LAMP PUBLICADOS (A. marginale msp1b, SARS-CoV-2 S,
    M. tuberculosis IS6110), que sao a base da concordancia do LAMPrime;
  - janelas amostradas do alvo AT-rico *P. falciparum* 18S (M19172.1) — item #9 —
    para checar o acordo tambem no regime de baixo GC.

Objetivo: um terceiro independente confirmar que o nucleo NN do LAMPrime esta
correto; diferencas residuais localizam-se no modelo de sal/Mg (esperado e
documentado), nao no NN.

Convencao de concentracao: Biopython com dnac1=dnac2=25 nM dá conc. efetiva
25 - 25/2 = 12,5 nM para oligo nao autocomplementar, IDENTICA a convencao do
LAMPrime (oligo_total 50 nM, termo [C]/4 = 12,5 nM). Assim a comparacao isola o
modelo de sal, nao a definicao de concentracao.

Biopython saltcorr: 5 (Owczarzy 2004, monovalente + Mg) é o mais proximo do
caminho de sal do LAMPrime (von Ahsen Mg->Na); 7 (Owczarzy 2008, correcao
divalente completa com dNTP) é um segundo ponto independente que mostra a
dispersao entre modelos de sal.

Conjuntos PrimerDigital: quando o autor fornecer os arquivos, adicione-os em
EXTRA_SETS (ver abaixo) — o benchmark passa a inclui-los automaticamente.

Refs: SantaLucia 1998; SantaLucia & Hicks 2004; Owczarzy 2004/2008; von Ahsen 2001.
Uso:  python tools/benchmark_independent.py [--windows 30] [--wlen 20]
"""
import os, sys, math, argparse, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import concordance as C  # motor + primers publicados (fonte unica)

sys.stdout.reconfigure(encoding='utf-8')

try:
    from Bio.SeqUtils import MeltingTemp as mt
    HAVE_BIO = True
except Exception:
    HAVE_BIO = False


def bio_tm(seq, saltcorr):
    """Tm independente (Biopython Tm_NN, tabela NN4 = SantaLucia & Hicks 2004),
    nas condicoes do LAMPrime (Na 50, Mg 8, dNTP 1,4; conc. efetiva 12,5 nM)."""
    return mt.Tm_NN(seq, nn_table=mt.DNA_NN4, Na=C.NA, Mg=C.MG, dNTPs=C.DNTP,
                    dnac1=25, dnac2=25, saltcorr=saltcorr)


# conjuntos publicados (fonte: concordance.py). Adicione conjuntos PrimerDigital aqui.
PUBLISHED = [
    ('A. marginale msp1b (Giglioti 2018)', C.amarginale_msp1b),
    ('SARS-CoV-2 gene S (Prakash 2023)', C.sarscov2),
    ('M. tuberculosis IS6110 (Bentaleb 2016)', C.mtb_is6110),
]
EXTRA_SETS = []  # ex.: ('PrimerDigital <alvo>', {'F3':..., 'B3':..., 'FIP':..., 'BIP':..., 'LF':..., 'LB':...})


def oligos_of(P):
    order = ['F3', 'F2', 'F1c', 'B1c', 'B2', 'B3', 'FIP', 'BIP', 'LF', 'LB']
    out = []
    for k in order:
        if k in ('FIP', 'BIP') and P.get('F2') and P.get('B2'):
            continue
        if P.get(k):
            out.append((k, P[k]))
    return out


def sample_windows(seq, wlen, n):
    """Amostra ~n janelas de comprimento wlen ao longo de seq (passo uniforme),
    so ATGC. Regime AT-rico do P. falciparum 18S (#9)."""
    seq = ''.join(c for c in seq.upper() if c in 'ATGC')
    if len(seq) < wlen:
        return []
    positions = list(range(0, len(seq) - wlen + 1))
    step = max(1, len(positions) // max(1, n))
    return [seq[p:p + wlen] for p in positions[::step][:n]]


def summarize(deltas, label):
    if not deltas:
        print('  %s: (sem dados)' % label); return
    ad = [abs(d) for d in deltas]
    print('  %s: n=%d  media d=%+.2f  |d|: media %.2f, mediana %.2f, max %.2f C'
          % (label, len(deltas), statistics.mean(deltas), statistics.mean(ad),
             statistics.median(ad), max(ad)))


def main():
    ap = argparse.ArgumentParser(description='Benchmark independente de Tm (LAMPrime vs Biopython).')
    ap.add_argument('--windows', type=int, default=30, help='n de janelas amostradas do alvo AT-rico (#9)')
    ap.add_argument('--wlen', type=int, default=20, help='comprimento das janelas amostradas')
    args = ap.parse_args()

    print('benchmark_independent.py — LAMPrime (concordance.tm) vs Biopython Tm_NN (NN4=SantaLucia&Hicks 2004)')
    if not HAVE_BIO:
        print('\n  Biopython AUSENTE. Instale com:  python -m pip install biopython')
        print('  (o benchmark precisa de Bio.SeqUtils.MeltingTemp). Abortando sem erro.')
        return
    print('  cond.: Na %.0f mM, Mg %.0f mM, dNTP %.1f mM, conc.efetiva 12,5 nM (dnac1=dnac2=25).' % (C.NA, C.MG, C.DNTP))
    print('  colunas: LAMP=concordance.tm ; bio5=saltcorr Owczarzy 2004 ; bio7=saltcorr Owczarzy 2008')

    d5_all, d7_all = [], []
    for name, P in (PUBLISHED + EXTRA_SETS):
        print('\n== %s ==' % name)
        print('  %-5s %3s %7s %8s %8s %8s %8s' % ('oligo', 'len', 'LAMP', 'bio5', 'bio7', 'd(LAMP-bio5)', 'd(LAMP-bio7)'))
        for nm, seq in oligos_of(P):
            lamp = C.tm(seq); b5 = bio_tm(seq, 5); b7 = bio_tm(seq, 7)
            d5, d7 = lamp - b5, lamp - b7
            d5_all.append(d5); d7_all.append(d7)
            print('  %-5s %3d %7.1f %8.1f %8.1f %+8.2f %+8.2f' % (nm, len(seq), lamp, b5, b7, d5, d7))

    # #9: regime AT-rico (P. falciparum 18S)
    pf = C.load_fasta('pfalciparum_18s_M19172.fasta')
    print('\n== janelas AT-ricas — P. falciparum 18S M19172.1 (#9) ==')
    if not pf:
        print('  FASTA ausente (tools/data/pfalciparum_18s_M19172.fasta).')
    else:
        wins = sample_windows(pf, args.wlen, args.windows)
        gcs = [100 * sum(c in 'GC' for c in w) / len(w) for w in wins]
        dp5, dp7 = [], []
        for w in wins:
            lamp = C.tm(w)
            if math.isnan(lamp):
                continue
            b5 = bio_tm(w, 5); b7 = bio_tm(w, 7)
            dp5.append(lamp - b5); dp7.append(lamp - b7)
        print('  %d janelas de %d nt | GC medio %.1f%% (faixa %.0f–%.0f%%) — regime AT-rico'
              % (len(wins), args.wlen, statistics.mean(gcs), min(gcs), max(gcs)))
        summarize(dp5, 'LAMP-bio5 (AT-rico)')
        summarize(dp7, 'LAMP-bio7 (AT-rico)')

    print('\n== resumo global (conjuntos publicados) ==')
    summarize(d5_all, 'LAMP-bio5 (Owczarzy 2004)')
    summarize(d7_all, 'LAMP-bio7 (Owczarzy 2008)')
    print('\n  interpretacao: bio5 (modelo de sal proximo ao von Ahsen do LAMPrime) deve concordar')
    print('  em ~0,1–1 C; bio7 (correcao divalente completa) diverge mais por escolha de modelo de sal,')
    print('  NAO pelo nucleo NN. Acordo tambem no regime AT-rico valida o motor fora da faixa usual de GC.')


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        try:
            sys.stdout.close()
        except Exception:
            pass
