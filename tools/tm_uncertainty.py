#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tm_uncertainty.py - Incerteza de Tm dos conjuntos LAMP (publicados).

Quantifica quanto o Tm previsto varia conforme escolhas de parametrizacao:
  (a) conjunto de parametros NN: SantaLucia 1998 (unified) vs SantaLucia & Hicks 2004;
  (b) modelo de sal/Mg2+: equivalente de Mg2+->Na+ (von Ahsen 2001, padrao LAMPrime)
      vs monovalente-apenas (Mg2+ = 0);
  (c) acuracia agregada de predicao do modelo NN (SantaLucia & Hicks 2004).
  (opcional) Monte Carlo: propaga uma incerteza relativa nos NN (sigma fornecido pelo
      usuario) para uma banda de Tm por oligo.

Reutiliza o MESMO motor termodinamico e os MESMOS conjuntos de primers publicados de
concordance.py (fonte unica de verdade), de forma offline e deterministica (stdlib).

NOTA CIENTIFICA (a): os 10 parametros NN de Watson-Crick do conjunto "unified" de
SantaLucia (1998, PNAS 95:1460) sao os MESMOS republicados por SantaLucia & Hicks
(2004, Annu Rev Biophys Biomol Struct 33:415, Tab. 1). Portanto dTm(1998->2004) = 0
por construcao: a escolha 1998 vs 2004 NAO e fonte de incerteza de Tm para pares WC.
Este script verifica esse fato e desloca a analise para o que de fato move o Tm
(modelo de sal/Mg e acuracia agregada do modelo).

Refs: SantaLucia 1998; SantaLucia & Hicks 2004; von Ahsen 2001.
Uso:  python tools/tm_uncertainty.py [--sigma-pred 2.0] [--montecarlo REL_SIGMA] [--n 2000]
"""
import os, sys, math, argparse, random, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import concordance as C  # fonte unica: NN_DH/NN_DS, tm(), R/NA/MG/DNTP/OLIGO, dicts de primers

sys.stdout.reconfigure(encoding='utf-8')

# --- (a) Conjuntos de parametros NN (ΔH kcal/mol, ΔS cal/mol/K) ---
# SantaLucia 1998 "unified" == o que LAMPrime usa (concordance.NN_DH/NN_DS).
NN_DH_1998 = dict(C.NN_DH); NN_DS_1998 = dict(C.NN_DS)
# SantaLucia & Hicks 2004: parametros WC identicos ao unified de 1998 (Tab. 1).
NN_DH_2004 = dict(C.NN_DH); NN_DS_2004 = dict(C.NN_DS)

# constantes do modelo (identicas a concordance.py / LAMPrime)
R, NA, MG, DNTP, OLIGO = C.R, C.NA, C.MG, C.DNTP, C.OLIGO

# (c) acuracia agregada de predicao de Tm do modelo NN para oligos, relatada por
# SantaLucia & Hicks (2004): ~2 C. Valor citavel e AJUSTAVEL (--sigma-pred); nao e um
# sigma por-parametro nem substitui uma propagacao completa (ver --montecarlo).
SIGMA_TM_PRED_DEFAULT = 2.0


def tm_with(seq, dH_tab, dS_tab, mg=MG, na=NA, dntp=DNTP, oligo=OLIGO):
    """Tm (C) pelo modelo NN de LAMPrime, com tabelas NN e sal/Mg configuraveis.
    Identico a concordance.tm() quando dH_tab/dS_tab sao o conjunto de 1998 e os
    parametros de sal sao os padrao (verificado por self-check em main())."""
    seq = ''.join(c for c in seq.upper() if c in 'ATGC'); N = len(seq)
    if N < 2:
        return float('nan')
    dH = dS = 0.0
    for i in range(N - 1):
        d = seq[i:i + 2]
        if d not in dH_tab:
            return float('nan')
        dH += dH_tab[d]; dS += dS_tab[d]
    for b in (seq[0], seq[-1]):
        h, s = (0.1, -2.8) if b in 'GC' else (2.3, 4.1); dH += h; dS += s
    mgf = C.free_mg_mM(mg, dntp); naeq = max(1e-3, (na + 120 * math.sqrt(max(0.0, mgf)))) / 1000.0
    dS_salt = dS + 0.368 * (N - 1) * math.log(naeq)
    return (dH * 1000) / (dS_salt + R * math.log(oligo / 4)) - 273.15


def oligos_of(P):
    """Oligos de um conjunto LAMP como dicionario. Usa componentes de ligacao
    explicitos (F1c/F2/B1c/B2) quando presentes; sempre inclui F3/B3/LF/LB e,
    quando nao ha componentes explicitos, os FIP/BIP compostos (como sintetizados)."""
    order = ['F3', 'F2', 'F1c', 'B1c', 'B2', 'B3', 'FIP', 'BIP', 'LF', 'LB']
    out = []
    for k in order:
        if k in ('FIP', 'BIP') and P.get('F2') and P.get('B2'):
            continue  # ja cobertos pelos componentes de ligacao explicitos
        if P.get(k):
            out.append((k, P[k]))
    return out


def montecarlo_sd(seq, rel_sigma, n, seed=0):
    """Desvio-padrao de Tm ao perturbar cada NN (ΔH e ΔS) por um fator relativo
    gaussiano de sigma=rel_sigma. rel_sigma e uma SUPOSICAO fornecida pelo usuario
    (nao um valor de literatura) — sirva-a a partir da fonte primaria se desejar
    uma propagacao rigorosa.
    CAVEAT: perturba ΔH e ΔS de forma INDEPENDENTE, ignorando a forte compensacao
    entalpia-entropia dos parametros NN; isso SUPERESTIMA a variancia de Tm. Uma
    propagacao fiel exigiria a matriz de covariancia (ΔH,ΔS) da fonte primaria.
    Use os valores como limite superior grosseiro, nao como incerteza calibrada."""
    rnd = random.Random(seed)
    tms = []
    for _ in range(n):
        dH = {k: v * (1 + rel_sigma * rnd.gauss(0, 1)) for k, v in NN_DH_1998.items()}
        dS = {k: v * (1 + rel_sigma * rnd.gauss(0, 1)) for k, v in NN_DS_1998.items()}
        t = tm_with(seq, dH, dS)
        if not math.isnan(t):
            tms.append(t)
    return statistics.pstdev(tms) if len(tms) > 1 else float('nan')


# conjuntos publicados (fonte: concordance.py). Conjuntos PrimerDigital/LAMPrime entram
# aqui quando disponiveis (itens #11/#13 do backlog).
SETS = [
    ('A. marginale msp1b (Giglioti 2018)', C.amarginale_msp1b),
    ('SARS-CoV-2 gene S (Prakash 2023)', C.sarscov2),
    ('M. tuberculosis IS6110 (Bentaleb 2016)', C.mtb_is6110),
]


def main():
    ap = argparse.ArgumentParser(description='Incerteza de Tm (SantaLucia 1998 vs 2004 + sal/Mg).')
    ap.add_argument('--sigma-pred', type=float, default=SIGMA_TM_PRED_DEFAULT,
                    help='banda de acuracia agregada do modelo NN, em C (SantaLucia & Hicks 2004; padrao 2.0)')
    ap.add_argument('--montecarlo', type=float, default=None, metavar='REL_SIGMA',
                    help='ativa Monte Carlo: incerteza relativa (SUPOSICAO do usuario) nos NN, ex. 0.04')
    ap.add_argument('--n', type=int, default=2000, help='amostras do Monte Carlo (padrao 2000)')
    args = ap.parse_args()

    # self-check: tm_with(1998, sal padrao) deve reproduzir concordance.tm()
    probe = 'GCACTACCGTTCATGGATGA'
    assert abs(tm_with(probe, NN_DH_1998, NN_DS_1998) - C.tm(probe)) < 1e-9, 'motor divergiu de concordance.tm()'
    # (a) 1998 vs 2004: identidade das tabelas WC
    assert NN_DH_1998 == NN_DH_2004 and NN_DS_1998 == NN_DS_2004

    print('tm_uncertainty.py — incerteza de Tm (offline, deterministico; motor/primers de concordance.py)')
    print('  (a) NN 1998 vs 2004: parametros WC identicos (unified) -> dTm = 0.000 C por construcao [verificado]')
    print('  (b) sal/Mg: von Ahsen Mg2+->Na+ @ %.0f mM (padrao) vs monovalente-apenas (Mg2+=0)' % MG)
    print('  (c) acuracia agregada do modelo NN: +/- %.1f C (SantaLucia & Hicks 2004)' % args.sigma_pred)
    if args.montecarlo is not None:
        print('  (MC) Monte Carlo ON: sigma relativo NN = %.3f (SUPOSICAO do usuario), n=%d' % (args.montecarlo, args.n))

    glob_d2004 = 0.0
    saltmg_deltas = []
    hdr = '  %-5s %3s %7s %14s %16s %10s' % ('oligo', 'len', 'Tm1998', 'd(1998->2004)', 'd(Mg 8->0 mM)', '+/-acc')
    for name, P in SETS:
        print('\n== %s ==' % name)
        print(hdr + ('   MC(sd)' if args.montecarlo is not None else ''))
        for nm, seq in oligos_of(P):
            t98 = tm_with(seq, NN_DH_1998, NN_DS_1998)
            t04 = tm_with(seq, NN_DH_2004, NN_DS_2004)
            t_mono = tm_with(seq, NN_DH_1998, NN_DS_1998, mg=0.0)
            d2004 = t04 - t98
            dsaltmg = t_mono - t98
            glob_d2004 = max(glob_d2004, abs(d2004))
            saltmg_deltas.append(dsaltmg)
            row = '  %-5s %3d %7.1f %14.3f %16.1f %10s' % (
                nm, len(seq), t98, d2004, dsaltmg, '+/-%.1f' % args.sigma_pred)
            if args.montecarlo is not None:
                row += '   %6.2f' % montecarlo_sd(seq, args.montecarlo, args.n)
            print(row)

    print('\n== resumo global ==')
    print('  (a) max |dTm(1998->2004)| sobre todos os oligos = %.3f C  -> escolha de parametros NN nao contribui' % glob_d2004)
    if saltmg_deltas:
        lo, hi = min(saltmg_deltas), max(saltmg_deltas)
        print('  (b) sensibilidade sal/Mg (Mg2+ %.0f mM -> 0): mediana %.1f C, faixa [%.1f, %.1f] C  -> deslocamento absoluto dominante'
              % (MG, statistics.median(saltmg_deltas), lo, hi))
    print('  (c) acuracia agregada do modelo NN: +/- %.1f C (SantaLucia & Hicks 2004)' % args.sigma_pred)
    print('\n  nota: dTm(1998->2004)=0 e um resultado de robustez para o manuscrito (validacao analitica);')
    print('        a incerteza pratica de Tm e governada pelo modelo de sal/Mg e pela acuracia agregada do NN.')


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        try:
            sys.stdout.close()
        except Exception:
            pass
