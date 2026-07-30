#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
concordance.py - Concordancia LAMPrime vs primers LAMP PUBLICADOS.
Para cada patogeno: localiza cada primer publicado na sequencia-alvo e calcula
Tm (modelo NN identico ao LAMPrime: SantaLucia 1998 + von Ahsen Mg->Na) + GC + posicao.
Reprodutivel e offline (sequencias-alvo lidas de data/; rede so se faltar arquivo). Saida: concordancia geometrica e termodinamica com ensaios LAMP publicados.
"""
import sys, math, os, urllib.request
sys.stdout.reconfigure(encoding='utf-8')

NN_DH={'AA':-7.9,'AT':-7.2,'TA':-7.2,'CA':-8.5,'GT':-8.4,'CT':-7.8,'GA':-8.2,'CG':-10.6,'GC':-9.8,'GG':-8.0,'AC':-8.4,'AG':-7.8,'TC':-8.2,'TG':-8.5,'TT':-7.9,'CC':-8.0}
NN_DS={'AA':-22.2,'AT':-20.4,'TA':-21.3,'CA':-22.7,'GT':-22.4,'CT':-21.0,'GA':-22.2,'CG':-27.2,'GC':-24.4,'GG':-19.9,'AC':-22.4,'AG':-21.0,'TC':-22.2,'TG':-22.7,'TT':-22.2,'CC':-19.9}
R=1.987; NA=50.0; MG=8.0; DNTP=1.4; OLIGO=50e-9
def free_mg_mM(tot_mM, dntp_mM, Ka=3e4):
    # Mg2+ livre pela quelacao de equilibrio 1:1 com dNTP (Ka=3e4 M^-1); MESMO modelo do
    # app.js tmNN (#6), no lugar da subtracao crua max(0, Mg-dNTP). Entrada/saida em mM.
    t=tot_mM/1000.0; d=dntp_mM/1000.0; B=Ka*d - Ka*t + 1.0
    return ((-(B)+math.sqrt(B*B+4*Ka*t))/(2*Ka))*1000.0
def comp(b): return {'A':'T','T':'A','G':'C','C':'G','N':'N'}.get(b,'N')
def rc(s): return ''.join(comp(c) for c in reversed(s))
def gc(s): 
    s=s.upper(); n=len(s) or 1; return 100.0*sum(c in 'GC' for c in s)/n
def tm(seq):
    seq=''.join(c for c in seq.upper() if c in 'ATGC'); N=len(seq)
    if N<2: return float('nan')
    dH=dS=0.0
    for i in range(N-1):
        d=seq[i:i+2]
        if d not in NN_DH: return float('nan')
        dH+=NN_DH[d]; dS+=NN_DS[d]
    for b in (seq[0],seq[-1]):
        h,s=(0.1,-2.8) if b in 'GC' else (2.3,4.1); dH+=h; dS+=s
    mgf=free_mg_mM(MG,DNTP); naeq=max(1e-3,(NA+120*math.sqrt(max(0.0,mgf))))/1000.0
    dS_salt=dS+0.368*(N-1)*math.log(naeq)
    return (dH*1000)/(dS_salt+R*math.log(OLIGO/4))-273.15

def find(target, primer):
    """retorna (start1based, end1based, fita) localizando primer (sense ou via revcomp)."""
    p=primer.upper()
    i=target.find(p)
    if i!=-1: return (i+1, i+len(p), '+')
    j=target.find(rc(p))
    if j!=-1: return (j+1, j+len(p), '-')
    return None

def split_fip(target, fip):
    """FIP=F1c+F2: F2 = maior sufixo 3' que é substring sense do alvo."""
    for L in range(min(28,len(fip)-10), 14, -1):
        f2=fip[-L:]
        if target.find(f2)!=-1:
            return fip[:-L], f2
    return None, fip
def split_bip(target, bip):
    """BIP=B1c+B2: B1c = maior prefixo 5' que é substring sense do alvo."""
    for L in range(min(28,len(bip)-10), 14, -1):
        b1c=bip[:L]
        if target.find(b1c)!=-1:
            return b1c, bip[L:]
    return bip, ''

def report(name, target, P):
    print(f"\n===== {name}  (alvo {len(target)} nt) =====")
    # componentes explicitos (F1c/F2/B1c/B2) tem prioridade; senao, split heuristico do FIP/BIP
    f1c,f2 = (P['F1c'],P['F2']) if P.get('F1c') and P.get('F2') else split_fip(target, P['FIP'])
    b1c,b2 = (P['B1c'],P['B2']) if P.get('B1c') and P.get('B2') else split_bip(target, P['BIP'])
    comps=[('F3',P['F3']),('F2',f2),('F1c',f1c),('B1c',b1c),('B2',b2),('B3',P['B3'])]
    if P.get('LF'): comps.append(('LF',P['LF']))
    if P.get('LB'): comps.append(('LB',P['LB']))
    print(f"{'prim':5} {'len':>3} {'Tm':>6} {'GC':>5}  pos(fita)        seq")
    locs={}
    for nm,seq in comps:
        if not seq: continue
        loc=find(target,seq); locs[nm]=loc
        ps=f"{loc[0]}-{loc[1]}({loc[2]})" if loc else "NAO ACHADO"
        print(f"{nm:5} {len(seq):>3} {tm(seq):6.1f} {gc(seq):5.0f}  {ps:16} {seq}")
    # amplicons
    if locs.get('F2') and locs.get('B2'):
        a=min(locs['F2'][0],locs['B2'][0]); b=max(locs['F2'][1],locs['B2'][1])
        print(f"  amplicon F2-B2 = {b-a+1} nt")
    if locs.get('F3') and locs.get('B3'):
        a=min(locs['F3'][0],locs['B3'][0]); b=max(locs['F3'][1],locs['B3'][1])
        print(f"  amplicon F3-B3 = {b-a+1} nt")

# ---- alvos de validacao (usados em report()/struct() abaixo) ----
# A. marginale: msp1b (Giglioti 2018); SARS-CoV-2: gene S (Prakash 2023);
# M. tuberculosis: IS6110 (Bentaleb et al. 2016).
# (O ensaio de msp5 NAO e usado para concordancia — primers de isolado divergente — removido.)

# Sequencias-alvo lidas de data/ (offline, deterministico). Fetch via rede so se faltar.
DATA=os.path.join(os.path.dirname(os.path.abspath(__file__)),'data')
def load_fasta(fname, fetch=None):
    p=os.path.join(DATA,fname)
    if os.path.exists(p):
        seq=''.join(l.strip() for l in open(p,encoding='utf-8') if not l.startswith('>'))
        return ''.join(c for c in seq.upper() if c in 'ATGCN')
    if fetch is not None:
        print('FASTA ausente (%s); baixando...'%fname)
        try: return fetch()
        except Exception as e: print('FALHA fetch:', e)
    return ''
def _fetch_sgene():
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text&seq_start=21563&seq_stop=25384'
    raw=urllib.request.urlopen(url, timeout=30).read().decode()
    seq=''.join(l.strip() for l in raw.splitlines() if not l.startswith('>'))
    return ''.join(c for c in seq.upper() if c in 'ATGCN')
def _fetch_mtb():
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=X17348.1&rettype=fasta&retmode=text'
    raw=urllib.request.urlopen(url, timeout=30).read().decode()
    seq=''.join(l.strip() for l in raw.splitlines() if not l.startswith('>'))
    return ''.join(c for c in seq.upper() if c in 'ATGCN')

# A. marginale msp1b: alvo sintetico (gBlocks) de Giglioti 2018 (base GenBank M59845.1).
# SARS-CoV-2 gene S: CDS de NC_045512.2 (21563-25384), Prakash 2023.
# M. tuberculosis IS6110: alvo GenBank X17348.1, Bentaleb et al. 2016 (BMC Infect Dis 16:517).
# (sequencias-alvo sao carregadas em main(); as definicoes/dicionarios sao importaveis por outros modulos.)

sarscov2={'F3':'TGGTGATATTGCTGCTAGA','B3':'GCACTATTAAATTGGTGGGC',
 'FIP':'AGGTCCAACCAGAAGTGATTCACCTTTGCTCACAGATG','BIP':'GCAGGTGCTGCATTACAATCTGTGTAACTCCAATACCA',
 'LF':'GCTAACAGTGCAGAAGTGTATT','LB':'GCTATGCAAATGGCTTATAGGT'}
amarginale_msp1b={'F3':'GCACTACCGTTCATGGATGA','B3':'TCCCCTGTGATATCTGTGCC',
 'FIP':'TGCCTTGCCAAATTCTTGCTCCCACCTGACACTGGTGAGAAG',
 'BIP':'AGCAGGCTTCAAGCGTACAGTTCCGCGAGCATGTGCA',
 'LF':'TCACCCGCTGGTACTTCAA','LB':'GCCTGGAGATGTTAGACCGA'}
mtb_is6110={'F3':'TCTCGTCCAGCGCCGCTT','B3':'GCGGGTCCAGATGGCTTG',
 'FIP':'ACGTAGGCGAACCCTGCCCCCAGCACCTAACCGGCTG',
 'BIP':'GTCACCGACGCCTACGCTCTCGCGTCGAGGACCATGG',
 'F1c':'ACGTAGGCGAACCCTGCCC','F2':'CCAGCACCTAACCGGCTG',
 'B1c':'GTCACCGACGCCTACGCTC','B2':'TCGCGTCGAGGACCATGG',
 'LF':'TCGACACATAGGTGAGGTC','LB':'TCGCTTCCACGATGGCCA'}

# ---- estruturas secundarias (hairpin/homo/heterodimero) @63C, modelo do LAMPrime ----
TREACT=63+273.15
def segdg(seg):
    seg=seg.upper(); N=len(seg)
    if N<2: return 0.0
    dH=dS=0.0
    for i in range(N-1):
        d=seg[i:i+2]
        if d in NN_DH: dH+=NN_DH[d]; dS+=NN_DS[d]
    for b in (seg[0],seg[-1]):
        h,s=(0.1,-2.8) if b in 'GC' else (2.3,4.1); dH+=h; dS+=s
    return dH-TREACT*dS/1000.0
def dimerdg(a,b):
    a=a.upper(); B=b.upper()[::-1]; best=0.0
    for off in range(-(len(a)-1), len(B)):
        i=max(0,-off)
        while i<len(a) and i+off<len(B):
            if i+off>=0 and a[i]==comp(B[i+off]):
                j=i; seg=''
                while j<len(a) and 0<=j+off<len(B) and a[j]==comp(B[j+off]): seg+=a[j]; j+=1
                if len(seg)>=3:
                    dg=segdg(seg);  best=min(best,dg)
                i=j
            else: i+=1
    return best
def hairpindg(a):
    a=a.upper(); best=0.0
    for i in range(len(a)):
        for j in range(len(a)-1, i+3, -1):
            if a[i]!=comp(a[j]): continue
            k=0; seg=''
            while i+k<j-k and (j-k)-(i+k)-1>=3 and a[i+k]==comp(a[j-k]): seg+=a[i+k]; k+=1
            if len(seg)>=3: best=min(best,segdg(seg))
    return best
def struct(name, P):
    prims=[('F3',P['F3']),('B3',P['B3']),('FIP',P['FIP']),('BIP',P['BIP'])]
    if P.get('LF'): prims.append(('LF',P['LF']))
    if P.get('LB'): prims.append(('LB',P['LB']))
    wh=min(hairpindg(s) for _,s in prims)
    wsd=min(dimerdg(s,s) for _,s in prims)
    wcd=0.0
    for x in range(len(prims)):
        for y in range(x+1,len(prims)):
            wcd=min(wcd, dimerdg(prims[x][1], prims[y][1]))
    print(f"  [estrutura @63C] pior hairpin ΔG {wh:.1f} | pior self-dimero {wsd:.1f} | pior hetero-dimero {wcd:.1f} kcal/mol (limiar -3.0)")

def main():
    # Sequencias-alvo lidas de data/ (offline, deterministico). Fetch via rede so se faltar arquivo.
    msp1b=load_fasta('amarginale_msp1b_synthetic.fasta')
    sgene=load_fasta('sarscov2_spike_NC045512.2_21563-25384.fasta', fetch=_fetch_sgene)
    if sgene: print('S gene:', len(sgene), 'nt')
    mtb=load_fasta('mtb_is6110_X17348.fasta', fetch=_fetch_mtb)
    if mtb: print('IS6110:', len(mtb), 'nt')
    report('A. marginale msp1b (Giglioti 2018, Exp Appl Acarol) - alvo gBlocks sintetico', msp1b, amarginale_msp1b)
    struct('A. marginale msp1b', amarginale_msp1b)
    if sgene:
        report('SARS-CoV-2 gene S (Prakash 2023, MethodsX)', sgene, sarscov2)
        struct('SARS-CoV-2 gene S', sarscov2)
    if mtb:
        report('M. tuberculosis IS6110 (Bentaleb et al. 2016, BMC Infect Dis) - alvo GenBank X17348.1', mtb, mtb_is6110)
        struct('M. tuberculosis IS6110', mtb_is6110)

if __name__=='__main__':
    main()
