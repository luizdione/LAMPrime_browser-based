# Changelog

## v2.0 — 2026-07-30
- Reproducibility test suite (`tools/tests/`): 14 deterministic checks that verify the
  nearest-neighbor constants and the salt/Mg²⁺ equations against the primary literature
  (SantaLucia 1998; Allawi & SantaLucia 1997; Peyret 1999; Owczarzy 2008; von Ahsen 2001)
  and against BioPython 1.86 as an independent transcription. Added a GitHub Actions CI
  workflow (`.github/workflows/tests.yml`) that runs the suite, `node --check app.js` and
  the engine self-test on every push.

## Housekeeping — 2026-07-17
- Reorganização do repositório (sem mudança de comportamento do app): removida a cópia
  obsoleta de `concordance.py` na raiz (a versão canônica é `tools/concordance.py`);
  referências de HTML/README e a allowlist do `.gitignore` atualizadas; `tools/` mantido
  intacto.

## v1.4.1 — 2026-07-16
- Consistência do motor: `tools/crossreact.js` passa a usar o mesmo Mg²⁺ livre por quelação (Ka=3×10⁴) de `app.js` (#6) e torna-se importável. Como o conjunto top-ranked (Set #1) é sensível ao empate de penalidade entre desenhos quase equivalentes, a triagem de reatividade cruzada passa a ser reportada em **ensemble** sobre os 10 melhores desenhos (`ensembleCrossReact`): *A. marginale* msp1b × *A. centrale* dá **0–0/8** (nenhum primer sinaliza em nenhum desenho); *B. bovis* 18S × *B. bigemina* dá **2–7/8** (mediana 3), com **todos os 8 primers** sinalizando em ao menos um desenho — a reatividade no locus 18S conservado é, portanto, essencialmente inevitável, um resultado robusto ao desempate.
- Novo `tools/pfalciparum_design.js`: desenho LAMP do alvo AT-rico *P. falciparum* 18S pelo motor do app (penalidade 1,2; F2–B2 175 nt); demonstra que o desafio de baixo GC recai sobre os primers de alça.

## v1.4 — 2026-07-16
- Motor de Tm: unificação do Mg²⁺ livre. O motor de desenho (`tmNN`) passa a computar o Mg²⁺ livre pelo mesmo equilíbrio de quelação 1:1 com dNTP (Ka=3×10⁴ M⁻¹) já usado pela varredura de especificidade, no lugar da subtração crua `max(0, Mg−dNTP)`. No ponto de operação padrão (Mg 8 mM, dNTP 1,4 mM) o efeito é desprezível (Mg livre 6,600→6,607 mM); a correção importa perto da estequiometria (Mg≈dNTP), onde o modelo cru zerava o Mg livre e derrubava o Tm. Espelhado em `tools/concordance.py` para manter a concordância. Ferramentas novas de validação: `tools/tm_uncertainty.py` (incerteza de Tm; SantaLucia 1998 vs 2004), `tools/benchmark_independent.py` (verificação independente do Tm vs Biopython; acordo ≤0,09 °C) e `tools/specificity_metrics.py` (sens/spec/ROC parcial da triagem de especificidade). Novo alvo AT-rico *P. falciparum* 18S (`tools/data/pfalciparum_18s_M19172.fasta`). `concordance.py` importável; correção do ano de Giglioti (2019→2018).

## v1.3 — 2026-06-24
- Specificity Mg sweep: heatmap now colors by the Mg-dependent off-target safety margin (T_rxn − Tm_off) instead of the Mg-invariant discrimination margin S; suggested Mg keeps off-targets below the reaction temperature while on-targets stay functional; clarified that Mg shifts absolute Tm, not the match/mismatch margin.

## v1.2 — 2026-06-23
- Specificity tab: Mg²⁺ sweep (2–10 mM) producing a per-primer on-target−off-target Tm-margin heatmap (Owczarzy 2008 divalent correction; Allawi/Peyret internal-mismatch NN); suggested working [Mg²⁺] by maximin. New offline validation script tools/mgspec.mjs.

## v1.1 — 2026-06-23
- Specificity tab: added one-click cross-reactivity example presets (A. marginale vs A. centrale CP001759.1; B. bovis 18S vs B. bigemina KP710228.1).
- Documented the in-browser specificity screen used for the published cross-reactivity benchmark (tools/crossreact.js).
