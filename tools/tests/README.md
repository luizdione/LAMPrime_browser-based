# LAMPrime reproducibility test suite

Independent, runnable verification that the LAMPrime thermodynamic engine
(`app.js`, mirrored offline by `tools/mgspec.mjs` and `tools/concordance.py`)
implements the published constants and equations **faithfully** — every value
checked against the primary literature and against BioPython as an independent
transcription. Nothing is asserted against our own code alone.

## Run

```bash
pip install -r tools/tests/requirements.txt   # biopython==1.86
python tools/tests/run_tests.py               # 14 checks, no pytest needed
# (pytest tools/tests also works)
```

Requires Python 3.10+ and Node.js 18+ (the engine math is imported from
`tools/mgspec.mjs`, the importable mirror of `app.js`).

## What each test proves

| File | Checks | Primary source |
|------|--------|----------------|
| `test_nn_constants.py` | Watson–Crick NN ΔH/ΔS == golden **and** BioPython `DNA_NN3`; Python (`concordance.py`) ≡ JS (`mgspec.mjs`) | SantaLucia 1998, *PNAS* 95:1460, Table 1 |
| `test_mismatch_imm.py` | all 51 internal-mismatch NN entries == BioPython `DNA_IMM1`; all mismatch families present | Allawi & SantaLucia 1997; Peyret et al. 1999 |
| `test_owczarzy_freemg.py` | Mg²⁺-corrected Tm == a from-the-paper re-implementation (max\|Δ\| < 1e-9 °C); low-Mg regime; 1:1 Mg–dNTP chelation (Ka=3e4) | Owczarzy et al. 2008, *Biochemistry* 47:5336 |
| `test_vonahsen_tm.py` | S1.1 design-engine Tm (NN + salt + Mg→Na) reproduces the validation-panel Tm values (SARS-CoV-2 spike; M. tuberculosis IS6110) | von Ahsen 2001; SantaLucia 1998 |
| `test_biopython_parity.py` | Mg²⁺-path duplex Tm == BioPython `MeltingTemp.Tm_NN` (saltcorr=7): perfect panel max\|Δ\| ≈ 0.000 °C, single-mismatch max\|Δ\| = **0.31 °C** | BioPython 1.86 (independent Tm) |

`golden_literature.py` holds the values typed from the papers plus helpers;
`_dump_constants.mjs` / `_dump_owcz.mjs` expose the engine's runtime constants to
the Python checks.
