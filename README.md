# General Instrument

**A machine for investigating machines.** Started 24 November 2026, following the GraphNN control audit and the Theorem operating-point microscope discussion.

A scientific instrument that proposes competing explanations, chooses discriminating probes, measures responses, updates beliefs, and reports unresolved ambiguity. Not a new neural architecture.

## Status

- **Gate 001 executed:** 144 balanced simulated trials per arm, 36 supplied models, three equal-energy probes. [Protocol](PROTOCOL.md) | [Code](gate001.mjs) | [Results](results/gate001.json).
- **Gate 002 scaffold:** conventional recorded-data prediction baselines for the measured Wiener–Hammerstein and Silverbox benchmarks. The physical benchmark runs have not yet been executed. No claim of real-world identification or improved benchmark performance.

### Gate 001

| Policy | Correct group | Wrong confident group | Abstain |
|---|---:|---:|---:|
| Constant | 19 | 1 | 124 |
| Random | 57 | 2 | 85 |
| Fixed prior-EIG ranked | 120 | 1 | 23 |
| Adaptive pair separation | 139 | 0 | 5 |
| Standard Bayesian EIG | 139 | 0 | 5 |

The plant-output and sensor-output nonlinearities are exactly equivalent in this setup. No input-only measurement separates them. The reported group intentionally does not resolve them; this is a known-library positive control, not discovery of equivalence in arbitrary models. Pair separation ties standard EIG; no novelty or general superiority established.

Run with Node.js 18+ (no dependencies):

```sh
node gate001.mjs
```

### Gate 002: measured systems, ordinary baselines first

```sh
python -m pip install numpy nonlinear-benchmarks
python benchmark.py smoke
python benchmark.py wh --output results/wh.json
python benchmark.py silverbox --output results/silverbox.json
```

The benchmark loader uses the published splits; first use may download the datasets. The runner fits a fixed 50-tap FIR, polynomial FIR, and a two-stage Wiener-style FIR using training data only; it reports held-out RMSE in native units after the specified warm-up. These are preliminary reference models, not full WH identification, a strong nonlinear state-space baseline, or an active probing test. If the loader's default download fails, consult the official dataset page.

## Scientific boundaries

1. Prediction error alone cannot establish where an internal component lives.
2. If competing mechanisms yield identical input-output maps, a different sensor, internal measurement or intervention is required.
3. A probability conditional on a small model library is not evidence that reality belongs to that library.
4. Never use benchmark test outputs to tune a model. Recorded data cannot validate a policy's choice of new physical interventions.
5. Next: out-of-library rejection and unknown-noise gates, then fit established nonlinear system identification methods on recorded hardware data, then evaluate a controllable physical or simulated testbed under a fixed intervention budget.

## Lineage

[Luotain2](https://github.com/anttiluode/Luotain2) — candidate laws and noisy measurements; [Theorem](https://github.com/anttiluode/Theorem) — distinguishability through interventions; [WonderingQwenSol](https://github.com/anttiluode/WonderingQwenSol) — informative observations; [BlackBoxSpotlight](https://github.com/anttiluode/BlackBoxSpotlight) — internal operators; [GraphNNSol](https://github.com/anttiluode/GraphNNSol) — ordinary baseline discipline.

## Prior art and data

- [Wiener–Hammerstein benchmark](https://www.nonlinearbenchmark.org/benchmarks/wiener-hammerstein)
- [Silverbox benchmark](https://www.nonlinearbenchmark.org/benchmarks/silverbox)
- [Official data loader](https://github.com/MaartenSchoukens/nonlinear_benchmarks)
- Schoukens & Tiels, *Identification of block-oriented nonlinear systems starting from linear approximations: A survey*, Automatica 85 (2017).
- [Parametric identification of parallel Wiener–Hammerstein systems](https://doi.org/10.1016/j.automatica.2014.10.105) — BLA across excitation levels and structural ambiguity.

Standard system identification and optimal experimental design are the relevant prior art. The project must earn any additional contribution through independent controls.
