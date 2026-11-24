# Gate 001 — Known-library identifiability smoke test

This is an initial simulated calibration of an instrument, not a new reservoir, a blind system-identification result, or a claim of novelty. The finite hypothesis library is deliberately supplied to the instrument. A physical circuit is the next external validation.

- 36 candidate models: two first-pole values, two second-pole values, linear or cubic nonlinearity (two strengths) at the input, between filters, or at output. Output plant and sensor alternatives are exact input/output twins. The two low-pass filters have unity DC gain. All runs start from zero state.
- 10 candidate inputs, each 80 samples, equal RMS 0.60 and peak <= 1.2. Three experiments per arm; independent Gaussian output noise SD 0.12; same noise draws across arms per trial. 144 deterministic trials, four per candidate, with a mixed-state seeded PRNG. Known noise and a uniform candidate prior; no process noise or unknown model class.
- Arms: fixed constant, fixed prior-EIG ranked sequence, random, adaptive posterior-weighted pair separation, and standard Monte Carlo expected information gain (EIG), targeting observable groups.
- Posterior mass >= 0.95 is required for a group call. Report correct calls, wrong calls, abstentions, group negative log probability, and the per-round ledger. Do not claim plant vs sensor localization for the exact twins. Equal posterior weights must stay equal for every tested trace.
- The finite-library group decision is conditional on the supplied models and noise assumption, not a universal confidence certificate.
- A future gate must introduce out-of-library models, model discrepancy, uncertain noise, independently measured hardware and conventional system-ID baselines. A failure to outperform fixed or standard EIG is a negative result for adaptive pair selection, not a reason to retune this run.

Protocol and results were generated together for this initial smoke test; this is not a prospective preregistration. Freeze this file for subsequent replication.

Audit note: the initial 120-trial run used poorly mixed consecutive seeds; random first probes were concentrated in two choices. It was superseded by this balanced smoke test; initial commits are preserved in history. Do not call either run prospective.
