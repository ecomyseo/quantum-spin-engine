# Changelog

## 0.1.1 — 2026-10-06

- `bloch-renderer`: on a short canvas (a phone held sideways) the sphere now shrinks and recentres in
  the room left above its caption, so the south-pole label no longer overlaps the `|⟨σ⟩|` readout and
  the north-pole label is never clipped at the top.

## 0.1.0 — 2026-10-05

First public version, extracted from the spin, Stern–Gerlach and one-qubit scenarios of
[Fundamental Universe](https://fundamental.gmartos.es/).

- `cmatrix`: small dense complex matrices, matrix exponential, Hermitian eigenvalues.
- `density`: kets, density matrices, Bloch vector, SU(2) rotations, purity, entropy, fidelity,
  trace distance, partial trace.
- `open-system`: Lindblad superoperator with an exact propagator, qubit Kraus channels and
  quantum-jump trajectories.
- `measurement`: projective measurement along any axis, shot statistics, one-qubit tomography.
- `spinor`: spin-½ wavefunction on a 1D grid for Stern–Gerlach magnets.
- `bloch-renderer`: Canvas 2D view with the Bloch sphere, a beam map, the analyser screen and a
  stage strip.
- Live demo, three runnable examples and a test suite that checks every piece against an exact law
  or an independent method.
