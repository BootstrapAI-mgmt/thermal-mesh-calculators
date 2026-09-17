# Copilot Instructions — thermal-mesh-calculators

## Project Summary

Pure Python thermal mesh sizing calculators for automotive CAE. Derives maximum finite element sizes from heat transfer physics (Fourier's law, Stefan-Boltzmann, Biot number, Newton-Raphson energy balance solvers).

## Critical Constraints

- **NO external dependencies.** Do not import numpy, scipy, pandas, or any non-stdlib package. The code must run in CAE pre-processor embedded Python (HyperMesh, ANSA).
- **All temperatures in Kelvin** internally. Output dicts may include `_C` convenience keys.
- **All lengths in meters** internally. Output mesh sizes are in `_mm` keys.
- **Dict returns** (not dataclasses) for Python 3.6+ compatibility.

## Code Style

- Type hints on all public methods
- NumPy-style docstrings (Parameters / Returns sections)
- Module docstrings explain the governing physics and equations
- Static methods preferred (no persistent state)
- Snake_case everywhere

## Testing

```bash
python -m pytest tests/ -q              # 239 passed
python -m examples.automotive_examples  # 11 worked scenarios
```

239 tests across 11 files: **221 physics** (hand-computed analytical solutions,
energy-balance closure on the Newton-Raphson shield solvers, scaling checks) plus
**18** for the open/closed map checker, which is governance rather than physics.
New physics tests must be verified against an independently derived analytical
result, never against the code's own current output.

Do not propagate a test count without checking which pin you mean. `v0.6.0` has
10 test files and runs **221**; `v0.6.1` has 11 and runs
**239**. Downstream repos that still pin `v0.6.0` are correct to say 221.
The wheel ships no tests at all — only the sdist does (`MANIFEST.in`).

## Key Physics

The conduction flux q'' is unknown before solving, so we substitute the surface energy balance:
`dx_max = k * dT_max / |h(Ts-Tf) + eps*sigma*(Ts^4 - Tsurr^4)|`

Heat shields use Newton-Raphson to solve for their equilibrium temperature before computing mesh sizes. Single-layer uses scalar N-R; multilayer uses 2x2 Jacobian with Cramer's rule.
