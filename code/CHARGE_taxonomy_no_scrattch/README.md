# CHARGE taxonomy without scrattch runtime dependencies

This codebase is derived from the user-supplied CHARGE taxonomy source. Functions added locally are tagged `## NEW FUNCTION`. The existing `plot_constellation` function is tagged `## EDITED FUNCTION` because its stale `scrattch.bigcat` import declaration was removed.

No `scrattch.taxonomy`, `scrattch.bigcat`, `scrattch.hicat`, `scrattch.io`, or `scrattch.mapping` package is imported or attached.

## Important validation note

The cluster sum, mean, squared-mean, and variance helpers use mathematically equivalent `Matrix` operations in place of the original compiled scrattch.bigcat routines. The local `bezier` function implements the quadratic three-control-point curve used by `edgeMaker`. Validate generated `CHARGE.RData` objects against a known reference before replacing a production workflow.
