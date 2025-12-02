# AI Coding Agent Instructions for treeamps

Purpose: Make AI agents immediately productive in this Rust workspace that generates tensor-structure bases and gauge-constraint matrices for tree-level scattering amplitudes.

## Big Picture
- Workspace: two crates
  - `treeamps-core` (library): physics/combinatorics + linear algebra
  - `treeamps-cli` (binary): command-line interface
- Flow:
  1. Enumerate tensor structures from dot products of momenta `p_i` and polarizations `e_i` under simple local rules.
  2. Build gauge-variation constraint rows (Ward identities) into a matrix.
  3. Compute nullspace to get gauge-invariant amplitude candidates.
- Current limitation: scalar products are treated as independent. Physics-consistent results require momentum conservation, on-shell, and Mandelstam variable relations; partial symbolic scaffolding exists and will evolve.

## Key Modules (treeamps-core/src)
- `types.rs`: `LegIndex`, `ScalarKind`, `Transversality`, `PolarizationPattern` enums and type aliases.
- `dot_product.rs`: `ScalarFactor` variants (`PP`, `PE`, `EE`), ordering, display.
- `tensor_structure.rs`: `TensorStructure` representation (multiset of factors), canonicalization.
- `generator.rs`: `GenConfig`, `CatalogCounts`, `count_valid_factors`, `generate_tensor_structures` (DFS with pruning by degree, EE count, and polarization pattern).
- `constraints.rs`: `GaugeVariation`, `apply_gauge_variation`, `build_gauge_matrix` (maps `e_r→p_r` rules over basis to linear combinations).
- `linear_system.rs`: `LinearSystem` wrapper using `nalgebra` SVD to compute `nullspace`.
- Symbolic work-in-progress lives in `symbolic.rs` and `symbolic_constraints.rs` for Mandelstam and momentum-conservation substitutions.

## CLI (treeamps-cli/src/main.rs)
- Subcommands:
  - `demo`: sanity example (3 legs)
  - `gen-ts`: generate tensor structures for fixed `--n`, `--deg`, `--ee`, `--elim` and optional `--one-pol-per-leg`
  - `solve`: build basis over degrees, construct gauge matrix for `--gauge <leg>` flags, compute and report nullspace
- Typical runs (fish shell):
  - `cargo run -p treeamps-cli -- demo`
  - `cargo run -p treeamps-cli -- gen-ts --n 4 --deg 3 --ee 1 --elim 4`
  - `cargo run -p treeamps-cli -- gen-ts --n 4 --deg 3 --ee 1 --elim 4 --one-pol-per-leg`
  - `cargo run -p treeamps-cli -- solve --n 4 --gauge 1 --gauge 2 --gauge 3`

## Build, Test, Debug
- Build: `cargo build`
- Run: see CLI examples above; prefer `-p treeamps-cli` from workspace root.
- There is no test suite yet; use CLI outputs and the printed sanity checks (expected counts) to validate behavior.
- If you modify `linear_system` or degree ranges, verify nullspace dimensions don’t regress using `solve` for `n=4`.

## Project Conventions
- Canonicalization: tensors store factors in sorted, stable order to avoid overcounting by permutation; generators select factors in non-decreasing index order.
- Polarization pattern:
  - Default: unconstrained
  - `OnePerLeg`: pruning enforces exactly one `e_i` appearance per leg
- Eliminated leg: `GenConfig.eliminate_leg` forbids any factor involving that momentum (e.g., `p_4` in 4-leg examples).
- Gauge variation rules:
  - `(p_j·e_r) → (p_j·p_r)` (PE→PP)
  - `(e_i·e_r) → (e_i·p_r)` (EE→PE)
  - `(e_r·e_r) → (p_r·p_r)=0` (on-shell)
- Linear algebra: `nalgebra` SVD on dense matrices; matrices are currently numeric over independent scalar products. Symbolic variants are being introduced but are not authoritative yet.

## Integration Points and Extensibility
- External crates: `nalgebra` (SVD), `num-{bigint,rational,traits}` for symbolic/scalar work.
- CLI uses `clap` with derive; extend subcommands in `treeamps-cli/src/main.rs` and wire to core APIs.
- Planned symbolic features in `symbolic.rs`/`symbolic_constraints.rs`: Mandelstam substitution (`u = -s - t` for 4-point), momentum conservation, sparse polynomials. Keep numeric and symbolic paths separate behind clear types.

## Practical Guidance for Agents
- When adding generator features, update `CatalogCounts` and pruning in `generator.rs`, plus sanity prints in CLI.
- When changing basis degree ranges for `solve`, keep thesis convention: union over `deg ∈ [ceil(n/2), max_deg]` with `ee = n - deg`.
- For new physics constraints, introduce transformations in symbolic modules and, if needed, a separate `LinearSystem` variant; don’t bake physics into factor enumeration.
- Prefer small, targeted changes; maintain ordering and dedup rules to avoid exponential blowups.
- Validate with `gen-ts` counts for `n=4, deg=3, ee=1, elim=4` (expect 468 unconstrained; 30 with `--one-pol-per-leg`).

## References in Repo
- Overview and workflows: `README.md`
- Workspace membership: root `Cargo.toml`
- Core APIs: `treeamps-core/src/*.rs`
- CLI entry: `treeamps-cli/src/main.rs`
