# treeamps

Rust rewrite of a tree-level amplitude tensor-structure generator and gauge-constraint solver.

The original implementation lives in C++ under `tree_amplitudes`. This Rust workspace (`treeamps`) mirrors and extends the same logic, with a clearer separation between physics/combinatorics and CLI.

## Vision
Unify amplitude bootstrap steps into a single CLI that, given particle content and configuration constraints, outputs:
- Tensor basis (Lorentz-invariant monomials built from dot products).
- Gauge-variation constraint matrix (numeric or symbolic).
- Nullspace (gauge-invariant ansatz) with coefficients as rational functions of Mandelstam invariants.
- Canonical amplitude representatives after imposing kinematic relations (momentum conservation, on-shell, Mandelstam sum rules).

Planned extension: support generic mixtures (gluons, scalars; later gravitons), polarization patterns, optional color-factor placeholders, and upgrade to double-copy–friendly graviton basis.

## Symbolic Roadmap (In Progress)
1. Symbolic gauge matrix builder producing rows in Q[s,t,u] (later reduced to Q[s,t]).
2. Momentum conservation + Mandelstam substitution (u = -s - t for 4-point; general s_{ij}).
3. Full automatic generation of all 4-gluon (min 1 metric contraction) gauge equations (53 independent rows) → nullspace dim = 1.
4. Generalization to higher n: selective elimination of dependent scalar products, sparse polynomial representation, rational function simplification.
5. CLI normalization options (fix coefficient index or impose sum rule) and pretty-print rational expressions.
6. Verification harness comparing nullspace dimensions to thesis / literature tables.
7. Optional export: JSON (machine), LaTeX (report), and symbolic code (e.g. for substitution into other CAS).

## Unified CLI Goal
Introduce a `solve` command accepting:
```
--legs 4 \
--types g,g,g,g          # particle types per leg (g=gluon, s=scalar, h=graviton future)\
--degree 2               # tensor structure factor count\
--ee-min 1 --ee-max 1    # EE contraction range\
--pol-pattern one-per-leg|unconstrained \
--gauge 1 --gauge 2 ...  # gauge legs\
--symbolic               # use symbolic algebra\
--mandelstam             # perform invariant reduction (e.g. u = -s - t)\
--normalize 27           # fix coefficient α_27 = 1\
--emit solution|basis|matrix|all
```
This supersedes example-specific subcommands once stable.

## Conceptual overview

### Objects

- **LegIndex**: 1-based external leg label `i = 1..n`.
- **Momenta** `p_i` and **polarizations** `e_i` per leg.
- **ScalarFactor**:
  - `PP`: `(p_i · p_j)` (pure momentum).
  - `PE`: `(p_i · e_j)` (one momentum, one polarization).
  - `EE`: `(e_i · e_j)` (pure polarization).
- **TensorStructure**:
  - Product of dot products, e.g. `(p1·p2) (p1·e3) (e2·e4)`.
  - Represents a single monomial in a tensor basis for the amplitude.

### Generation logic

Given:

- `GenConfig = { n_legs, eliminate_leg, transversality, pol_pattern }`.
- `target_degree = D` (number of factors).
- `ee_contractions = K` (number of EE factors).

The generator:

1. Builds a catalog of allowed dot products:
   - `PP`: `(p_i·p_j)` with `1 ≤ i < j ≤ n`, excluding any factor involving `eliminate_leg`.
  - `PE`: `(p_i·e_j)` with `i ≠ eliminate_leg`; if `transversality = ForbidPiDotEi`, also `i ≠ j`.
   - `EE`: `(e_i·e_j)` with `1 ≤ i < j ≤ n`.
2. Runs a DFS over **multisets** of these factors of size `D`:
   - Factors are chosen with repetition in non-decreasing index order (no overcounting by permutation).
   - Tracks current degree, current EE count, and per-leg polarization counts.
3. Applies pruning:
   - If degree exceeds `D` or EE count exceeds `K` → prune.
  - If `pol_pattern = OnePerLeg`:
     - If any leg already has >1 polarization → prune.
     - If remaining factors can no longer give each leg exactly 1 polarization → prune.
4. At leaves (`deg == D`):
   - Accept structures with `ee_contractions == K` and, when enforced, exactly one polarization per leg.
   - Canonicalize (sort factors) and deduplicate.

From a physics point of view, the generator enumerates all Lorentz tensors built from scalar products of `p_i` and `e_i` that obey simple local rules (no `p_elim`, optional transversality, optional one-polarization-per-leg).

### Gauge constraints and the `solve` command

The program implements **gauge-invariance constraints** via Ward identities. For each gauge leg `r`, we impose:

\[
A(\epsilon_r \to p_r) = 0
\]

Physically: replacing a polarization vector with its momentum (a gauge transformation) must leave the amplitude unchanged.

**Implementation:**

For each tensor structure $T_i$ in a basis, applying $\epsilon_r \to p_r$ produces a linear combination of other structures:

- $(p_j \cdot \epsilon_r) \to (p_j \cdot p_r)$ (changes `PE` to `PP`)
- $(\epsilon_i \cdot \epsilon_r) \to (\epsilon_i \cdot p_r) = (p_r \cdot \epsilon_i)$ (changes `EE` to `PE`)
- $(\epsilon_r \cdot \epsilon_r) \to (p_r \cdot p_r) = 0$ (on-shell massless particle)

This builds a constraint matrix $M$ where each row corresponds to one gauge variation on one basis structure. The **nullspace** of $M$ gives all gauge-invariant amplitude candidates.

The `solve` command:

1. Builds a union basis over degrees `deg ∈ [ceil(n/2), max_deg]` with `ee = n - deg` (thesis convention).
2. Constructs the gauge constraint matrix for specified gauge legs.
3. Computes the nullspace via SVD.
4. Reports: matrix dimensions, rank, nullspace dimension, and a sample solution.

**Important limitation**: The current implementation treats all scalar products as independent variables. For quantitative agreement with physics results, you need to:
- Implement momentum conservation: $\sum_{i=1}^n p_i = 0$
- Apply Mandelstam variable relations: $(p_i \cdot p_j) = \frac{1}{2}s_{ij}$ with dependencies
- Handle on-shell conditions: $p_i^2 = 0$

Without these, the nullspace dimensions will not match theoretical predictions. The gauge variation logic is correct, but the linear system needs additional physics constraints.

## Workspace layout

- `Cargo.toml` (workspace root): defines a Cargo workspace with members:
  - `treeamps-core`: library crate with physics/combinatorics and linear algebra.
  - `treeamps-cli`: binary crate with the command-line interface.

- `treeamps-core/src`:
  - `types.rs`: `LegIndex`, `ScalarKind`, `Transversality`, `PolarizationPattern`.
  - `dot_product.rs`: `ScalarFactor` struct, ordering, string formatting.
  - `tensor_structure.rs`: `TensorStructure` representation and canonicalization.
  - `generator.rs`: `GenConfig`, `CatalogCounts`, `count_valid_factors`, `generate_tensor_structures`.
  - `constraints.rs`: `GaugeVariation`, `apply_gauge_variation`, `build_gauge_matrix`.
  - `linear_system.rs`: `LinearSystem` wrapper with SVD-based `nullspace` method.

- `treeamps-cli/src/main.rs`:
  - CLI entrypoint using `clap`:
    - `demo`: Small 3-leg example.
    - `gen-ts`: Generate tensor structures for fixed degree and EE count.
    - `solve`: Find gauge-invariant structures via nullspace computation.

## How to build and run

From the workspace root `~/Projects/treeamps`:

```fish
cargo build
```

Run the small 3-leg demo:

```fish
cargo run -p treeamps-cli -- demo
```

Generate tensor structures for chosen parameters, e.g. the canonical 4-leg case (unrestricted polarizations):

```fish
cargo run -p treeamps-cli -- gen-ts --n 4 --deg 3 --ee 1 --elim 4
```

This prints all monomials and performs a built-in sanity check.

To enforce exactly one polarization per leg, add the flag:

```fish
cargo run -p treeamps-cli -- gen-ts --n 4 --deg 3 --ee 1 --elim 4 --one-pol-per-leg
```

Solve for gauge-invariant structures by specifying gauge legs:

```fish
cargo run -p treeamps-cli -- solve --n 4 --gauge 1 --gauge 2 --gauge 3
```

This builds a basis over degrees `[2, 4]`, constructs the gauge constraint matrix, and reports the nullspace dimension (number of independent gauge-invariant structures).

## Built-in physics sanity check

For `n = 4`, `deg = 3`, `ee = 1`, `elim = 4`, the generator uses the counts

- `PP = num_pp`, `PE = num_pe`, `EE = num_ee`

and a combinatorial formula

\[
\text{expected").append("_count} = E \\cdot \\binom{M+1}{2}, \\
M = \text{num_pp} + \text{num_pe},\\quad E = \text{num_ee}.
\]

Interpretation:

- Each structure has exactly one `EE` factor and two more factors from the multiset `{PP ∪ PE}`.
- There are `E` choices for the `EE` factor, and `C(M+1, 2)` multisets of size 2 from `M` other factors.

In the **unrestricted polarization pattern** (default, no extra CLI flag), the Rust program prints, for this canonical case:

- `PP=3, PE=9, EE=6  => expected count=468  (OK)`

matching the explicit enumeration of the 468 generated tensor structures.

In the **one-pol-per-leg regime** (use `--one-pol-per-leg`), the same parameter point instead yields 30 allowed tensor structures. The CLI reports this via a separate sanity line:

- `[Sanity-one-pol-per-leg] expected count=30  (OK)`

so you have an explicit check for both regimes.

## How to check if the program is behaving physically

You don’t need to read the Rust; treat `treeamps-cli` as a black box and use it like a calculator:

1. **Check basic pruning rules**
   - Use `demo` and `gen-ts` and inspect outputs:
     - No factors should involve the eliminated momentum leg.
     - With transversality enforced, there should be no `(p_i·e_i)` terms.
     - With one-polarization-per-leg (once exposed via config/flags), each leg should appear exactly once among all E’s.

2. **Check scaling and patterns**
   - Vary `n`, `deg`, `ee` in `gen-ts` and see if counts move as expected:
     - Increasing `n` increases the number of possible tensor structures.
     - Enforcing stronger constraints (transversality, one-polarization-per-leg) reduces counts.

3. **Compare against independent small-`n` calculations**
   - For a few special cases (like the 4-leg `deg=3`, `ee=1` example), derive your own counts and compare them with `gen-ts`.
   - Discrepancies indicate a bug in the combinatorics/pruning logic.

4. **Gauge invariance via `solve`**
   - Use the `solve` subcommand to find gauge-invariant structures:
     ```fish
     cargo run -p treeamps-cli -- solve --n 4 --gauge 1 --gauge 2 --gauge 3
     ```
   - Inspect the nullspace dimension (number of independent gauge-invariant structures).
   - Compare this against physics expectations for given `n` and degree range.
   - For 4-leg Yang-Mills amplitudes with all legs gauged, theory predicts specific nullspace dimensions you can check against.

If you know additional parameter points where counts or nullspace dimensions are fixed by theory, they can be added as explicit checks or tests to guard against regressions.
