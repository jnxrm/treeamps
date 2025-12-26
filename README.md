# treeamps

**Tensor-structure generator for tree-level scattering amplitudes in quantum field theory**

A high-performance Rust implementation that systematically enumerates Lorentz-invariant basis tensors built from momenta and polarization vectors, with built-in validation and benchmarking capabilities.

## Quick Start

```fish
# Build the project
cargo build --release

# Generate tensor structures for 4 gluons (canonical example)
cargo run -p treeamps-cli --release -- gen-ts --n 4 --deg 3 --ee 1

# Run benchmarks
cargo bench

# Profile with flamegraph (requires cargo-flamegraph)
cargo flamegraph -p treeamps-cli -- gen-ts --n 9 --ee 2
```

## What This Project Does

Given a configuration of external particles (currently gluons), this tool systematically generates all possible Lorentz-invariant tensor structures that can appear in tree-level scattering amplitudes.

### Core Physics Concepts

- **Tensor Structures**: Products of dot products like `(p₁·p₂)(p₁·e₃)(e₂·e₄)`, where:
  - `pᵢ` are external momenta (4-vectors)
  - `eᵢ` are polarization vectors for gauge bosons
  
- **Scalar Factors**: The building blocks of tensor structures:
  - `PP`: Momentum-momentum dot products `(pᵢ·pⱼ)`
  - `PE`: Momentum-polarization dot products `(pᵢ·eⱼ)`
  - `EE`: Polarization-polarization dot products `(eᵢ·eⱼ)`

- **Gluon Basis**: For `n` gluons with "one polarization per leg" constraint:
  - Each tensor structure has exactly `n` total factors
  - Constraint: `2×EE + PE = n` (each polarization appears exactly once)
  - This implies: `deg = n - ee` where `deg` is total degree

### Why This Matters

In quantum field theory, scattering amplitudes can be decomposed into a basis of Lorentz tensors multiplied by scalar coefficients. This tool automates the enumeration of valid basis elements, which is combinatorially challenging for higher particle multiplicities.

## Project Status (December 2025)

### ✅ Fully Implemented
- **Tensor structure generation** with efficient DFS enumeration
- **One-polarization-per-leg constraint** for gluon bases
- **Transversality enforcement** (`pᵢ·eᵢ = 0`)
- **Momentum elimination** (last leg by convention)
- **Canonical ordering** to avoid permutation redundancy
- **Built-in validation** with known combinatorial formulas
- **Criterion benchmarks** for performance tracking
- **Flamegraph profiling** support

### 🚧 Previously Implemented, Now Removed
The following features were part of earlier development but have been **intentionally removed** to streamline the codebase:
- Gauge constraint matrix construction
- Nullspace computation via SVD
- Symbolic algebra scaffolding
- `solve` subcommand

These may be re-introduced in the future with improved architecture.

### 📋 Future Roadmap
1. **Symbolic algebra integration** with momentum conservation and Mandelstam relations
2. **Gauge constraint solver** rebuilt with better symbolic handling
3. **Extended particle types**: scalars, gravitons, mixed configurations
4. **Double-copy methods** for graviton amplitudes
5. **Export capabilities**: JSON, LaTeX, Computer Algebra System formats

## Architecture Overview

### Workspace Structure

```
treeamps/
├── Cargo.toml              # Workspace definition
├── treeamps-core/          # Core library crate
│   ├── src/
│   │   ├── lib.rs          # Public API
│   │   ├── types.rs        # Core types (LegIndex, ScalarKind, etc.)
│   │   ├── dot_product.rs  # ScalarFactor representation
│   │   ├── tensor_structure.rs  # TensorStructure type
│   │   └── generator.rs    # Main generation algorithm
│   ├── benches/
│   │   └── generator.rs    # Criterion benchmarks
│   └── Cargo.toml
└── treeamps-cli/           # Binary crate
    ├── src/
    │   └── main.rs         # CLI with clap
    └── Cargo.toml
```

### Key Types (from `types.rs`)

```rust
pub struct LegIndex(pub u8);  // 1-based external leg labels

pub enum ScalarKind {
    PP,  // Momentum × Momentum
    PE,  // Momentum × Polarization  
    EE,  // Polarization × Polarization
}

pub enum Transversality {
    Allow,          // Allow all PE factors
    ForbidPiDotEi,  // Enforce pᵢ·eᵢ = 0 (physical transversality)
}

pub enum PolarizationPattern {
    Unconstrained,  // Any number of polarizations per leg
    OnePerLeg,      // Exactly one polarization per leg (gluon basis)
}
```

### ScalarFactor (from `dot_product.rs`)

Represents a single dot product with canonical ordering:

```rust
pub struct ScalarFactor {
    pub kind: ScalarKind,
    pub a: LegIndex,
    pub b: LegIndex,
}
```

Examples:
- `ScalarFactor::pp(1, 2)` → `(p₁·p₂)`
- `ScalarFactor::pe(1, 3)` → `(p₁·e₃)`
- `ScalarFactor::ee(2, 4)` → `(e₂·e₄)`

Always maintains `a < b` for PP and EE, or `a` is the momentum leg for PE.

### TensorStructure (from `tensor_structure.rs`)

Represents a complete tensor monomial as an **ordered multiset** of ScalarFactors:

```rust
pub struct TensorStructure {
    pub factors: Vec<ScalarFactor>,
    pub ee_contractions: u32,
}
```

Example: `(p₁·p₂)(p₁·e₃)(e₂·e₄)` has 3 factors, 1 EE contraction.

### GenConfig (from `generator.rs`)

Configuration for generation:

```rust
pub struct GenConfig {
    pub n_legs: u8,
    pub transversality: Transversality,
    pub pol_pattern: PolarizationPattern,
}
```

Convention: The last leg (`n_legs`) is always the **eliminated leg** — its momentum `pₙ` never appears in any factor due to momentum conservation.

## Generation Algorithm

The core function `generate_tensor_structures(cfg, deg, ee)` uses **depth-first search with aggressive pruning**:

### Phase 1: Build Factor Catalog

For `n` legs, generate all allowed factors:

1. **PP factors**: `(pᵢ·pⱼ)` with `1 ≤ i < j < n` (excludes last leg)
2. **PE factors**: `(pᵢ·eⱼ)` where:
   - `i < n` (no eliminated momentum)
   - `i ≠ j` if transversality enforced
3. **EE factors**: `(eᵢ·eⱼ)` with `1 ≤ i < j ≤ n`

### Phase 2: DFS with Pruning

Enumerate multisets of exactly `deg` factors with `ee` EE contractions:

```
dfs(current_structure, next_factor_index):
    // Pruning conditions
    if current_degree > target_deg: return
    if current_ee_count > target_ee: return
    
    if pol_pattern == OnePerLeg:
        if any leg has > 1 polarization: return
        if impossible to reach 1 per leg with remaining factors: return
    
    // Accept completed structures
    if current_degree == target_deg && current_ee_count == target_ee:
        if pol_pattern == OnePerLeg:
            if all legs have exactly 1 polarization:
                emit(canonicalize(current_structure))
        else:
            emit(canonicalize(current_structure))
    
    // Recurse: try adding each valid factor
    for factor in factors[next_factor_index..]:
        add factor to current_structure
        dfs(current_structure, index_of(factor))
        remove factor from current_structure
```

**Key optimizations**:
- Non-decreasing factor selection (avoids permutation overcounting)
- Early pruning on degree and EE count
- Polarization counting for one-per-leg constraint
- BTreeSet deduplication for canonical results

### Phase 3: Validation

Built-in sanity checks for known combinatorial formulas:

**4-gluon mixed basis** (`n=4, deg=3, ee=1, one-per-leg`):
- Expected: 24 tensor structures
- Formula: specific counting based on factor catalog

**4-gluon pure-EE basis** (`n=4, deg=2, ee=2, one-per-leg`):
- Expected: 3 tensor structures  
- Formula: `C(4,2) - 3 = 3` (pairs minus on-shell forbidden)

## CLI Usage Guide

### Available Commands

Currently, the CLI provides one main command:

```bash
cargo run -p treeamps-cli -- gen-ts [OPTIONS]
```

### Options

- `--n <NUMBER>`: Number of external legs (default: 3)
  - Minimum: 1
  - Recommended practical maximum: ~10 for reasonable performance
  
- `--deg <DEGREE>`: Total number of factors in each tensor structure
  - Default: 0 (auto-computed from `--n` and `--ee`)
  - For gluons with one-per-leg: `deg = n - ee`
  
- `--ee <COUNT>`: Number of EE (polarization×polarization) contractions
  - Default: 0 (auto-computed from `--n` and `--deg`)
  - Must satisfy: `ee ≤ deg`
  - For gluons with one-per-leg: `ee = n - deg`

### Gluon Basis Conventions

The CLI **always enforces** the one-polarization-per-leg pattern, which means:

1. **Consistency requirement**: `deg + ee = n`
2. **Parameter flexibility**: You can specify any two of {`n`, `deg`, `ee`}, and the third is computed
3. **Default behavior**: If only `--n` is given, defaults to pure PE basis (`deg=n, ee=0`)

### Examples

**Example 1**: Pure PE basis (no EE contractions)
```fish
cargo run -p treeamps-cli -- gen-ts --n 4
# Equivalent to: --n 4 --deg 4 --ee 0
# Generates structures like (p₁·e₂)(p₁·e₃)(p₁·e₄)(p₂·e₃)
```

**Example 2**: Mixed basis (canonical 4-gluon case)
```fish
cargo run -p treeamps-cli -- gen-ts --n 4 --deg 3 --ee 1
# Generates 24 structures with exactly 1 EE and 2 PE factors
# Example: (p₁·p₂)(p₁·e₃)(e₂·e₄)
```

**Example 3**: Specify `n` and `ee`, let `deg` be computed
```fish
cargo run -p treeamps-cli -- gen-ts --n 6 --ee 2
# Automatically sets deg = 6 - 2 = 4
# Generates structures with 2 EE and 4 PE factors
```

**Example 4**: Pure EE basis
```fish
cargo run -p treeamps-cli -- gen-ts --n 4 --deg 2 --ee 2
# Generates 3 structures: all pairs of different polarizations
# (e₁·e₂)(e₃·e₄), (e₁·e₃)(e₂·e₄), (e₁·e₄)(e₂·e₃)
```

**Example 5**: Large-scale generation
```fish
# 9-gluon pure EE basis (for benchmarking)
cargo run -p treeamps-cli --release -- gen-ts --n 9 --ee 4

# Warning: n≥10 can be very slow without --release
```

### Output Format

The program prints:
1. Configuration summary (n, deg, ee, eliminated leg)
2. Total count of generated structures
3. Complete list of structures in human-readable form
4. Sanity check results (for known cases)

Example output:
```
Tensor structures (n=4, deg=3, ee=1, elim=p4, one_pol_per_leg=true) count=24
  1) (p1·p2) (p1·e3) (e2·e4)
  2) (p1·p2) (p1·e4) (e2·e3)
  ...
  24) (p2·p3) (p2·e4) (e1·e3)

[Sanity-one-pol-per-leg] expected count=24  (OK)
```

### Performance Tips

1. **Always use `--release` for n ≥ 7**:
   ```fish
   cargo run -p treeamps-cli --release -- gen-ts --n 8 --ee 3
   ```

2. **Benchmark with Criterion**:
   ```fish
   cargo bench
   # Results in target/criterion/report/index.html
   ```

3. **Profile with flamegraph**:
   ```fish
   cargo install flamegraph  # Once
   cargo flamegraph -p treeamps-cli -- gen-ts --n 9 --ee 2
   # Opens flamegraph.svg
   ```

## Development Guide for AI Agents

### Code Organization Principles

1. **Separation of concerns**:
   - `treeamps-core`: Pure Rust library, no I/O
   - `treeamps-cli`: User interface, depends on core

2. **Type safety**:
   - Newtype wrappers (`LegIndex`) prevent index errors
   - Enums for mutually exclusive states
   - No unwrap() in production code paths

3. **Performance**:
   - Prune early in DFS (degree, EE count, polarization constraints)
   - BTreeSet for canonical ordering and deduplication
   - Zero-copy factor catalog references

### Common Modification Patterns

**Adding a new polarization pattern**:
1. Add variant to `PolarizationPattern` enum in [types.rs](treeamps-core/src/types.rs)
2. Update pruning logic in `dfs_emit()` in [generator.rs](treeamps-core/src/generator.rs#L130)
3. Add validation in CLI parameter processing in [main.rs](treeamps-cli/src/main.rs)

**Adding a new particle type**:
1. Create new ScalarKind variants (e.g., `PM` for metric)
2. Update factor catalog generation
3. Modify DFS state tracking for new constraints
4. Add validation formulas for new cases

**Changing eliminated leg convention**:
Currently hardcoded to last leg. To make configurable:
1. Add field to `GenConfig`
2. Update filter in `generate_valid_factors()` [generator.rs](treeamps-core/src/generator.rs#L45)
3. Update CLI to accept `--elim <LEG>` flag

### Testing Strategy

**Current approach**: Built-in sanity checks for known formulas
- 4-gluon mixed basis: 24 structures
- 4-gluon pure-EE basis: 3 structures

**Future needs**:
- Unit tests for individual functions
- Property-based tests (e.g., all structures are canonical)
- Regression tests for nullspace dimensions
- Golden file tests for specific configurations

### Benchmarking Best Practices

Benchmarks live in [treeamps-core/benches/generator.rs](treeamps-core/benches/generator.rs):

```rust
use criterion::{Criterion, criterion_group, criterion_main};
use treeamps_core::generator::{GenConfig, generate_tensor_structures};

fn bench_gen_ts(c: &mut Criterion) {
    let cfg = GenConfig {
        n_legs: 6,
        transversality: Transversality::ForbidPiDotEi,
        pol_pattern: PolarizationPattern::OnePerLeg,
    };
    c.bench_function("gen_ts n6 deg4 ee2 one_per_leg", |b| {
        b.iter(|| generate_tensor_structures(&cfg, 4, 2))
    });
}
```

**To add a benchmark**:
1. Add function to `benches/generator.rs`
2. Add to `criterion_group!` macro
3. Run with `cargo bench`
4. Check results in `target/criterion/`

### Dependencies

**Production**:
- `nalgebra` (0.33): Linear algebra (currently unused, retained for future)
- `num-bigint`, `num-rational`, `num-traits`: Symbolic algebra support (currently unused)
- `clap` (4.x): CLI argument parsing with derive macros

**Development**:
- `criterion` (0.8.1): Benchmarking framework

### Release Checklist

Before marking a version stable:
1. All sanity checks pass
2. Benchmarks show no regression
3. `cargo clippy` clean
4. `cargo fmt` applied
5. Flamegraph shows no obvious bottlenecks
6. Documentation updated

## Physics Background for LLMs

### Why Enumerate Tensor Structures?

In quantum field theory, tree-level scattering amplitudes for gauge bosons (gluons, photons) can be written as:

```
A(1,2,...,n) = Σᵢ cᵢ × Tᵢ(p,e)
```

Where:
- `Tᵢ` are Lorentz-invariant tensor structures (what this code generates)
- `cᵢ` are scalar coefficients (functions of Mandelstam invariants)

### Constraints from Physics

**Gauge invariance**: Replacing any polarization with its momentum must give zero:
```
A(..., eᵢ → pᵢ, ...) = 0
```

This provides linear constraints on the coefficients `cᵢ`.

**Momentum conservation**: 
```
Σᵢ pᵢ = 0
```

This eliminates one momentum (conventionally the last) from independent variables.

**On-shell conditions** (massless particles):
```
pᵢ² = 0  for all i
```

This eliminates factors like `(pᵢ·pᵢ)` and `(eᵢ·eᵢ)`.

**Transversality** (physical polarizations):
```
pᵢ·eᵢ = 0  for all i
```

This eliminates factors like `(pᵢ·eᵢ)`.

### Current Limitations

The generator correctly enumerates **kinematically allowed** tensors but treats all dot products as independent. It does **not** currently:
- Enforce momentum conservation relations
- Apply Mandelstam identities (e.g., `s + t + u = 0` for 4-particle scattering)
- Compute gauge-invariant subspaces
- Determine minimal amplitude bases

These features require symbolic algebra and were removed in the current version for architectural clarity.

## Project History

- **Original C++ implementation**: `tree_amplitudes` directory (reference implementation)
- **Rust rewrite**: This workspace, with improved type safety and performance
- **December 2025**: Streamlined to focus on tensor generation; solver features deferred

## References

For the physics context, see the thesis documents in `/thesis`:
- [sec1.tex](thesis/sec1.tex): Introduction
- [sec2.tex](thesis/sec2.tex): Scattering amplitudes basics  
- [sec3.tex](thesis/sec3.tex): Tensor decomposition methods
- [sec4.tex](thesis/sec4.tex): Applications

## Contributing Guidelines for AI Agents

When modifying this codebase:

1. **Preserve performance**: Verify with `cargo bench` before and after
2. **Maintain conventions**:
   - 1-based leg indexing
   - Last leg always eliminated
   - Canonical factor ordering (sorted)
3. **Update validation**: Add sanity checks for new cases with known counts
4. **Document physics**: Explain the physical meaning of new constraints
5. **Keep separation**: Core library should never depend on CLI

## Contact & License

This is research code developed for a physics thesis. See [main.tex](thesis/main.tex) for academic context.

**Current maintainer**: jonykrem  
**Workspace**: `/Users/jonykrem/Projects/treeamps`

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
