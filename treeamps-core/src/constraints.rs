use crate::{KinematicConfig, LegIndex, ScalarFactor, ScalarKind, TensorStructure};
use std::collections::HashMap;

/// Represents a gauge variation: ε_r → p_r for specified legs.
#[derive(Debug, Clone)]
pub struct GaugeVariation {
    pub legs: Vec<LegIndex>,
}

impl GaugeVariation {
    pub fn new(legs: Vec<LegIndex>) -> Self {
        Self { legs }
    }
}

/// Represents a linear combination of tensor structures with coefficients.
/// Used to express the result of gauge substitutions.
#[derive(Debug, Clone)]
pub struct LinearCombination {
    /// Maps tensor structure index to coefficient
    pub terms: HashMap<usize, f64>,
}

impl LinearCombination {
    pub fn new() -> Self {
        Self {
            terms: HashMap::new(),
        }
    }

    pub fn add_term(&mut self, basis_idx: usize, coeff: f64) {
        *self.terms.entry(basis_idx).or_insert(0.0) += coeff;
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

/// Apply gauge variation ε_r → p_r to a single scalar factor.
/// Returns a vector of (coefficient, new_factor) pairs.
fn apply_gauge_to_factor(factor: &ScalarFactor, gauge_leg: LegIndex) -> Vec<(f64, ScalarFactor)> {
    match factor.kind {
        ScalarKind::PP => {
            // (p_i · p_j) is unchanged
            vec![(1.0, factor.clone())]
        }
        ScalarKind::PE => {
            // (p_i · ε_j): if j == gauge_leg, becomes (p_i · p_j)
            if factor.b == gauge_leg {
                vec![(1.0, ScalarFactor::pp(factor.a, gauge_leg))]
            } else {
                vec![(1.0, factor.clone())]
            }
        }
        ScalarKind::EE => {
            // (ε_i · ε_j):
            // - if i == gauge_leg: (ε_i · ε_j) → (p_i · ε_j)
            // - if j == gauge_leg: (ε_i · ε_j) → (ε_i · p_j)
            // - if both: (p_i · p_j)
            let i_is_gauge = factor.a == gauge_leg;
            let j_is_gauge = factor.b == gauge_leg;

            if i_is_gauge && j_is_gauge {
                // (ε_r · ε_r) → (p_r · p_r) = 0 (on-shell)
                vec![]
            } else if i_is_gauge {
                vec![(1.0, ScalarFactor::pe(gauge_leg, factor.b))]
            } else if j_is_gauge {
                vec![(1.0, ScalarFactor::pe(gauge_leg, factor.a))]
            } else {
                vec![(1.0, factor.clone())]
            }
        }
    }
}

/// Apply gauge variation ε_r → p_r to a tensor structure.
/// Returns a list of resulting tensor structures (before looking up in basis).
pub fn apply_gauge_variation(ts: &TensorStructure, gauge_leg: LegIndex) -> Vec<TensorStructure> {
    // Quick check: if no factors involve gauge_leg, return clone unchanged
    let mut has_gauge_leg = false;
    for factor in &ts.factors {
        match factor.kind {
            ScalarKind::PP => {}
            ScalarKind::PE => {
                if factor.b == gauge_leg {
                    has_gauge_leg = true;
                    break;
                }
            }
            ScalarKind::EE => {
                if factor.a == gauge_leg || factor.b == gauge_leg {
                    has_gauge_leg = true;
                    break;
                }
            }
        }
    }

    if !has_gauge_leg {
        return vec![ts.clone()];
    }

    // For each factor in the tensor structure, apply the gauge variation
    // This creates a product of sums, which we expand into a sum of products
    let mut results: Vec<Vec<ScalarFactor>> = vec![vec![]];

    for factor in &ts.factors {
        let substitutions = apply_gauge_to_factor(factor, gauge_leg);

        if substitutions.is_empty() {
            // Factor vanishes (e.g., (ε_r · ε_r) → 0), entire product vanishes
            return vec![];
        }

        // Optimization: if substitution is identity, just append to all results
        if substitutions.len() == 1 && substitutions[0].1 == *factor {
            for r in &mut results {
                r.push(factor.clone());
            }
        } else {
            let mut new_results = Vec::new();
            for current_factors in results {
                for (coeff, new_factor) in &substitutions {
                    if *coeff == 0.0 {
                        continue;
                    }
                    let mut extended = current_factors.clone();
                    extended.push(new_factor.clone());
                    new_results.push(extended);
                }
            }
            results = new_results;
        }
    }

    // Convert factor lists to TensorStructures
    results
        .into_iter()
        .map(|factors| {
            let ee_count = factors.iter().filter(|f| f.kind == ScalarKind::EE).count() as u32;
            let mut ts = TensorStructure {
                factors,
                ee_contractions: ee_count,
            };
            ts.canonicalize();
            ts
        })
        .collect()
}

/// Build a gauge constraint matrix.
/// Each row represents: for a specific gauge leg, one resulting structure that appears
/// when gauge-varying the basis.
/// The constraint is: sum over basis of (coeff_ij * alpha_j) = 0 for each row i.
pub fn build_gauge_matrix(basis: &[TensorStructure], variation: &GaugeVariation) -> Vec<Vec<f64>> {
    let n_basis = basis.len();

    // For each gauge leg separately, collect all resulting structures
    // and build one constraint row per distinct post-variation structure,
    // in line with the bootstrap construction in the thesis.
    let mut all_rows = Vec::new();

    for &gauge_leg in &variation.legs {
        // Map: resulting structure -> list of (basis_idx, coefficient)
        let mut result_map: HashMap<TensorStructure, Vec<(usize, f64)>> = HashMap::new();

        for (basis_idx, ts) in basis.iter().enumerate() {
            let results = apply_gauge_variation(ts, gauge_leg);
            for result_ts in results {
                result_map
                    .entry(result_ts)
                    .or_insert_with(Vec::new)
                    .push((basis_idx, 1.0));
            }
        }

        // Each unique resulting structure gives one constraint row
        for (_result_ts, contributions) in result_map {
            let mut row = vec![0.0; n_basis];
            for (basis_idx, coeff) in contributions {
                row[basis_idx] += coeff;
            }
            // Only add non-trivial rows
            if row.iter().any(|&x| x.abs() > 1e-14) {
                all_rows.push(row);
            }
        }
    }

    all_rows
}

/// Build a gauge constraint matrix with kinematic simplification.
/// This version applies momentum conservation and Mandelstam substitutions.
pub fn build_gauge_matrix_with_kinematics(
    basis: &[TensorStructure],
    variation: &GaugeVariation,
    kinematic_config: &KinematicConfig,
) -> Vec<Vec<f64>> {
    let n_basis = basis.len();

    // For each gauge leg, build constraints with kinematic simplification
    let mut all_rows = Vec::new();

    for &gauge_leg in &variation.legs {
        // Map: canonical simplified key -> list of (basis_idx, coefficient)
        let mut result_map: HashMap<Vec<ScalarFactor>, Vec<(usize, f64)>> = HashMap::new();

        for (basis_idx, ts) in basis.iter().enumerate() {
            let results = apply_gauge_variation(ts, gauge_leg);

            for result_ts in results {
                // Simplify each factor in the result using kinematics
                let mut simplified_factors = Vec::new();
                let mut zero_structure = false;

                for factor in &result_ts.factors {
                    let symbolic = kinematic_config.simplify_factor(factor);

                    if symbolic.is_zero() {
                        zero_structure = true;
                        break;
                    }

                    // For now, use the dominant term as representative
                    // (Full symbolic expansion would be more complex)
                    if let Some((sf, _coeff)) = symbolic.terms.iter().next() {
                        simplified_factors.push(sf.clone());
                    }
                }

                if zero_structure {
                    continue;
                }

                // Canonicalize: sort factors
                simplified_factors.sort();

                result_map
                    .entry(simplified_factors)
                    .or_insert_with(Vec::new)
                    .push((basis_idx, 1.0));
            }
        }

        // Each unique simplified result gives one constraint row
        for (_key, contributions) in result_map {
            let mut row = vec![0.0; n_basis];
            for (basis_idx, coeff) in contributions {
                row[basis_idx] += coeff;
            }
            // Only add non-zero rows
            if row.iter().any(|&x| x.abs() > 1e-14) {
                all_rows.push(row);
            }
        }
    }

    all_rows
}
