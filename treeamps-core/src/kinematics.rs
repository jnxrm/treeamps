use crate::{LegIndex, ScalarFactor, ScalarKind};
use std::collections::HashMap;

/// Represents a symbolic linear combination of scalar factors with rational coefficients.
/// Used to express kinematic relations like p_4 = -(p_1 + p_2 + p_3).
#[derive(Debug, Clone)]
pub struct SymbolicSum {
    /// Maps ScalarFactor to its coefficient (rational as f64)
    pub terms: HashMap<ScalarFactor, f64>,
}

impl PartialEq for SymbolicSum {
    fn eq(&self, other: &Self) -> bool {
        if self.terms.len() != other.terms.len() {
            return false;
        }
        for (k, v) in &self.terms {
            if let Some(&other_v) = other.terms.get(k) {
                if (v - other_v).abs() > 1e-12 {
                    return false;
                }
            } else {
                return false;
            }
        }
        true
    }
}

impl Eq for SymbolicSum {}

impl std::hash::Hash for SymbolicSum {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // Hash based on sorted keys for deterministic ordering
        let mut items: Vec<_> = self.terms.iter().collect();
        items.sort_by_key(|(k, _)| *k);
        for (factor, coeff) in items {
            factor.hash(state);
            // Hash coefficient as bits to avoid floating point issues
            coeff.to_bits().hash(state);
        }
    }
}

impl SymbolicSum {
    pub fn new() -> Self {
        Self {
            terms: HashMap::new(),
        }
    }

    pub fn single(factor: ScalarFactor, coeff: f64) -> Self {
        let mut sum = Self::new();
        sum.add_term(factor, coeff);
        sum
    }

    pub fn add_term(&mut self, factor: ScalarFactor, coeff: f64) {
        if coeff.abs() < 1e-14 {
            return;
        }
        *self.terms.entry(factor).or_insert(0.0) += coeff;
        // Clean up near-zero entries
        self.terms.retain(|_, &mut v| v.abs() >= 1e-14);
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Substitute momentum for eliminated leg using momentum conservation.
    /// For n legs with leg k eliminated: p_k = -sum_{i != k} p_i
    pub fn apply_momentum_conservation(&mut self, n_legs: u8, eliminate_leg: LegIndex) {
        let mut new_terms = HashMap::new();

        for (factor, coeff) in &self.terms {
            match factor.kind {
                ScalarKind::PP => {
                    // Handle (p_i · p_j) with momentum conservation
                    let (i, j) = (factor.a, factor.b);
                    
                    if i == eliminate_leg && j == eliminate_leg {
                        // (p_k · p_k) = 0 (on-shell, massless)
                        continue;
                    } else if i == eliminate_leg {
                        // (p_k · p_j) = -(sum_{m != k} p_m) · p_j = -sum_{m != k} (p_m · p_j)
                        for m in 1..=n_legs {
                            let leg_m = LegIndex(m);
                            if leg_m != eliminate_leg {
                                let new_factor = ScalarFactor::pp(leg_m.min(j), leg_m.max(j));
                                *new_terms.entry(new_factor).or_insert(0.0) -= coeff;
                            }
                        }
                    } else if j == eliminate_leg {
                        // (p_i · p_k) = same as above by symmetry
                        for m in 1..=n_legs {
                            let leg_m = LegIndex(m);
                            if leg_m != eliminate_leg {
                                let new_factor = ScalarFactor::pp(i.min(leg_m), i.max(leg_m));
                                *new_terms.entry(new_factor).or_insert(0.0) -= coeff;
                            }
                        }
                    } else {
                        // (p_i · p_j) unchanged
                        *new_terms.entry(factor.clone()).or_insert(0.0) += coeff;
                    }
                }
                ScalarKind::PE => {
                    // (p_i · ε_j)
                    if factor.a == eliminate_leg {
                        // (p_k · ε_j) = -sum_{m != k} (p_m · ε_j)
                        for m in 1..=n_legs {
                            let leg_m = LegIndex(m);
                            if leg_m != eliminate_leg {
                                let new_factor = ScalarFactor::pe(leg_m, factor.b);
                                *new_terms.entry(new_factor).or_insert(0.0) -= coeff;
                            }
                        }
                    } else {
                        *new_terms.entry(factor.clone()).or_insert(0.0) += coeff;
                    }
                }
                ScalarKind::EE => {
                    // (ε_i · ε_j) unchanged by momentum conservation
                    *new_terms.entry(factor.clone()).or_insert(0.0) += coeff;
                }
            }
        }

        self.terms = new_terms;
        self.terms.retain(|_, &mut v| v.abs() >= 1e-14);
    }
}

/// Configuration for kinematic simplification.
#[derive(Debug, Clone)]
pub struct KinematicConfig {
    pub n_legs: u8,
    pub eliminate_leg: LegIndex,
}

impl KinematicConfig {
    pub fn new(n_legs: u8, eliminate_leg: LegIndex) -> Self {
        Self {
            n_legs,
            eliminate_leg,
        }
    }

    /// Apply momentum conservation to a scalar factor.
    /// Returns a SymbolicSum representing the factor after applying p_k = -sum_{i!=k} p_i.
    pub fn simplify_factor(&self, factor: &ScalarFactor) -> SymbolicSum {
        let mut sum = SymbolicSum::single(factor.clone(), 1.0);
        sum.apply_momentum_conservation(self.n_legs, self.eliminate_leg);
        sum
    }

    /// Get all independent PP scalar factors (Mandelstam variables).
    /// For n particles with one eliminated, there are C(n-1, 2) independent s_ij.
    pub fn independent_pp_factors(&self) -> Vec<ScalarFactor> {
        let mut factors = Vec::new();
        for i in 1..=self.n_legs {
            if LegIndex(i) == self.eliminate_leg {
                continue;
            }
            for j in (i + 1)..=self.n_legs {
                if LegIndex(j) == self.eliminate_leg {
                    continue;
                }
                factors.push(ScalarFactor::pp(LegIndex(i), LegIndex(j)));
            }
        }
        factors
    }

    /// Express a PP factor in terms of independent Mandelstam variables using momentum conservation.
    /// For 4 legs with leg 4 eliminated, we have s, t, u where s+t+u=0, so only 2 are independent.
    pub fn express_pp_factor(&self, factor: &ScalarFactor) -> SymbolicSum {
        assert_eq!(factor.kind, ScalarKind::PP);
        self.simplify_factor(factor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_momentum_conservation_pp() {
        let config = KinematicConfig::new(4, LegIndex(4));
        
        // (p_1 · p_4) should become -(p_1·p_1) - (p_1·p_2) - (p_1·p_3)
        // = 0 - (p_1·p_2) - (p_1·p_3)  [since p_1·p_1 = 0 on-shell]
        let factor = ScalarFactor::pp(LegIndex(1), LegIndex(4));
        let result = config.simplify_factor(&factor);
        
        assert_eq!(result.terms.len(), 2);
        assert_eq!(result.terms.get(&ScalarFactor::pp(LegIndex(1), LegIndex(2))), Some(&-1.0));
        assert_eq!(result.terms.get(&ScalarFactor::pp(LegIndex(1), LegIndex(3))), Some(&-1.0));
    }

    #[test]
    fn test_momentum_conservation_pe() {
        let config = KinematicConfig::new(4, LegIndex(4));
        
        // (p_4 · ε_1) should become -(p_1·ε_1) - (p_2·ε_1) - (p_3·ε_1)
        let factor = ScalarFactor::pe(LegIndex(4), LegIndex(1));
        let result = config.simplify_factor(&factor);
        
        assert_eq!(result.terms.len(), 3);
        assert_eq!(result.terms.get(&ScalarFactor::pe(LegIndex(1), LegIndex(1))), Some(&-1.0));
        assert_eq!(result.terms.get(&ScalarFactor::pe(LegIndex(2), LegIndex(1))), Some(&-1.0));
        assert_eq!(result.terms.get(&ScalarFactor::pe(LegIndex(3), LegIndex(1))), Some(&-1.0));
    }
}
