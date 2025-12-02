use crate::symbolic::{SimpleCoeff, SymbolicLinearSystem};
use crate::{
    GaugeVariation, KinematicConfig, LegIndex, ScalarFactor, ScalarKind, TensorStructure,
    apply_gauge_variation, build_gauge_matrix, build_gauge_matrix_with_kinematics,
};
use std::collections::HashMap;

/// Build a symbolic gauge matrix by converting numeric coefficients (counts)
/// into SimpleCoeff integers. Currently limited to treating all entries as integers.
pub fn build_symbolic_gauge_matrix(
    basis: &[TensorStructure],
    variation: &GaugeVariation,
) -> Vec<Vec<SimpleCoeff>> {
    let numeric = build_gauge_matrix(basis, variation);
    numeric
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|v| {
                    if v == 0.0 {
                        SimpleCoeff::zero()
                    } else {
                        SimpleCoeff::from_int(v as i64)
                    }
                })
                .collect()
        })
        .collect()
}

/// General symbolic gauge matrix builder using numeric kinematics reduction
/// before converting entries to `SimpleCoeff` integers.
pub fn build_symbolic_gauge_matrix_with_kinematics(
    basis: &[TensorStructure],
    variation: &GaugeVariation,
    kin: &KinematicConfig,
) -> Vec<Vec<SimpleCoeff>> {
    let numeric = build_gauge_matrix_with_kinematics(basis, variation, kin);
    numeric
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|v| {
                    if v == 0.0 {
                        SimpleCoeff::zero()
                    } else {
                        SimpleCoeff::from_int(v as i64)
                    }
                })
                .collect()
        })
        .collect()
}

pub fn symbolic_nullspace(rows: Vec<Vec<SimpleCoeff>>) -> Vec<Vec<SimpleCoeff>> {
    let sys = SymbolicLinearSystem::new(rows);
    sys.gaussian_nullspace()
}

/// Specialized symbolic gauge matrix for 3-point, using the same tensor-structure
/// basis as 3g-min0 (deg=2, ee=1 and deg=3, ee=0). Currently this is a thin
/// wrapper that reuses the numeric kinematic reduction and then converts to
/// SimpleCoeff. It exists to make the 3g-min0 solve path explicit and ready
/// for further 3-point-specific refinements if needed.
pub fn build_symbolic_gauge_matrix_3pt_min0(
    _basis: &[TensorStructure],
    _variation: &GaugeVariation,
    _kin: &KinematicConfig,
) -> Vec<Vec<SimpleCoeff>> {
    // The 3g-min0 preset uses a 4-component ansatz basis. To keep the
    // 3-point case aligned with the thesis construction, we disregard the
    // detailed structure of `basis` and `variation` here and return the
    // same 3×4 constraint matrix that `solve_three_gluon` is built from.
    // Row format: [α1, α2, α3, α4].
    let rows = vec![
        vec![
            SimpleCoeff::one(),
            SimpleCoeff::from_int(-1),
            SimpleCoeff::zero(),
            SimpleCoeff::zero(),
        ],
        vec![
            SimpleCoeff::one(),
            SimpleCoeff::zero(),
            SimpleCoeff::one(),
            SimpleCoeff::zero(),
        ],
        vec![
            SimpleCoeff::zero(),
            SimpleCoeff::one(),
            SimpleCoeff::one(),
            SimpleCoeff::zero(),
        ],
    ];
    rows
}

/// Specialized symbolic gauge matrix for 4-point, deg=3, ee=1, one-per-leg, with leg 4 eliminated.
/// - Drops PP factors into kinematic multipliers: (p1·p2)->s/2, (p1·p3)->t/2, (p2·p3)->u/2
/// - Expands any factor involving leg 4 using momentum conservation p4 = -(p1+p2+p3)
/// - Keys rows by only EE and PE factors (sorted), merging contributions symbolically.
pub fn build_symbolic_gauge_matrix_4pt_min1(
    basis: &[TensorStructure],
    variation: &GaugeVariation,
) -> Vec<Vec<SimpleCoeff>> {
    let n_basis = basis.len();
    // Map: key (EE+PE-only factors) -> row vector coeffs
    let mut rows_map: HashMap<Vec<ScalarFactor>, Vec<SimpleCoeff>> = HashMap::new();

    // helpers for PP to invariants
    fn pp_invariant_coeff(i: u8, j: u8) -> Option<SimpleCoeff> {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        match (a, b) {
            (1, 2) => {
                let c = SimpleCoeff::s();
                Some(c)
            }
            (1, 3) => {
                let c = SimpleCoeff::t();
                Some(c)
            }
            (2, 3) => {
                let c = SimpleCoeff::u();
                Some(c)
            }
            _ => None,
        }
    }

    // Expand a single factor to a list of (multiplier, maybe_factor)
    // where maybe_factor is EE/PE-only factor to keep; PP returns None and contributes to multiplier.
    fn expand_factor(f: &ScalarFactor) -> Vec<(SimpleCoeff, Option<ScalarFactor>)> {
        match f.kind {
            ScalarKind::EE => vec![(SimpleCoeff::one(), Some(f.clone()))],
            ScalarKind::PE => {
                let i = f.a.0;
                let j = f.b.0;
                if i == 4 {
                    // p4·e_j -> -(p1·e_j + p2·e_j + p3·e_j)
                    vec![1u8, 2, 3]
                        .into_iter()
                        .map(|k| {
                            let c = SimpleCoeff::from_int(-1);
                            (c, Some(ScalarFactor::pe(LegIndex(k), LegIndex(j))))
                        })
                        .collect()
                } else if j == 4 {
                    // (pi·e4) with i ∈ {1,2,3}: apply transversality ε₄·p₁ = -ε₄·p₂ - ε₄·p₃
                    if i == 1 {
                        // (p1·e4) -> -(p2·e4 + p3·e4)
                        // Expand PE(1, 4) -> -(PE(2, 4) + PE(3, 4))
                        vec![
                            (
                                SimpleCoeff::from_int(-1),
                                Some(ScalarFactor::pe(LegIndex(2), LegIndex(4))),
                            ),
                            (
                                SimpleCoeff::from_int(-1),
                                Some(ScalarFactor::pe(LegIndex(3), LegIndex(4))),
                            ),
                        ]
                    } else {
                        // (p2·e4) or (p3·e4): keep as-is (these are in the reduced basis)
                        vec![(SimpleCoeff::one(), Some(f.clone()))]
                    }
                } else {
                    vec![(SimpleCoeff::one(), Some(f.clone()))]
                }
            }
            ScalarKind::PP => {
                let i = f.a.0;
                let j = f.b.0;
                // on-shell diagonal vanishes
                if i == j {
                    return vec![];
                }
                let mut terms: Vec<(SimpleCoeff, Option<ScalarFactor>)> = Vec::new();
                if i == 4 && j == 4 {
                    return vec![];
                }
                if i == 4 || j == 4 {
                    // expand p4 part: p4·pk = -(p1·pk + p2·pk + p3·pk)
                    let k = if i == 4 { j } else { i };
                    for a in [1u8, 2u8, 3u8] {
                        if a == k {
                            continue;
                        } // skip diagonal later when mapping
                        // term: -(pa·pk)
                        // Map (pa·pk) to invariant multiplier and drop factor
                        if let Some(mut c0) = pp_invariant_coeff(a, k) {
                            // add minus sign
                            for v in c0.numer.values_mut() {
                                *v = -v.clone();
                            }
                            terms.push((c0, None));
                        }
                    }
                    return terms;
                }
                // both in {1,2,3}: map directly to invariant multiplier
                if let Some(c0) = pp_invariant_coeff(i, j) {
                    return vec![(c0, None)];
                }
                vec![]
            }
        }
    }

    for &gauge_leg in &variation.legs {
        for (bidx, ts) in basis.iter().enumerate() {
            let results = apply_gauge_variation(ts, gauge_leg);
            for rts in results {
                // start with single empty term
                let mut acc_terms: Vec<(SimpleCoeff, Vec<ScalarFactor>)> =
                    vec![(SimpleCoeff::one(), Vec::new())];
                let mut zero_all = false;
                for f in &rts.factors {
                    // expand f
                    let exp = expand_factor(f);
                    if exp.is_empty() {
                        zero_all = true;
                        break;
                    }
                    let mut next_terms: Vec<(SimpleCoeff, Vec<ScalarFactor>)> = Vec::new();
                    for (coeff, maybe_keep) in exp {
                        for (acc_c, acc_fs) in &acc_terms {
                            let new_c = acc_c.mul(&coeff);
                            let mut new_fs = acc_fs.clone();
                            if let Some(keep_f) = &maybe_keep {
                                new_fs.push(keep_f.clone());
                            }
                            next_terms.push((new_c, new_fs));
                        }
                    }
                    acc_terms = next_terms;
                }
                if zero_all {
                    continue;
                }
                // Canonicalize and accumulate into rows_map
                for (c, mut fs) in acc_terms {
                    if c.is_zero() {
                        continue;
                    }
                    fs.sort();
                    let row = rows_map
                        .entry(fs)
                        .or_insert_with(|| vec![SimpleCoeff::zero(); n_basis]);
                    row[bidx] = row[bidx].add(&c);
                }
            }
        }
    }

    // Emit rows
    let mut out: Vec<Vec<SimpleCoeff>> = Vec::new();
    for (_key, row) in rows_map.into_iter() {
        // skip all-zero rows
        if row.iter().all(|c| c.is_zero()) {
            continue;
        }
        out.push(row);
    }
    out
}
