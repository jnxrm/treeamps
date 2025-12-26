use std::collections::BTreeSet;

use crate::{
    dot_product::ScalarFactor,
    tensor_structure::TensorStructure,
    types::{LegIndex, PolarizationPattern, ScalarKind, Transversality},
};

/// High-level configuration describing which tensors are allowed.
#[derive(Clone, Debug)]
pub struct GenConfig {
    pub n_legs: u8,
    pub transversality: Transversality,
    pub pol_pattern: PolarizationPattern,
}

impl Default for GenConfig {
    fn default() -> Self {
        Self {
            n_legs: 3,
            transversality: Transversality::ForbidPiDotEi,
            pol_pattern: PolarizationPattern::OnePerLeg,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CatalogCounts {
    pub num_pp: usize,
    pub num_pe: usize,
    pub num_ee: usize,
}

fn generate_valid_factors(
    cfg: &GenConfig,
) -> (Vec<ScalarFactor>, Vec<ScalarFactor>, Vec<ScalarFactor>) {
    let n = cfg.n_legs as u8;

    let mut pp = Vec::new();
    let mut pe = Vec::new();
    let mut ee = Vec::new();

    // PP factors: forbid any factor involving p_n (momentum of leg n)
    for i in 1..=n {
        for j in (i + 1)..=n {
            if i == n || j == n {
                continue;
            }
            pp.push(ScalarFactor::pp(LegIndex(i), LegIndex(j)));
        }
    }

    // PE factors: forbid p_n as momentum, and forbid p_1·e_n
    for i in 1..=n {
        if i == n {
            continue;
        }
        for j in 1..=n {
            if matches!(cfg.transversality, Transversality::ForbidPiDotEi) && i == j {
                continue;
            }
            if j == n && i == 1 {
                continue;
            }
            pe.push(ScalarFactor::pe(LegIndex(i), LegIndex(j)));
        }
    }

    // EE factors
    for i in 1..=n {
        for j in (i + 1)..=n {
            ee.push(ScalarFactor::ee(LegIndex(i), LegIndex(j)));
        }
    }

    pp.sort();
    pe.sort();
    ee.sort();

    (pp, pe, ee)
}

pub fn count_valid_factors(cfg: &GenConfig) -> CatalogCounts {
    let (pp, pe, ee) = generate_valid_factors(cfg);
    CatalogCounts {
        num_pp: pp.len(),
        num_pe: pe.len(),
        num_ee: ee.len(),
    }
}

#[derive(Clone)]
struct DfsState<'a> {
    target_deg: u32,
    ee_needed: u32,
    nlegs: u8,
    enforce_one_pol: bool,
    catalog: &'a [ScalarFactor],
    cur: TensorStructure,
    pe_so_far: u32,
    pol_count: Vec<u32>,
    // out: BTreeSet<TensorStructure>,
}

fn add_polarizations(pc: &mut [u32], f: &ScalarFactor) {
    match f.kind {
        ScalarKind::PE => {
            pc[f.b.0 as usize] += 1;
        }
        ScalarKind::EE => {
            pc[f.a.0 as usize] += 1;
            pc[f.b.0 as usize] += 1;
        }
        ScalarKind::PP => {}
    }
}

fn remove_polarizations(pc: &mut [u32], f: &ScalarFactor) {
    match f.kind {
        ScalarKind::PE => {
            pc[f.b.0 as usize] -= 1;
        }
        ScalarKind::EE => {
            pc[f.a.0 as usize] -= 1;
            pc[f.b.0 as usize] -= 1;
        }
        ScalarKind::PP => {}
    }
}

fn dfs_emit(s: &mut DfsState, idx_start: usize, out: &mut BTreeSet<TensorStructure>) {
    let deg_so_far = s.cur.factors.len() as u32;
    let ee_so_far = s.cur.ee_contractions;

    if deg_so_far > s.target_deg || ee_so_far > s.ee_needed {
        return;
    }

    if s.enforce_one_pol {
        for r in 1..=s.nlegs as usize {
            if s.pol_count[r] > 1 {
                return;
            }
        }

        let remain = s.target_deg - deg_so_far;
        let mut missing = 0;
        for r in 1..=s.nlegs as usize {
            if s.pol_count[r] == 0 {
                missing += 1;
            }
        }

        let max_addable = remain * 2;
        if max_addable < missing {
            return;
        }

        let pol_so_far = 2 * ee_so_far + s.pe_so_far;
        if pol_so_far > s.nlegs as u32 {
            return;
        }
    }

    if deg_so_far == s.target_deg {
        if ee_so_far == s.ee_needed {
            if !s.enforce_one_pol {
                let mut t = s.cur.clone();
                t.canonicalize();
                out.insert(t);
            } else {
                let pol_total = 2 * ee_so_far + s.pe_so_far;
                if pol_total == s.nlegs as u32 {
                    let ok = (1..=s.nlegs as usize).all(|r| s.pol_count[r] == 1);
                    if ok {
                        let mut t = s.cur.clone();
                        t.canonicalize();
                        out.insert(t);
                    }
                }
            }
        }
        return;
    }

    for i in idx_start..s.catalog.len() {
        let f = &s.catalog[i];
        s.cur.factors.push(f.clone());

        if matches!(f.kind, ScalarKind::EE) {
            s.cur.ee_contractions += 1;
        }

        if s.enforce_one_pol {
            if matches!(f.kind, ScalarKind::PE) {
                s.pe_so_far += 1;
            }
            add_polarizations(&mut s.pol_count, f);
        }

        // let out_after = dfs_emit(s.clone(), i);
        // s.out = out_after;

        dfs_emit(s, i, out);

        if s.enforce_one_pol {
            remove_polarizations(&mut s.pol_count, f);
            if matches!(f.kind, ScalarKind::PE) {
                s.pe_so_far -= 1;
            }
        }

        if matches!(f.kind, ScalarKind::EE) {
            s.cur.ee_contractions -= 1;
        }

        s.cur.factors.pop();
    }
}

pub fn generate_tensor_structures(
    cfg: &GenConfig,
    target_degree: u32,
    ee_contractions: u32,
) -> Vec<TensorStructure> {
    if target_degree == 0 {
        return Vec::new();
    }
    if ee_contractions > target_degree {
        return Vec::new();
    }

    let (pp, pe, ee) = generate_valid_factors(cfg);
    let mut catalog = Vec::with_capacity(pp.len() + pe.len() + ee.len());
    catalog.extend(pp);
    catalog.extend(pe);
    catalog.extend(ee);

    let nlegs = cfg.n_legs;
    let mut s = DfsState {
        target_deg: target_degree,
        ee_needed: ee_contractions,
        nlegs,
        enforce_one_pol: matches!(cfg.pol_pattern, PolarizationPattern::OnePerLeg),
        catalog: &catalog,
        cur: TensorStructure::new(),
        pe_so_far: 0,
        pol_count: vec![0; nlegs as usize + 1],
    };

    let mut out_set = BTreeSet::new();
    dfs_emit(&mut s, 0, &mut out_set);
    out_set.into_iter().collect()
}
