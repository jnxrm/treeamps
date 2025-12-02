use crate::types::{LegIndex, ScalarKind};

/// A single scalar factor (dot product) in the tensor basis.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct ScalarFactor {
    pub kind: ScalarKind,
    pub a: LegIndex,
    pub b: LegIndex,
}

impl ScalarFactor {
    pub fn pp(i: LegIndex, j: LegIndex) -> Self {
        Self { kind: ScalarKind::PP, a: i, b: j }
    }
    pub fn pe(i: LegIndex, j: LegIndex) -> Self {
        Self { kind: ScalarKind::PE, a: i, b: j }
    }
    pub fn ee(i: LegIndex, j: LegIndex) -> Self {
        Self { kind: ScalarKind::EE, a: i, b: j }
    }

    pub fn to_string(&self) -> String {
        match self.kind {
            ScalarKind::PP => format!("(p{}·p{})", self.a.0, self.b.0),
            ScalarKind::PE => format!("(p{}·e{})", self.a.0, self.b.0),
            ScalarKind::EE => format!("(e{}·e{})", self.a.0, self.b.0),
        }
    }
}

impl Ord for ScalarFactor {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering::*;
        match self.kind.cmp(&other.kind) {
            Equal => match self.a.cmp(&other.a) {
                Equal => self.b.cmp(&other.b),
                non_eq => non_eq,
            },
            non_eq => non_eq,
        }
    }
}

impl PartialOrd for ScalarFactor {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
