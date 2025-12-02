use crate::dot_product::ScalarFactor;

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct TensorStructure {
    pub factors: Vec<ScalarFactor>,
    pub ee_contractions: u32,
}

impl TensorStructure {
    pub fn new() -> Self {
        Self { factors: Vec::new(), ee_contractions: 0 }
    }

    pub fn canonicalize(&mut self) {
        self.factors.sort();
    }

    pub fn to_string(&self) -> String {
        if self.factors.is_empty() {
            return "1".to_string();
        }
        self.factors
            .iter()
            .map(|f| f.to_string())
            .collect::<Vec<_>>()
            .join(" · ")
    }
}

impl Ord for TensorStructure {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.factors.cmp(&other.factors)
    }
}

impl PartialOrd for TensorStructure {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
