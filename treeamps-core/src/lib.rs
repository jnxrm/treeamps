// Tensor-structure (TS) subsystem: combinatorics and canonical representation
pub mod dot_product;
pub mod generator;
pub mod tensor_structure;
pub mod types;

// Public TS API only
pub use crate::dot_product::ScalarFactor;
pub use crate::generator::{CatalogCounts, GenConfig, generate_tensor_structures};
pub use crate::tensor_structure::TensorStructure;
pub use crate::types::{LegIndex, PolarizationPattern, ScalarKind, Transversality};
