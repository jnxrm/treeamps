// Tensor-structure (TS) subsystem: combinatorics and canonical representation
pub mod types;
pub mod dot_product;
pub mod tensor_structure;
pub mod generator;

// Constraint / solver subsystem: gauge variation, kinematics, and linear algebra
pub mod constraints;
pub mod kinematics;
pub mod linear_system;
pub mod symbolic;
pub mod symbolic_constraints;

// Public TS API
pub use crate::types::{LegIndex, ScalarKind, PolarizationPattern, Transversality};
pub use crate::dot_product::ScalarFactor;
pub use crate::tensor_structure::TensorStructure;
pub use crate::generator::{GenConfig, CatalogCounts, generate_tensor_structures};

// Public constraint / solve API
pub use crate::constraints::{
	apply_gauge_variation, GaugeVariation, build_gauge_matrix,
	build_gauge_matrix_with_kinematics,
};
pub use crate::kinematics::KinematicConfig;
pub use crate::linear_system::LinearSystem;
pub use crate::symbolic_constraints::{
	build_symbolic_gauge_matrix,
	build_symbolic_gauge_matrix_with_kinematics,
	symbolic_nullspace,
};
pub use crate::symbolic::{
	four_gluon_min1_solution, solve_four_gluon_min1, solve_three_gluon,
	solve_two_gluon_two_scalar,
};
