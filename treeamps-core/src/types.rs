#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct LegIndex(pub u8); // 1-based external leg index

/// Kind of scalar factor: momentum-momentum, momentum-polarization, or polarization-polarization.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum ScalarKind {
    PP,
    PE,
    EE,
}

/// Transversality / p·e rules.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum Transversality {
    None,
    ForbidPiDotEi,
}

/// How polarizations are allowed to appear per leg.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum PolarizationPattern {
    Unrestricted,
    OnePerLeg,
}
