use clap::{ArgAction, Parser, Subcommand};
use treeamps_core::{
    GenConfig, PolarizationPattern, four_gluon_min1_solution, generate_tensor_structures,
    solve_four_gluon_min1, solve_three_gluon, solve_two_gluon_two_scalar,
};

fn main() {
    let cli = Cli::parse();
    match cli.cmd {
        Command::Demo => run_demo(),
        Command::FourGluonSymbolic => run_four_gluon_symbolic(),
        Command::ThreeGluonSymbolic => run_three_gluon_symbolic(),
        Command::TwoGluonTwoScalarSymbolic => run_two_gluon_two_scalar_symbolic(),
        Command::GenTs { n, deg, ee } => run_gen_ts(n, deg, ee),
        Command::Solve {
            n,
            deg,
            min_deg,
            max_deg,
            min_ee,
            gauge,
            kinematics,
            tol,
            symbolic,
            show_unreduced,
            one_pol_per_leg,
            preset,
        } => {
            // Apply preset overrides
            let (n, deg, min_deg, max_deg, min_ee, gauge_legs, use_kinematics) =
                if let Some(preset_name) = preset {
                    match preset_name.as_str() {
                        // 3-gluon min-0: use degrees 2 and 3 with ee=1,0 respectively.
                        // We encode this via degree range [2,3] and min_ee=0; the
                        // basis builder in run_solve will special-case n=3.
                        "3g-min0" => (3, None, Some(2), Some(3), 0, vec![1, 2, 3], kinematics),
                        "3g-min1" => (3, Some(2), None, Some(2), 1, vec![1, 2, 3], kinematics),
                        "4g-min1" => (4, Some(3), None, Some(3), 1, vec![1, 2, 3, 4], kinematics),
                        _ => (n, deg, min_deg, max_deg, min_ee, gauge, kinematics),
                    }
                } else {
                    (n, deg, min_deg, max_deg, min_ee, gauge, kinematics)
                };

            run_solve(
                n,
                deg,
                min_deg,
                max_deg,
                min_ee,
                gauge_legs,
                use_kinematics,
                tol,
                symbolic,
                show_unreduced,
                one_pol_per_leg,
            );
        }
    }
}

#[derive(Parser, Debug)]
#[command(
    name = "treeamps",
    about = "Tree amplitude tensor-structure explorer (Rust rewrite)"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Small demo basis for n=3
    Demo,

    /// Symbolic 4-gluon (min 1 metric contraction) example from thesis (deprecated; use Solve --preset)
    FourGluonSymbolic,

    /// Symbolic 3-gluon example (deprecated; use Solve --preset)
    ThreeGluonSymbolic,

    /// Symbolic 2-gluon 2-scalar example
    TwoGluonTwoScalarSymbolic,

    /// Generate tensor structures for fixed degree and EE count
    GenTs {
        /// Number of external legs
        #[arg(long, default_value_t = 3)]
        n: u8,

        /// Total number of factors (degree); leave 0 to infer from n and ee
        #[arg(long, default_value_t = 0)]
        deg: u32,

        /// Number of EE contractions; leave 0 to infer from n and deg
        #[arg(long, default_value_t = 0)]
        ee: u32,
    },

    /// Solve for gauge-invariant structures via nullspace
    Solve {
        /// Number of external legs
        #[arg(long, default_value_t = 3)]
        n: u8,

        /// Total number of factors (degree)
        #[arg(long)]
        deg: Option<u32>,

        /// Minimum degree (inclusive)
        #[arg(long)]
        min_deg: Option<u32>,

        /// Maximum degree (inclusive)
        #[arg(long)]
        max_deg: Option<u32>,

        /// Minimum number of EE contractions
        #[arg(long, default_value_t = 0)]
        min_ee: u32,

        /// Gauge legs (multiple allowed)
        #[arg(long, num_args = 1.., value_delimiter = ',')]
        gauge: Vec<u8>,

        /// Use kinematic reduction (momentum conservation, Mandelstam)
        #[arg(long, action = ArgAction::SetTrue)]
        kinematics: bool,

        /// Nullspace tolerance (numeric mode)
        #[arg(long, default_value_t = 1e-10)]
        tol: f64,

        /// Use symbolic mode (integer/rational coefficients)
        #[arg(long, action = ArgAction::SetTrue)]
        symbolic: bool,

        /// Show unreduced 4-gluon solution (thesis cross-check)
        #[arg(long, action = ArgAction::SetTrue)]
        show_unreduced: bool,

        /// (Currently unused) placeholder; one-pol-per-leg is treated as standard
        #[arg(long, hide = true)]
        one_pol_per_leg: bool,

        /// Use preset for canonical cases (3g-min0, 3g-min1, 4g-min1)
        #[arg(long)]
        preset: Option<String>,
    },
}

fn run_demo() {
    let mut cfg = GenConfig::default();
    cfg.n_legs = 3;
    let ts = generate_tensor_structures(&cfg, 2, 1);
    println!("Basis (size={})", ts.len());
    for (i, t) in ts.iter().enumerate() {
        println!("  {}) {}", i + 1, t.to_string());
    }
}

fn run_four_gluon_symbolic() {
    println!(
        "4-gluon symbolic example (min 1 metric contraction)\nOrdering matches thesis (27 structures).\nCoefficients shown relative to overall normalisation α_27 = 1."
    );
    let solved = solve_four_gluon_min1();
    let expected_dim = 1usize; // thesis: unique solution up to normalisation
    if solved.is_empty() {
        println!("No symbolic solution found (unexpected). Using hard-coded thesis vector.");
        let sol = four_gluon_min1_solution();
        for (i, expr) in sol.iter().enumerate() {
            println!("α_{:>2} = {}", i + 1, expr);
        }
        println!(
            "\n[Regression] Expected nullspace dimension={}  => observed={}  (MISMATCH)",
            expected_dim, 0
        );
    } else {
        println!("Computed solution:");
        let sol = &solved[0];
        for (i, coeff) in sol.iter().enumerate() {
            println!("α_{:>2} = {}", i + 1, coeff.to_string());
        }
        println!(
            "\n[Regression] Expected nullspace dimension={}  => observed={}{}",
            expected_dim,
            solved.len(),
            if solved.len() == expected_dim {
                "  (OK)"
            } else {
                "  (MISMATCH)"
            }
        );
        println!("\nThesis (cross-check):");
        let thesis = four_gluon_min1_solution();
        for (i, expr) in thesis.iter().enumerate() {
            println!("α_{:>2} = {}", i + 1, expr);
        }
    }
    println!(
        "\nVariables: s, t, u with s + t + u = 0 (4-point kinematics).\nYou can substitute u = -s - t if desired."
    );
}

fn run_three_gluon_symbolic() {
    println!("3-gluon symbolic system (nullspace basis):");
    let sols = solve_three_gluon();
    let expected_dim = 2usize; // thesis: two solutions for min 0 contractions
    for (k, sol) in sols.iter().enumerate() {
        println!("Solution {}:", k + 1);
        for (i, coeff) in sol.iter().enumerate() {
            println!("  α{} = {}", i + 1, coeff.to_string());
        }
    }
    println!(
        "\n[Regression] Expected nullspace dimension={}  => observed={}{}",
        expected_dim,
        sols.len(),
        if sols.len() == expected_dim {
            "  (OK)"
        } else {
            "  (MISMATCH)"
        }
    );
}

fn run_two_gluon_two_scalar_symbolic() {
    println!("2 gluons + 2 scalars symbolic system (nullspace basis):");
    let sols = solve_two_gluon_two_scalar();
    for (k, sol) in sols.iter().enumerate() {
        println!("Solution {}:", k + 1);
        for (i, coeff) in sol.iter().enumerate() {
            println!("  α{} = {}", i + 1, coeff.to_string());
        }
    }
    println!("Note: 4-point kinematics relation u = -s - t may be applied.");
}

fn run_gen_ts(n: u8, mut deg: u32, mut ee: u32) {
    if n == 0 {
        eprintln!("--n must be >= 1");
        std::process::exit(1);
    }
    // For gluon bases we always enforce "one polarization per leg".
    // The constraint 2*EE + PE = n and deg = EE + PE implies
    // deg = n - ee and ee = n - deg. Enforce consistency and
    // allow one to be inferred from the other when left as zero.
    {
        let implied_deg = n as u32 - ee;
        let implied_ee = n as u32 - deg;

        // If both deg and ee are nonzero, require mutual consistency.
        if deg != 0 && ee != 0 {
            if deg != implied_deg || ee != implied_ee {
                eprintln!(
                    "Inconsistent inputs for one-pol-per-leg: n = {}, deg = {}, ee = {}. Expected deg = n - ee = {} and ee = n - deg = {}.",
                    n, deg, ee, implied_deg, implied_ee,
                );
                std::process::exit(1);
            }
        } else if deg == 0 && ee != 0 {
            deg = implied_deg;
        } else if ee == 0 && deg != 0 {
            ee = implied_ee;
        } else if deg == 0 && ee == 0 {
            // Default: pure PE basis with no EE contractions (min-0 case)
            deg = n as u32;
            ee = 0;
        }
    }

    if ee > deg {
        eprintln!("--ee must be <= --deg");
        std::process::exit(1);
    }

    let mut cfg = GenConfig::default();
    cfg.n_legs = n;
    cfg.pol_pattern = PolarizationPattern::OnePerLeg;

    let ts = generate_tensor_structures(&cfg, deg, ee);
    println!(
        "Tensor structures (n={}, deg={}, ee={}, elim=p{}, one_pol_per_leg={}) count={}",
        n,
        deg,
        ee,
        n,
        true,
        ts.len()
    );
    for (i, t) in ts.iter().enumerate() {
        println!("  {}) {}", i + 1, t.to_string());
    }

    // Canonical sanity checks for the 4-leg case, mirroring the C++ tool
    if n == 4 {
        // Mixed (EE)(PE)(PE) basis with one polarization per leg
        if deg == 3 && ee == 1 {
            let expected_one_pol = 24i64;
            println!(
                "\n[Sanity-one-pol-per-leg] expected count={}{}",
                expected_one_pol,
                if expected_one_pol == ts.len() as i64 {
                    "  (OK)"
                } else {
                    "  (MISMATCH)"
                }
            );
        }

        // Pure EE basis with one polarization per leg: three structures
        if deg == 2 && ee == 2 {
            let expected_pure_ee = 3i64;
            println!(
                "[Sanity-4g-pure-EE-one-pol] expected count={}{}",
                expected_pure_ee,
                if expected_pure_ee == ts.len() as i64 {
                    "  (OK)"
                } else {
                    "  (MISMATCH)"
                }
            );
        }
    }
}

fn run_solve(
    n: u8,
    deg: Option<u32>,
    min_deg: Option<u32>,
    max_deg: Option<u32>,
    min_ee: u32,
    _gauge_legs: Vec<u8>,
    _use_kinematics: bool,
    _tol: f64,
    _symbolic_mode: bool,
    _show_unreduced: bool,
    one_pol_per_leg: bool,
    // Optional: enforce thesis-known zero coefficients for 4-pt min1
    // Future: expose via CLI flag if desired
) {
    // Temporary stripped-down solver: only explore tensor-structure bases.
    // Gauge legs, kinematics, and nullspace computation are intentionally
    // ignored here so we can redesign the solver from scratch on top of a
    // trusted `gen-ts` basis.

    // Determine degree range just for reporting.
    let (min_deg, max_deg) = if let Some(d) = deg {
        (d, d)
    } else {
        let min = min_deg.unwrap_or((n as u32 + 1) / 2);
        let max = max_deg.unwrap_or(n as u32);
        (min, max)
    };

    let mut cfg = GenConfig::default();
    cfg.n_legs = n;
    cfg.pol_pattern = if one_pol_per_leg {
        PolarizationPattern::OnePerLeg
    } else {
        PolarizationPattern::Unrestricted
    };

    // Build union basis over degree range using the same generator as `gen-ts`.
    println!(
        "[Explore] Building basis for n={}, elim=p{}, deg ∈ [{}, {}]",
        n, n, min_deg, max_deg
    );

    let mut basis = Vec::new();
    if n == 3 && min_deg == 2 && max_deg == 3 && min_ee == 0 {
        // 3g-min0 preset: union basis using the same generic generator that
        // `gen-ts` uses. We do not hard-code counts here; instead, we
        // mirror whatever `treeamps-core` considers valid for n=3.
        let mut s_ee1 = generate_tensor_structures(&cfg, 2, 1);
        println!("  deg=2, ee=1: {} structures", s_ee1.len());
        basis.append(&mut s_ee1);

        let mut s_ee0 = generate_tensor_structures(&cfg, 3, 0);
        println!("  deg=3, ee=0: {} structures", s_ee0.len());
        basis.append(&mut s_ee0);
    } else {
        for deg_val in min_deg..=max_deg {
            // For each degree, use ee = n - deg (thesis convention: deg + ee = n)
            let ee = if n as u32 >= deg_val {
                n as u32 - deg_val
            } else {
                0
            };

            // Apply min_ee filter
            if ee < min_ee {
                continue;
            }

            let mut structures = generate_tensor_structures(&cfg, deg_val, ee);
            println!(
                "  deg={}, ee={}: {} structures",
                deg_val,
                ee,
                structures.len()
            );
            basis.append(&mut structures);
        }
    }

    println!("\n[Explore] Total basis size: {}", basis.len());

    println!("[Explore] Basis elements:");
    for (i, t) in basis.iter().enumerate() {
        println!("  {:>3}) {}", i + 1, t.to_string());
    }

    println!(
        "\n[Note] Solver functionality (gauge matrices, nullspaces) is currently disabled.\n       This subcommand now serves only as a tensor-structure explorer over\n       degree ranges, ready for a fresh solver design."
    );
}
