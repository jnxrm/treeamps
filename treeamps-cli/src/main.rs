use clap::{Parser, Subcommand};
use treeamps_core::{GenConfig, generate_tensor_structures};

fn main() {
    let cli = Cli::parse();
    match cli.cmd {
        Command::GenTs { n, deg, ee } => run_gen_ts(n, deg, ee),
        // All solver/symbolic functionality has been removed for now; `solve`
        // is intentionally omitted to keep this CLI focused on tensor-structure
        // generation via `gen-ts`.
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
