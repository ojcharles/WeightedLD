//! Joseph Roberts & Oscar Charles 2021
use human_format::Formatter;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, log_enabled, Level};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    time::Instant,
};
use structopt::StructOpt;

use weighted_ld::*;

#[derive(Debug, StructOpt)]
#[structopt(
    name = "weighted_ld",
    about = "A tool for computing sequence weighted linkage disequilibrium"
)]
struct Opt {
    #[structopt(long, help = "The source file to load in FASTA or VCF format")]
    fasta_input: PathBuf,

    #[structopt(
        long,
        default_value = "0.8",
        help = "Minimum fraction of A,C,G & T required for a site to be considered in LD and weighting calculations. Increase to account for poor sequence coverage"
    )]
    min_acgt: f32,

    #[structopt(
        long,
        default_value = "0.02",
        help = "The minimum (dominant) minor allele fraction for a site to be considered in LD calculations"
    )]
    min_variability: f32,

    #[structopt(
        long,
        default_value = "0.1",
        help = "Minimum value of R2 to be included in the output"
    )]
    r2_threshold: f32,

    #[structopt(
        long,
        help = "Filename to write the per-sequence weights to, in Tab Separated Value format"
    )]
    weights_output: Option<PathBuf>,

    #[structopt(
        long,
        help = "Filename to write the per-pair weighted LD figures to, in Tab Separated Value format"
    )]
    pair_output: PathBuf,

    #[structopt(
        long,
        help = "Use unit weights instead of Henikoff weights"
    )]
    unweighted: bool,
}

fn write_henikoff_weights(path: &PathBuf, weights: &[f32]) -> Result<(), std::io::Error> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);

    writeln!(w, "Sequence_index\thk_weight")?;
    for (idx, weight) in weights.iter().enumerate() {
        writeln!(w, "{}\t{:.3}", idx, weight)?;
    }

    Ok(())
}

fn write_pair_stats(
    path: &PathBuf,
    pairs: &PairStore<LdStats>,
) -> Result<(), std::io::Error> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);

    let pb = if log_enabled!(Level::Info) {
        let bar = ProgressBar::new(pairs.len() as u64);
        bar.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({per_sec} {eta_precise})")
            .progress_chars("#>-"));
        Some(bar)
    } else {
        None
    };

    let mut written = 0u64;

    writeln!(w, "site_a\tsite_b\tD\tD'\tr2")?;

    for (first_idx, second_idx, ld_stat) in pairs.iter() {
        writeln!(
            w,
            "{}\t{}\t{:.3}\t{:.3}\t{:.3}",
            first_idx, second_idx, ld_stat.d, ld_stat.d_prime, ld_stat.r2
        )?;

        written += 1;
        if written % 5000 == 0 {
            if let Some(pb) = &pb {
                pb.set_position(written);
            }
        }
    }

    Ok(())
}

fn main() -> Result<(), std::io::Error> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let opt = Opt::from_args();

    debug!("{:?}", opt);

    let sw = Instant::now();
    let multiseq = read_fasta(opt.fasta_input)?;
    let siteset = SiteSet::from_multiseq(&multiseq);
    info!("Loaded fasta file in {:?}", sw.elapsed());
    info!(
        "    {} sequences, {} sites",
        siteset.n_seqs(),
        siteset.n_sites(),
    );

    let sw = Instant::now();
    let min_acgt = (opt.min_acgt * siteset.n_seqs() as f32).ceil() as usize;
    let min_minor = opt.min_variability;
    let filtered_siteset =
        siteset.filter_by(|s| is_site_of_interest(s, min_acgt, min_minor));
    info!(
        "Computed + filtered sites of interest in {:?}",
        sw.elapsed()
    );
    info!("    Found {} sites of interest", filtered_siteset.n_sites());

    let weights = if opt.unweighted {
        std::iter::repeat(1f32)
            .take(siteset.n_seqs())
            .collect::<Vec<_>>()
    } else {
        let sw = Instant::now();
        let weights = henikoff_weights(&filtered_siteset);
        info!("Computed Henikoff weights in {:?}", sw.elapsed());
        weights
    };

    if let Some(weights_filepath) = opt.weights_output {
        info!("Writing weights to {:?}", weights_filepath);
        write_henikoff_weights(&weights_filepath, &weights)?;
    }

    info!("Beginning pairwise weighted LD computation");
    let sw = Instant::now();
    let total_pairs = (filtered_siteset.n_sites() - 1) * (filtered_siteset.n_sites() - 2) / 2;
    let pair_store = {
        let pb = if log_enabled!(Level::Info) {
            let bar = ProgressBar::new(total_pairs as u64);
            bar.set_style(ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({per_sec} {eta_precise})")
                .progress_chars("#>-"));
            Some(bar)
        } else {
            None
        };

        all_weighted_ld_pairs(
            &filtered_siteset,
            &weights,
            opt.r2_threshold,
            |computed| {
                if let Some(pb) = &pb {
                    pb.set_position(computed as u64);
                }
            },
        )
    };
    let pair_calc_duration = sw.elapsed();
    info!(
        "Finished computing pairwise weighted LD stats in {:?}",
        pair_calc_duration
    );
    info!(
        "    {} pairs computed at ~{}, {} passed threshold",
        Formatter::new()
            .format(total_pairs as f64),
        Formatter::new()
            .with_units("pairs/s")
            .format(total_pairs as f64 / pair_calc_duration.as_secs_f64()),
        Formatter::new()
            .format(pair_store.len() as f64),
    );

    info!("Writing output to {:?}", opt.pair_output);
    let sw = Instant::now();
    write_pair_stats(&opt.pair_output, &pair_store)?;
    info!("Finshed writing output in {:?}", sw.elapsed());

    Ok(())
}
