use human_format::Formatter;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info};
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
    about = "A tool for computing weighted linkage disequilibrium"
)]
struct Opt {
    #[structopt(long, help = "The source file to load")]
    fasta_input: PathBuf,

    #[structopt(
        long,
        default_value = "0.8",
        help = "Minimum fractions of ACTG for a site to be considered"
    )]
    min_acgt: f32,

    #[structopt(
        long,
        default_value = "0.02",
        help = "Minimum fraction of minor symbols for a site to be considered"
    )]
    min_minor: f32,

    #[structopt(
        long,
        default_value = "0.5",
        help = "Maximum fraction of minor symbols for a site to be considered"
    )]
    max_minor: f32,

    #[structopt(
        long,
        help = "Filename to write the per-sequence Henikoff weights to, in Tab Separated Value format"
    )]
    henikoff_output: Option<PathBuf>,

    #[structopt(
        long,
        help = "Filename to write the per-pair weighted LD figures to, in Tab Separated Value format"
    )]
    pair_output: PathBuf,
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
    pairs: &PairStore,
    source_set: &SiteSet,
) -> Result<(), std::io::Error> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);

    let pb = ProgressBar::new(pairs.len() as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({per_sec} {eta_precise})")
        .progress_chars("#>-"));

    let mut written = 0u64;

    writeln!(w, "site_a\tsite_b\td\td'\tr2")?;

    for first_idx in 0..(source_set.n_sites() - 1) {
        let parent_first_idx = source_set.parent_site_index(first_idx);
        for second_idx in (first_idx + 1)..source_set.n_sites() {
            let parent_second_idx = source_set.parent_site_index(second_idx);
            let ld_stat = match pairs.get_pair(first_idx, second_idx) {
                Some(x) => x,
                None => continue,
            };

            writeln!(
                w,
                "{}\t{}\t{:.3}\t{:.3}\t{:.3}",
                parent_first_idx, parent_second_idx, ld_stat.d, ld_stat.d_prime, ld_stat.r2
            )?;

            written += 1;
            if written % 5000 == 0 {
                pb.set_position(written);
            }
        }
    }

    Ok(())
}

fn main() -> Result<(), std::io::Error> {
    env_logger::init();
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
    let min_acgt = (opt.min_acgt * siteset.n_seqs() as f32).ceil() as u32;
    let min_minor = opt.min_minor;
    let max_minor = opt.max_minor;
    let filtered_siteset = siteset.filter_by(|s| is_site_of_interest(s, min_acgt, min_minor, max_minor));
    info!(
        "Computed + filtered sites of interest in {:?}",
        sw.elapsed()
    );
    info!("    Found {} sites of interest", filtered_siteset.n_sites());

    let sw = Instant::now();
    let weights_hk = henikoff_weights(&filtered_siteset);
    info!("Computed Henikoff weights in {:?}", sw.elapsed());
    if let Some(hk_filepath) = opt.henikoff_output {
        info!("Writing Henikoff weights to {:?}", hk_filepath);
        write_henikoff_weights(&hk_filepath, &weights_hk)?;
    }

    let mut pair_store = PairStore::new(filtered_siteset.n_sites());

    info!("Beginning pairwise weighted LD computation");
    let sw = Instant::now();
    {
        let pb = ProgressBar::new(pair_store.len() as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent}% ({per_sec} {eta_precise})")
            .progress_chars("#>-"));

        all_weighted_ld_pairs(
            &filtered_siteset,
            &weights_hk,
            &mut pair_store,
            |computed| {
                pb.set_position(computed as u64);
            },
        );
    }
    let pair_calc_duration = sw.elapsed();
    info!(
        "Finished computing pairwise weighted LD stats in {:?}",
        pair_calc_duration
    );
    info!(
        "    ~{}",
        Formatter::new()
            .with_units("pairs/s")
            .format(pair_store.len() as f64 / pair_calc_duration.as_secs_f64())
    );

    info!("Writing output to {:?}", opt.pair_output);
    let sw = Instant::now();
    write_pair_stats(&opt.pair_output, &pair_store, &filtered_siteset)?;
    info!("Finshed writing output in {:?}", sw.elapsed());

    Ok(())
}
