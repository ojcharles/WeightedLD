use std::{
    fs::File,
    io::{BufRead, BufReader},
    iter::Sum,
    ops::{Index, IndexMut},
    path::{Path, PathBuf},
    sync::{
        Mutex,
        atomic::{AtomicUsize, Ordering},
    }
};

use human_format::{Formatter, Scales};
use log::{info, log_enabled, Level};
use ndarray::prelude::*;
use num_derive::FromPrimitive;
use rayon::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq, FromPrimitive)]
#[repr(u8)]
pub enum Symbol {
    A = 0,
    C,
    G,
    T,
    Missing,
    Unknown,
}

impl Symbol {
    pub fn is_acgt(&self) -> bool {
        match self {
            Symbol::A | Symbol::C | Symbol::G | Symbol::T => true,
            Symbol::Missing | Symbol::Unknown => false,
        }
    }
}

impl Default for Symbol {
    fn default() -> Self {
        Self::Unknown
    }
}

impl From<char> for Symbol {
    fn from(c: char) -> Self {
        match c {
            'a' | 'A' => Self::A,
            'c' | 'C' => Self::C,
            'g' | 'G' => Self::G,
            't' | 'T' => Self::T,
            '-' => Self::Missing,
            _ => Self::Unknown,
        }
    }
}

pub struct SymbolHistogram<T: num::Integer> {
    data: [T; 6],
}

impl<T: num::Integer> Index<Symbol> for SymbolHistogram<T> {
    type Output = T;

    fn index(&self, index: Symbol) -> &Self::Output {
        &self.data[index as usize]
    }
}

impl<T: num::Integer> IndexMut<Symbol> for SymbolHistogram<T> {
    fn index_mut(&mut self, index: Symbol) -> &mut Self::Output {
        &mut self.data[index as usize]
    }
}

impl<T: num::Integer + Copy + Sum> SymbolHistogram<T> {
    fn zero() -> Self {
        Self {
            data: [T::zero(); 6],
        }
    }

    fn from_slice(symbols: &[Symbol]) -> Self {
        let mut hist = Self::zero();
        for symbol in symbols {
            hist[*symbol] = hist[*symbol].add(T::one());
        }
        hist
    }

    fn from_slice_mask(symbols: &[Symbol], mask: &[bool]) -> Self {
        debug_assert_eq!(symbols.len(), mask.len());

        let mut hist = Self::zero();
        for (symbol, mask) in symbols.iter().zip(mask.iter()) {
            if !mask {
                continue;
            }
            hist[*symbol] = hist[*symbol].add(T::one());
        }
        hist
    }

    fn acgt(&self) -> T {
        use Symbol::*;
        self[A] + self[C] + self[G] + self[T]
    }

    fn major_symbol(&self) -> Option<Symbol> {
        let mut major = None;

        for sym in &[Symbol::A, Symbol::C, Symbol::G, Symbol::T] {
            if self[*sym] > major.map(|m| self[m]).unwrap_or(T::zero()) {
                major = Some(*sym);
            }
        }

        major
    }
}

pub struct Sequence {
    pub name: Option<String>,
    pub symbols: Vec<Symbol>,
}

pub enum MultiSequenceSource {
    FastaFile(PathBuf),
}

pub struct MultiSequence {
    pub source: MultiSequenceSource,
    pub sequences: Vec<Sequence>,
}

/// Stores symbol data such that all data for a given site is contiguous
pub struct SiteSet {
    n_sites: usize,
    n_seqs: usize,

    buffer: Vec<Symbol>,

    /// Vector of length n_sites, mapping a site index in this set to some other index
    ///
    /// Used when this SiteSet is only a subset of the total number of sites.
    /// The site with index N in this set has an index `site_map[N]` in the parent set
    site_map: Option<Vec<usize>>,
}

impl Index<usize> for SiteSet {
    type Output = [Symbol];

    fn index(&self, index: usize) -> &Self::Output {
        let start = index * self.n_seqs;
        &self.buffer[start..(start + self.n_seqs)]
    }
}

impl SiteSet {
    pub fn from_multiseq(ms: &MultiSequence) -> Self {
        let n_seqs = ms.sequences.len();
        let n_sites = ms.sequences[0].symbols.len();

        if ms.sequences.iter().any(|s| s.symbols.len() != n_sites) {
            panic!("Not all sequences have the same number of symbols");
        }

        let buffer_len = n_seqs * n_sites;
        let mut buffer = vec![Symbol::Missing; buffer_len];

        for site in 0..n_sites {
            for seq in 0..n_seqs {
                buffer[site * n_seqs + seq] = ms.sequences[seq].symbols[site];
            }
        }

        Self {
            n_sites,
            n_seqs,
            buffer,
            site_map: None,
        }
    }

    pub fn filter_by(&self, filter_func: impl Fn(&[Symbol]) -> bool) -> Self {
        let mut new_buffer = Vec::with_capacity(self.buffer.len());
        let mut site_map = Vec::new();

        for site in 0..self.n_sites {
            let site_slice = &self[site];
            if filter_func(site_slice) {
                site_map.push(site);
                new_buffer.extend_from_slice(site_slice);
            }
        }

        Self {
            n_sites: site_map.len(),
            n_seqs: self.n_seqs,
            buffer: new_buffer,
            site_map: Some(site_map),
        }
    }

    /// Get a reference to the site set's n sites.
    pub fn n_sites(&self) -> usize {
        self.n_sites
    }

    /// Get a reference to the site set's n seqs.
    pub fn n_seqs(&self) -> usize {
        self.n_seqs
    }

    pub fn parent_site_index(&self, idx: usize) -> usize {
        self.site_map.as_ref().map(|m| m[idx]).unwrap_or(idx)
    }
}

pub fn read_fasta<P: AsRef<Path>>(path: P) -> Result<MultiSequence, std::io::Error> {
    let source = MultiSequenceSource::FastaFile(path.as_ref().to_path_buf());
    let mut sequences = Vec::new();

    let file = File::open(path.as_ref())?;
    let mut reader = BufReader::new(file);

    let mut line = String::new();
    let mut name = None;

    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break;
        }

        if let Some(new_name) = line.strip_prefix(">") {
            name = Some(new_name.to_string());
        } else {
            let symbols = line.chars().map(Symbol::from).collect::<Vec<_>>();

            sequences.push(Sequence {
                name: name.take(),
                symbols,
            })
        }
    }

    Ok(MultiSequence { source, sequences })
}

/// Given a slice of all symbols in a site, should the site be considered for further computations
pub fn is_site_of_interest(site: &[Symbol], min_acgt: u32, min_minor: f32, max_minor: f32) -> bool {
    let hist = SymbolHistogram::<u32>::from_slice(site);

    let acgt_count = hist.acgt();

    if acgt_count <= min_acgt {
        // There aren't enough ACGT symbols at this site
        false
    } else {
        let major_sym = match hist.major_symbol() {
            Some(m) => m,
            _ => return false,
        };
        
        let minor_frac = (acgt_count as f32 - hist[major_sym] as f32) / hist.acgt() as f32;

        if minor_frac < min_minor || minor_frac > max_minor{
            // The fraction of minor symbols is either too low or too high
            false
        } else {
            true
        }
    }
}

pub fn henikoff_weights(data: &SiteSet) -> Vec<f32> {
    let mut contributions = Array2::<f32>::zeros((data.n_sites(), data.n_seqs()));

    for site in 0..data.n_sites() {
        let site_slice = &data[site];

        let mut contrib_row = contributions.row_mut(site);
        let contrib_slice = contrib_row
            .as_slice_mut()
            .expect("Expected contributions to be dense in memory");

        henikoff_site_contributions(site_slice, contrib_slice);
    }

    let mut weights = contributions.sum_axis(Axis(0));
    weights /= weights.fold(0.0f32, |a, b| a.max(*b));
    assert!(weights.is_standard_layout());
    weights.into_raw_vec()
}

pub fn henikoff_site_contributions(site: &[Symbol], contributions: &mut [f32]) {
    let hist = SymbolHistogram::<u32>::from_slice(site);

    let acgt_count = hist.acgt() as f32;
    let mut total_site_contrib = 0f32;

    for (idx, sym) in site.iter().enumerate() {
        if sym.is_acgt() {
            contributions[idx] = 1f32 / (acgt_count * hist[*sym] as f32);
            total_site_contrib += contributions[idx];
        }
    }

    let mean_site_contrib = total_site_contrib / acgt_count;

    for (idx, sym) in site.iter().enumerate() {
        if !sym.is_acgt() {
            contributions[idx] = mean_site_contrib;
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct LdStats {
    pub r2: f32,
    pub d: f32,
    pub d_prime: f32,
}

#[allow(non_snake_case)]
pub fn single_weighted_ld_pair(a: &[Symbol], b: &[Symbol], weights: &[f32]) -> Option<LdStats> {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), weights.len());

    let good_mask = a
        .iter()
        .zip(b.iter())
        .map(|(a, b)| a.is_acgt() & b.is_acgt())
        .collect::<Vec<bool>>();

    let a_hist = SymbolHistogram::<u32>::from_slice_mask(a, &good_mask);
    let b_hist = SymbolHistogram::<u32>::from_slice_mask(b, &good_mask);

    let a_maj_sym = match a_hist.major_symbol() {
        Some(maj) => maj,
        _ => return None,
    };

    let b_maj_sym = match b_hist.major_symbol() {
        Some(maj) => maj,
        _ => return None,
    };

    let mut PA = 0f32;
    let mut Pa = 0f32;
    let mut PB = 0f32;
    let mut Pb = 0f32;
    let mut ld_obs = [0f32; 4];

    for seq in 0..a.len() {
        if !good_mask[seq] {
            continue;
        }

        if a[seq] == a_maj_sym {
            PA += weights[seq];
        } else {
            Pa += weights[seq];
        }

        if b[seq] == b_maj_sym {
            PB += weights[seq];
        } else {
            Pb += weights[seq];
        }

        match (a[seq] == a_maj_sym, b[seq] == b_maj_sym) {
            (false, false) => ld_obs[0] += weights[seq],
            (true, true) => ld_obs[3] += weights[seq],
            (false, true) => ld_obs[1] += weights[seq],
            (true, false) => ld_obs[2] += weights[seq],
        }
    }

    let total_weight: f32 = weights
        .iter()
        .zip(good_mask.iter())
        .filter(|(_w, m)| **m)
        .map(|(w, _m)| w)
        .sum();

    PA /= total_weight;
    Pa /= total_weight;
    PB /= total_weight;
    Pb /= total_weight;

    let PAB = PA * PB;
    let PAb = PA * Pb;
    let PaB = Pa * PB;
    let Pab = Pa * Pb;

    ld_obs[0] /= total_weight;
    ld_obs[1] /= total_weight;
    ld_obs[2] /= total_weight;
    ld_obs[3] /= total_weight;

    let d = ((PAB - ld_obs[3]) + (Pab - ld_obs[0]) + (ld_obs[2] - PAb) + (ld_obs[1] - PaB)) / 4f32;

    let mut denominator: f32;
    if d < 0f32 {
        denominator = (-ld_obs[0]).max(-ld_obs[3]);
        if denominator == 0f32 {
            denominator = (-ld_obs[0]).min(-ld_obs[3]);
        }
    } else {
        denominator = (-ld_obs[1]).min(-ld_obs[2]);
        if denominator == 0f32 {
            denominator = (-ld_obs[1]).max(-ld_obs[2]);
        }
    }
    let d_prime = d / denominator;

    let r2 = d * d / (PA * Pa * PB * Pb);

    Some(LdStats { r2, d, d_prime })
}

pub struct PairStore {
    n_sites: usize,
    buffer: Vec<Option<LdStats>>,
}

impl PairStore {
    pub fn new(n_sites: usize) -> Self {
        if log_enabled!(Level::Info) {
            let expected_mem_usage = n_sites * n_sites * std::mem::size_of::<Option<LdStats>>();
            let expected_mem_usage = Formatter::new()
                .with_scales(Scales::Binary())
                .with_units("B")
                .format(expected_mem_usage as f64);
            info!(
                "Creating a new in-memory store for pairwise LD data, expected memory usage: {}",
                expected_mem_usage
            )
        }

        Self {
            n_sites,
            buffer: vec![None; n_sites * n_sites],
        }
    }
    
    pub fn len(&self) -> usize {
        (self.n_sites - 1) * (self.n_sites - 2) / 2
    }

    fn linear_idx(&self, first_idx: usize, second_idx: usize) -> usize {
        debug_assert_ne!(first_idx, second_idx);
        let (first_idx, second_idx) = if first_idx < second_idx {
            (first_idx, second_idx)
        } else {
            (second_idx, first_idx)
        };

        // TODO: Inefficient - uses ~twice as much memory as needed
        first_idx + second_idx * self.n_sites
    }

    pub fn get_pair(&self, first_idx: usize, second_idx: usize) -> Option<LdStats> {
        self.buffer[self.linear_idx(first_idx, second_idx)]
    }

    pub fn write_pair(&mut self, first_idx: usize, second_idx: usize, ld_stat: Option<LdStats>) {
        let linear_idx = self.linear_idx(first_idx, second_idx);
        self.buffer[linear_idx] = ld_stat;
    }

    /// As write_pair, but with an immutable receiver.
    ///
    /// It is up to the caller to ensure that the given pair of indices is not being read from or
    /// written to concurrently by any other thread of execution.
    unsafe fn write_pair_unsafe(&self, first_idx: usize, second_idx: usize, ld_stat: Option<LdStats>) {
        let linear_idx = self.linear_idx(first_idx, second_idx);
        let ptr = &self.buffer[linear_idx] as *const _ as *mut Option<LdStats>;
        ptr.write(ld_stat);
    }
}

pub fn all_weighted_ld_pairs(
    site_set: &SiteSet,
    weights: &[f32],
    store: &mut PairStore,
    mut progress_report: impl FnMut(usize) + Send,
) {
    progress_report(0);

    let computed_pair_count = AtomicUsize::new(0);
    let progress_report = Mutex::new(progress_report);

    for first_idx in 0..(site_set.n_sites() - 1) {
        ((first_idx + 1)..site_set.n_sites()).into_par_iter()
            .for_each(|second_idx| {
                let a = &site_set[first_idx];
                let b = &site_set[second_idx];
                let ld_stat = single_weighted_ld_pair(a, b, weights);

                unsafe { store.write_pair_unsafe(first_idx, second_idx, ld_stat); }

                let computed = computed_pair_count.fetch_add(1, Ordering::Relaxed);
                // Avoid spamming the progress bar every single time
                if computed % 5_000 == 0 {
                    (progress_report.lock().expect("Failed to get lock on progress indicator callback"))(computed);
                }
            });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_histogram_from_slice() {
        use Symbol::*;
        let a = [A, A, A, C, C, G, T, T, T, T];
        let hist = SymbolHistogram::<u8>::from_slice(&a);

        assert_eq!(hist[A], 3);
        assert_eq!(hist[C], 2);
        assert_eq!(hist[G], 1);
        assert_eq!(hist[T], 4);
    }

    #[test]
    fn test_hist_major_minor() {
        use Symbol::*;

        let hist = SymbolHistogram {
            data: [0, 1, 10, 2, 0, 0],
        };
        let maj = hist.major_symbol();
        assert_eq!(maj, Some(G));

        let hist = SymbolHistogram {
            data: [10, 1, 9, 2, 0, 0],
        };
        let maj = hist.major_symbol();
        assert_eq!(maj, Some(A));
    }
}
