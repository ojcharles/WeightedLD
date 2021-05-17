use std::{
    fs::File,
    io::{BufRead, BufReader},
    iter::Sum,
    ops::{Index, IndexMut},
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicUsize, Ordering},
        Mutex,
    },
};

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

    pub fn is_acgtm(&self) -> bool {
        match self {
            Symbol::A | Symbol::C | Symbol::G | Symbol::T | Symbol::Missing => true,
            Symbol::Unknown => false,
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

pub struct SymbolHistogram {
    data: [usize; 6],
}

impl Index<Symbol> for SymbolHistogram {
    type Output = usize;

    fn index(&self, index: Symbol) -> &Self::Output {
        &self.data[index as usize]
    }
}

impl IndexMut<Symbol> for SymbolHistogram {
    fn index_mut(&mut self, index: Symbol) -> &mut Self::Output {
        &mut self.data[index as usize]
    }
}

impl SymbolHistogram {
    fn zero() -> Self {
        Self {
            data: [0; 6],
        }
    }

    fn from_slice(symbols: &[Symbol]) -> Self {
        let mut hist = Self::zero();
        for symbol in symbols {
            hist[*symbol] += 1;
        }
        hist
    }

    fn acgt(&self) -> usize {
        use Symbol::*;
        self[A] + self[C] + self[G] + self[T]
    }

    fn acgtm(&self) -> usize {
        use Symbol::*;
        self[A] + self[C] + self[G] + self[T] + self[Missing]
    }
    
    fn distinct_known_count(&self) -> usize {
        use Symbol::*;

        (if self[A] > 0 { 1usize } else { 0 }) +
        if self[C] > 0 { 1 } else { 0 } +
        if self[G] > 0 { 1 } else { 0 } +
        if self[T] > 0 { 1 } else { 0 } +
        if self[Missing] > 0 { 1 } else { 0 }
    }

    fn major_minor_symbols(&self) -> (Option<Symbol>, Option<Symbol>) {
        let mut major = None;
        let mut minor = None;

        for sym in &[Symbol::A, Symbol::C, Symbol::G, Symbol::T, Symbol::Missing] {
            if self[*sym] > major.map(|m| self[m]).unwrap_or(0) {
                minor = major;
                major = Some(*sym);
            } else if self[*sym] > minor.map(|m| self[m]).unwrap_or(0) {
                minor = Some(*sym);
            }
        }

        (major, minor)
    }
}

pub struct Sequence {
    pub name: Option<String>,
    pub symbols: Vec<Symbol>,
}

pub enum MultiSequenceSource {
    FastaFile(PathBuf),
    None,
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

    #[cfg(test)]
    pub fn from_strs(raw: &[&str]) -> Self {
        let sequences = raw
            .iter()
            .map(|seq_str| {
                seq_str
                    .chars()
                    .map(Symbol::from)
                    .collect::<Vec<_>>()
            })
            .map(|symbols| Sequence {
                name: None,
                symbols,
            })
            .collect::<Vec<Sequence>>();

        Self::from_multiseq(&MultiSequence {
            source: MultiSequenceSource::None,
            sequences,
        })
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
pub fn is_site_of_interest(site: &[Symbol], min_acgt: usize, min_minor: f32, max_minor: f32) -> bool {
    let hist = SymbolHistogram::from_slice(site);

    let acgt_count = hist.acgt();

    if acgt_count <= min_acgt {
        // There aren't enough ACGT symbols at this site
        false
    } else {
        let (major_sym, minor_sym) = match hist.major_minor_symbols() {
            (Some(maj), Some(min)) => (maj, min),
            _ => return false,
        };

        let maj_count = hist[major_sym] as f32;
        let min_count = hist[minor_sym] as f32;

        // let minor_frac = min_count / acgt_count as f32;
        let minor_frac = min_count / (min_count + maj_count);
        // let minor_frac = (acgt_count as f32 - min_count) / acgt_count as f32;

        if minor_frac < min_minor || minor_frac > max_minor {
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
    let hist = SymbolHistogram::from_slice(site);

    let distinct_count = hist.distinct_known_count() as f32;
    let mut total_site_contrib = 0f32;

    for (idx, sym) in site.iter().enumerate() {
        if sym.is_acgtm() {
            contributions[idx] = 1f32 / (distinct_count * hist[*sym] as f32);
            total_site_contrib += contributions[idx];
        }
    }

    let mean_site_contrib = total_site_contrib / distinct_count;

    for (idx, sym) in site.iter().enumerate() {
        if !sym.is_acgtm() {
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

    // TODO: It is wasteful to compute these histograms for every pair - should
    // lift their computations outside the pairwise loop.
    let a_hist = SymbolHistogram::from_slice(a);
    let b_hist = SymbolHistogram::from_slice(b);

    let (a_maj_sym, a_min_sym) = match a_hist.major_minor_symbols() {
        (Some(maj), Some(min)) => (maj, min),
        _ => return None,
    };

    let (b_maj_sym, b_min_sym) = match b_hist.major_minor_symbols() {
        (Some(maj), Some(min)) => (maj, min),
        _ => return None,
    };

    let good_mask = a
        .iter()
        .zip(b.iter())
        .map(|(&a, &b)| (a == a_maj_sym || a == a_min_sym) && (b == b_maj_sym || b == b_min_sym))
        .collect::<Vec<bool>>();

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
            // good_mask has already filtered down the sequences to ones where
            // a is either a_maj_sym or a_min_sym -> if a isn't a_maj_sym it
            // must be a_min_sym.
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
        denominator = (ld_obs[1]).min(ld_obs[2]);
        if denominator == 0f32 {
            denominator = (ld_obs[1]).max(ld_obs[2]);
        }
    }
    let d_prime = d / denominator;

    let r2 = d * d / (PA * Pa * PB * Pb);

    Some(LdStats { r2, d, d_prime })
}

struct PairData<T> {
    first_idx: usize,
    second_idx: usize,
    data: T,
}

pub struct PairStore<T> {
    chunks: Vec<Vec<PairData<T>>>,
}

impl<T> PairStore<T> {
    pub fn iter(&self) -> PairStoreIter<T> {
        PairStoreIter {
            store: self,
            chunk_idx: 0,
            pair_idx: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.chunks.iter().map(|c| c.len()).sum()
    }
}

pub struct PairStoreIter<'a, T> {
    store: &'a PairStore<T>,

    /// Index of the next chunk
    chunk_idx: usize,

    /// Index of the next pair within the next chunk
    pair_idx: usize,
}

impl<'a, T> Iterator for PairStoreIter<'a, T> {
    type Item = (usize, usize, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        let ret = if self.chunk_idx < self.store.chunks.len() {
            let elem = &self.store.chunks[self.chunk_idx][self.pair_idx];
            (elem.first_idx, elem.second_idx, &elem.data)
        } else {
            return None;
        };

        self.pair_idx += 1;
        if self.pair_idx >= self.store.chunks[self.chunk_idx].len() {
            self.pair_idx = 0;
            self.chunk_idx += 1;
        }

        Some(ret)
    }
}

pub fn all_weighted_ld_pairs(
    site_set: &SiteSet,
    weights: &[f32],
    r2_threshold: f32,
    mut progress_report: impl FnMut(usize) + Send,
) -> PairStore<LdStats> {
    progress_report(0);

    let computed_pair_count = AtomicUsize::new(0);
    let progress_report = Mutex::new(progress_report);

    let data_chunks = (0..(site_set.n_sites() - 1))
        .into_par_iter()
        .map(|first_idx| {
            let a = &site_set[first_idx];
            let mut results_chunk = Vec::new();
            for second_idx in (first_idx + 1)..site_set.n_sites() {
                let b = &site_set[second_idx];
                if let Some(ld_stat) = single_weighted_ld_pair(a, b, weights) {
                    if ld_stat.r2 > r2_threshold {
                        results_chunk.push(PairData {
                            first_idx: site_set.parent_site_index(first_idx),
                            second_idx: site_set.parent_site_index(second_idx),
                            data: ld_stat,
                        });
                    }
                }
            }

            let computed = computed_pair_count
                .fetch_add(site_set.n_sites() - first_idx - 1, Ordering::Relaxed);
            (progress_report
                .lock()
                .expect("Failed to get lock on progress indicator callback"))(computed);

            results_chunk
        })
        .filter(|chunk| chunk.len() > 0)
        .collect::<Vec<Vec<PairData<LdStats>>>>();

    PairStore {
        chunks: data_chunks,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::*;

    #[test]
    fn test_histogram_from_slice() {
        use Symbol::*;
        let a = [A, A, A, C, C, G, T, T, T, T, Missing, Missing, Unknown];
        let hist = SymbolHistogram::from_slice(&a);

        assert_eq!(hist[A], 3);
        assert_eq!(hist[C], 2);
        assert_eq!(hist[G], 1);
        assert_eq!(hist[T], 4);
        assert_eq!(hist[Missing], 2);
        assert_eq!(hist[Unknown], 1);
    }
    #[test]
    fn test_hist_major_minor() {
        use Symbol::*;

        let hist = SymbolHistogram {
            data: [0, 1, 10, 2, 0, 0],
        };
        let (maj, min) = hist.major_minor_symbols();
        assert_eq!(maj, Some(G));
        assert_eq!(min, Some(T));

        let hist = SymbolHistogram {
            data: [1, 9, 10, 2, 0, 0],
        };
        let (maj, min) = hist.major_minor_symbols();
        assert_eq!(maj, Some(G));
        assert_eq!(min, Some(C));

        let hist = SymbolHistogram {
            data: [1, 1, 40, 2, 4, 0],
        };
        let (maj, min) = hist.major_minor_symbols();
        assert_eq!(maj, Some(G));
        assert_eq!(min, Some(Missing));
    }

    #[test]
    fn test_henikoff_weights_1() {
        // as in the Henikoff and henikoff 1994 paper
        let siteset = SiteSet::from_strs(&["AAAAA", "AAAAA", "CCCCC", "CCCCC", "TTTTT"]);
        assert_ulps_eq!(henikoff_weights(&siteset)[..], [0.5, 0.5, 0.5, 0.5, 1.0]);
    }

    #[test]
    fn test_henikoff_weights_2() {
        // as in S.F. Altschul NIH
        let siteset = SiteSet::from_strs(&["GCGTTAGC", "GAGTTGGA", "CGGACTAA"]);
        assert_abs_diff_eq!(henikoff_weights(&siteset)[..], [0.769, 0.692, 1.0], epsilon = 0.001); // todo ever so slightly off could be float point error
    }

    #[test]
    fn test_henikoff_weights_3() {
        // ensure that indels are treated equally - seq 2 is most unique
        // 0.9166 , 1.25, 0.9166, 0.9166 -> 0.7333, 1, 0.7333, 0.7333
        let siteset = SiteSet::from_strs(&["AAGA", "AA-A", "GGGG", "GGGG"]);
        assert_abs_diff_eq!(henikoff_weights(&siteset)[..], [0.733, 1.0, 0.733, 0.733], epsilon = 0.001); // todo handle assert nearly equal
    }

    #[test]
    fn test_ld_pair_unweighted_ld0() {
        use Symbol::*;
        let a = [A, A, A, A, T, T, T, T ];
        let b = [T, T, A, A, A, A, T, T ];
        let weights = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0];

        let ld_stats = single_weighted_ld_pair(&a, &b, &weights)
            .expect("Expected test case to have LD statistics available");

        assert_abs_diff_eq!(ld_stats.d, 0.0, epsilon = 1e-5);
        assert_abs_diff_eq!(ld_stats.d_prime, 0.0, epsilon = 1e-5);
        assert_abs_diff_eq!(ld_stats.r2, 0.0, epsilon = 1e-5);
    }

    #[test]
    fn test_ld_pair_unweighted_ld1() {
        use Symbol::*;
        let a = [A, A, A, A, T, T, T, T ];
        let b = [T, T, T, T, A, A, A, A ];
        let weights = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0];

        let ld_stats = single_weighted_ld_pair(&a, &b, &weights)
            .expect("Expected test case to have LD statistics available");
        println!("{}",&ld_stats.d);
        assert_abs_diff_eq!(ld_stats.d, 0.25, epsilon = 1e-5);
        assert_abs_diff_eq!(ld_stats.d_prime, 0.5, epsilon = 1e-5);
        assert_abs_diff_eq!(ld_stats.r2, 1.0, epsilon = 1e-5);
    }

    #[test]
    fn test_single_weighted_ld_pair() { //todo paper and pen this example as a sanity check
        use Symbol::*;
        let a = [A, A, A, A, C, A, C];
        let b = [A, A, A, G, T, A, A];
        let weights = [1.0, 1.0, 0.4, 0.2, 0.5, 0.8, 0.2];

        let ld_stats = single_weighted_ld_pair(&a, &b, &weights)
            .expect("Expected test case to have LD statistics available");

        assert_abs_diff_eq!(ld_stats.d, 0.00308, epsilon = 1e-5);
        assert_abs_diff_eq!(ld_stats.d_prime, 0.05555, epsilon = 1e-5);
        assert_abs_diff_eq!(ld_stats.r2, 0.00346, epsilon = 1e-5);
    }
}
