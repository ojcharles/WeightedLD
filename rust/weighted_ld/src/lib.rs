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
use packed_simd::{f32x8, u8x8};
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

unsafe impl bytemuck::Zeroable for Symbol {}

unsafe impl bytemuck::Pod for Symbol {}

#[derive(Clone)]
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

    pub fn from_slice(symbols: &[Symbol]) -> Self {
        let mut hist = Self::zero();
        for symbol in symbols {
            hist[*symbol] = hist[*symbol].add(T::one());
        }
        hist
    }

    pub fn acgt(&self) -> T {
        use Symbol::*;
        self[A] + self[C] + self[G] + self[T]
    }

    pub fn major_minor_symbols(&self) -> (Option<Symbol>, Option<Symbol>) {
        let mut major = None;
        let mut minor = None;

        for sym in &[Symbol::A, Symbol::C, Symbol::G, Symbol::T] {
            if self[*sym] > major.map(|m| self[m]).unwrap_or(T::zero()) {
                minor = major;
                major = Some(*sym);
            } else if self[*sym] > minor.map(|m| self[m]).unwrap_or(T::zero()) {
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

    /// Vector of length n_sites, containing the symbol histogram for each site
    histograms: Vec<SymbolHistogram<u32>>,
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
        let mut histograms = Vec::new();

        for site in 0..n_sites {
            let buff_start_idx = site * n_seqs;
            for seq in 0..n_seqs {
                buffer[buff_start_idx + seq] = ms.sequences[seq].symbols[site];
            }

            let hist =
                SymbolHistogram::from_slice(&buffer[buff_start_idx..(buff_start_idx + n_seqs)]);
            histograms.push(hist);
        }

        Self {
            n_sites,
            n_seqs,
            buffer,
            histograms,
            site_map: None,
        }
    }

    pub fn filter_by(&self, filter_func: impl Fn(&[Symbol]) -> bool) -> Self {
        let mut new_buffer = Vec::with_capacity(self.buffer.len());
        let mut site_map = Vec::new();
        let mut new_histograms = Vec::new();

        for site in 0..self.n_sites {
            let site_slice = &self.site_symbols(site);
            if filter_func(site_slice) {
                site_map.push(site);
                new_buffer.extend_from_slice(site_slice);
                new_histograms.push(self.site_histogram(site).clone());
            }
        }

        Self {
            n_sites: site_map.len(),
            n_seqs: self.n_seqs,
            buffer: new_buffer,
            site_map: Some(site_map),
            histograms: new_histograms,
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

    pub fn site_symbols(&self, index: usize) -> &[Symbol] {
        let start = index * self.n_seqs;
        &self.buffer[start..(start + self.n_seqs)]
    }

    pub fn site_histogram(&self, index: usize) -> &SymbolHistogram<u32> {
        &self.histograms[index]
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
        let site_slice = data.site_symbols(site);

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
pub fn single_weighted_ld_pair(
    a: &[Symbol],
    a_hist: &SymbolHistogram<u32>,
    b: &[Symbol],
    b_hist: &SymbolHistogram<u32>,
    weights: &[f32],
) -> Option<LdStats> {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), weights.len());

    let (a_maj_sym, a_min_sym) = match a_hist.major_minor_symbols() {
        (Some(maj), Some(min)) => (maj, min),
        _ => return None,
    };

    let (b_maj_sym, b_min_sym) = match b_hist.major_minor_symbols() {
        (Some(maj), Some(min)) => (maj, min),
        _ => return None,
    };

    let mut total_weight = f32x8::splat(0.0);
    let mut PA = f32x8::splat(0.0);
    let mut PB = f32x8::splat(0.0);
    let mut ld_3 = f32x8::splat(0.0);

    let simd_end = (a.len() / 8) * 8;
    for seq in (0..simd_end).step_by(8) {
        let a = bytemuck::cast_slice::<Symbol, u8>(&a[seq..(seq + 8)]);
        let a = u8x8::from_slice_unaligned(a);

        let b = bytemuck::cast_slice::<Symbol, u8>(&b[seq..(seq + 8)]);
        let b = u8x8::from_slice_unaligned(b);

        let a_maj = a.eq(u8x8::splat(a_maj_sym as u8));
        let a_min = a.eq(u8x8::splat(a_min_sym as u8));
        let b_maj = b.eq(u8x8::splat(b_maj_sym as u8));
        let b_min = b.eq(u8x8::splat(b_min_sym as u8));

        let mask = (a_maj | a_min) & (b_maj | b_min);
        let weight = f32x8::from_slice_unaligned(&weights[seq..(seq + 8)]);

        let zero = f32x8::splat(0.0);
        let weight = mask.select(weight, zero);

        total_weight += weight;
        PA += a_maj.select(weight, zero);
        PB += b_maj.select(weight, zero);
        ld_3 += (a_maj & b_maj).select(weight, zero);
    }

    let mut total_weight = total_weight.sum();
    let mut PA = PA.sum();
    let mut PB = PB.sum();
    let mut ld_obs = [0.0, 0.0, 0.0, ld_3.sum()];

    for seq in simd_end..a.len() {
        if !(a[seq] == a_maj_sym || a[seq] == a_min_sym) {
            continue;
        }
        if !(b[seq] == b_maj_sym || b[seq] == b_min_sym) {
            continue;
        }

        total_weight += weights[seq];
        match (a[seq] == a_maj_sym, b[seq] == b_maj_sym) {
            (true, true) => {
                PA += weights[seq];
                ld_obs[3] += weights[seq];
            }
            (true, false) => PA += weights[seq],
            (false, true) => PB += weights[seq],
            (false, false) => (),
        }
    }

    PA /= total_weight;
    PB /= total_weight;
    ld_obs[3] /= total_weight;

    let Pa = 1.0 - PA;
    let Pb = 1.0 - PB;
    ld_obs[2] = PA - ld_obs[3];
    ld_obs[1] = PB - ld_obs[3];
    ld_obs[0] = Pa - ld_obs[1];

    let PAB = PA * PB;
    let PAb = PA * Pb;
    let PaB = Pa * PB;
    let Pab = Pa * Pb;

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
            let a = site_set.site_symbols(first_idx);
            let a_hist = site_set.site_histogram(first_idx);
            let mut results_chunk = Vec::new();
            for second_idx in (first_idx + 1)..site_set.n_sites() {
                let b = site_set.site_symbols(second_idx);
                let b_hist = site_set.site_histogram(second_idx);
                if let Some(ld_stat) = single_weighted_ld_pair(a, a_hist, b, b_hist, weights) {
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
        let (maj, min) = hist.major_minor_symbols();
        assert_eq!(maj, Some(G));
        assert_eq!(min, Some(T));

        let hist = SymbolHistogram {
            data: [1, 9, 10, 2, 0, 0],
        };
        let (maj, min) = hist.major_minor_symbols();
        assert_eq!(maj, Some(G));
        assert_eq!(min, Some(C));
    }
}
