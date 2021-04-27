use std::{
    fs::File,
    io::{BufRead, BufReader},
    marker::PhantomData,
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, FromPrimitive)]
#[repr(u8)]
pub enum MajMin {
    Major = 0,
    Minor,
    Other,
}

const SYMBOL_HISTOGRAM_CAPACITY: usize = 16;

pub unsafe trait Histogrammable: Sized + Copy {
    const VARIANT_COUNT: usize;
    fn to_index(&self) -> usize;
}

unsafe impl Histogrammable for Symbol {
    const VARIANT_COUNT: usize = 6;

    fn to_index(&self) -> usize {
        *self as usize
    }
}

unsafe impl Histogrammable for MajMin {
    const VARIANT_COUNT: usize = 3;

    fn to_index(&self) -> usize {
        *self as usize
    }
}

pub struct SymbolHistogram<T: Histogrammable> {
    data: [usize; SYMBOL_HISTOGRAM_CAPACITY],
    _phantom: PhantomData<T>,
}

impl<T: Histogrammable> Index<T> for SymbolHistogram<T> {
    type Output = usize;

    fn index(&self, sym: T) -> &Self::Output {
        &self.data[sym.to_index()]
    }
}

impl<T: Histogrammable> IndexMut<T> for SymbolHistogram<T> {
    fn index_mut(&mut self, sym: T) -> &mut Self::Output {
        &mut self.data[sym.to_index()]
    }
}

impl<T: Histogrammable> SymbolHistogram<T> {
    fn zero() -> Self {
        Self {
            data: [0; SYMBOL_HISTOGRAM_CAPACITY],
            _phantom: PhantomData,
        }
    }

    fn from_slice(symbols: &[T]) -> Self {
        let mut hist = Self::zero();
        for symbol in symbols {
            hist[*symbol] += 1;
        }
        hist
    }

    fn total(&self) -> usize {
        self.data.iter().sum()
    }
}

impl SymbolHistogram<Symbol> {
    fn acgt(&self) -> usize {
        use Symbol::*;
        self[A] + self[C] + self[G] + self[T]
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
}

pub struct MultiSequence {
    pub source: MultiSequenceSource,
    pub sequences: Vec<Sequence>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SiteIndex {
    pub index: usize,
    pub strand: usize,
}

/// Stores symbol data such that all data for a given site is contiguous
pub struct SiteSet<T> {
    n_sites: usize,
    n_seqs: usize,

    buffer: Vec<T>,

    /// Vector of length n_sites, mapping a site index in this set to some other index
    ///
    /// The site with index N in this set has an index `site_map[N]`
    site_map: Vec<SiteIndex>,
}

impl<T> Index<usize> for SiteSet<T> {
    type Output = [T];

    fn index(&self, index: usize) -> &Self::Output {
        let start = index * self.n_seqs;
        &self.buffer[start..(start + self.n_seqs)]
    }
}

impl<T: Clone> SiteSet<T> {
    pub fn filter_by(&self, filter_func: impl Fn(&[T]) -> bool) -> Self {
        let mut new_buffer = Vec::with_capacity(self.buffer.len());
        let mut site_map = Vec::new();

        for site in 0..self.n_sites {
            let site_slice = &self[site];
            if filter_func(site_slice) {
                site_map.push(self.site_map[site]);
                new_buffer.extend_from_slice(site_slice);
            }
        }

        Self {
            n_sites: site_map.len(),
            n_seqs: self.n_seqs,
            buffer: new_buffer,
            site_map: site_map,
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

    pub fn site_index(&self, idx: usize) -> SiteIndex {
        self.site_map[idx]
    }
}

impl SiteSet<Symbol> {
    pub fn from_multiseq(ms: &MultiSequence) -> SiteSet<Symbol> {
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
        
        let site_map = (0..n_sites)
            .map(|x| SiteIndex { index: x, strand: 0 })
            .collect();

        Self {
            n_sites,
            n_seqs,
            buffer,
            site_map,
        }
    }

    pub fn to_maj_min(&self) -> SiteSet<MajMin> {
        let mut new_buffer = vec![MajMin::Other; self.buffer.len()];
        for site in 0..self.n_sites() {
            let hist = SymbolHistogram::from_slice(&self[site]);
            let (maj, min) = match hist.major_minor_symbols() {
                (Some(maj), Some(min)) => (maj, min),
                _ => continue,
            };

            for seq in 0..self.n_seqs() {
                new_buffer[site * self.n_seqs() + seq] = match self[site][seq] {
                    x if x == maj => MajMin::Major,
                    x if x == min => MajMin::Minor,
                    _ => MajMin::Other,
                };
            }
        }

        SiteSet {
            n_sites: self.n_sites(),
            n_seqs: self.n_seqs(),
            buffer: new_buffer,
            site_map: self.site_map.clone(),
        }
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

/// Given a slice of all symbols in a site, should the site be considered for LD computations
pub fn is_site_of_interest(site: &[Symbol], min_acgt: f32, min_minor: f32, max_minor: f32) -> bool {
    let hist = SymbolHistogram::from_slice(site);

    let acgt_count = hist.acgt() as f32;
    let total_count = hist.total() as f32;

    if (acgt_count / total_count) <= min_acgt {
        // There aren't enough ACGT symbols at this site
        false
    } else {
        let (major_sym, minor_sym) = match hist.major_minor_symbols() {
            (Some(maj), Some(min)) => (maj, min),
            _ => return false,
        };

        let maj_count = hist[major_sym] as f32;
        let min_count = hist[minor_sym] as f32;

        // let minor_frac = min_count / acgt_count;
        let minor_frac = min_count / (min_count + maj_count);
        // let minor_frac = (acgt_count - min_count) / acgt_count;

        if minor_frac < min_minor || minor_frac > max_minor {
            // The fraction of minor symbols is either too low or too high
            false
        } else {
            true
        }
    }
}

pub trait Henikoffable: Histogrammable + Eq + 'static {
    fn exclude_syms() -> &'static [Self];
}

impl Henikoffable for Symbol {
    fn exclude_syms() -> &'static [Self] {
        &[Self::Unknown]
    }
}

impl Henikoffable for MajMin {
    fn exclude_syms() -> &'static [Self] {
        &[Self::Other]
    }
}

pub fn henikoff_weights<T: Henikoffable>(data: &SiteSet<T>) -> Vec<f32> {
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

pub fn henikoff_site_contributions<T: Henikoffable>(site: &[T], contributions: &mut [f32]) {
    fn is_excluded<T: Henikoffable>(s: &T) -> bool {
        for x in T::exclude_syms().iter() {
            if s == x {
                return true;
            }
        }

        false
    }

    let hist = SymbolHistogram::from_slice(site);

    let nuc_count = {
        let mut x = hist.total();
        for exclude in T::exclude_syms() {
            x -= hist[*exclude];
        }
        x as f32
    };

    let mut total_site_contrib = 0f32;

    for (idx, sym) in site.iter().enumerate() {
        if !is_excluded(sym) {
            contributions[idx] = 1f32 / (nuc_count * hist[*sym] as f32);
            total_site_contrib += contributions[idx];
        }
    }

    let mean_site_contrib = total_site_contrib / nuc_count;

    for (idx, sym) in site.iter().enumerate() {
        if is_excluded(sym) {
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
pub fn single_weighted_ld_pair(a: &[MajMin], b: &[MajMin], weights: &[f32]) -> Option<LdStats> {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), weights.len());

    let good_mask = a
        .iter()
        .zip(b.iter())
        .map(|(&a, &b)| a != MajMin::Other && b != MajMin::Other)
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

        if a[seq] == MajMin::Major {
            PA += weights[seq];
        } else {
            // good_mask has already filtered down the sequences to ones where
            // a is either a_maj_sym or a_min_sym -> if a isn't a_maj_sym it
            // must be a_min_sym.
            Pa += weights[seq];
        }

        if b[seq] == MajMin::Major {
            PB += weights[seq];
        } else {
            Pb += weights[seq];
        }

        match (a[seq] == MajMin::Major, b[seq] == MajMin::Major) {
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

struct PairData<T> {
    first_idx: SiteIndex,
    second_idx: SiteIndex,
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
    type Item = (SiteIndex, SiteIndex, &'a T);

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
    siteset: &SiteSet<MajMin>,
    weights: &[f32],
    r2_threshold: f32,
    mut progress_report: impl FnMut(usize) + Send,
) -> PairStore<LdStats> {
    progress_report(0);

    let computed_pair_count = AtomicUsize::new(0);
    let progress_report = Mutex::new(progress_report);

    let data_chunks = (0..(siteset.n_sites() - 1))
        .into_par_iter()
        .map(|first_idx| {
            let a = &siteset[first_idx];
            let mut results_chunk = Vec::new();
            for second_idx in (first_idx + 1)..siteset.n_sites() {
                let b = &siteset[second_idx];
                if let Some(ld_stat) = single_weighted_ld_pair(a, b, weights) {
                    if ld_stat.r2 > r2_threshold {
                        results_chunk.push(PairData {
                            first_idx: siteset.site_index(first_idx),
                            second_idx: siteset.site_index(second_idx),
                            data: ld_stat,
                        });
                    }
                }
            }

            let computed =
                computed_pair_count.fetch_add(siteset.n_sites() - first_idx - 1, Ordering::Relaxed);
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
        let hist = SymbolHistogram::from_slice(&a);

        assert_eq!(hist[A], 3);
        assert_eq!(hist[C], 2);
        assert_eq!(hist[G], 1);
        assert_eq!(hist[T], 4);
    }

    #[test]
    fn test_hist_major_minor() {
        use Symbol::*;

        let mut hist = SymbolHistogram::zero();
        hist[A] = 0;
        hist[C] = 1;
        hist[G] = 10;
        hist[T] = 2;
        let (maj, min) = hist.major_minor_symbols();
        assert_eq!(maj, Some(G));
        assert_eq!(min, Some(T));

        let mut hist = SymbolHistogram::zero();
        hist[A] = 1;
        hist[C] = 9;
        hist[G] = 10;
        hist[T] = 2;
        let (maj, min) = hist.major_minor_symbols();
        assert_eq!(maj, Some(G));
        assert_eq!(min, Some(C));
    }
}
