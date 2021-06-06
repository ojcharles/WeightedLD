
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rand::prelude::*;
use num_traits::FromPrimitive;

use weighted_ld::*;

fn make_sequence(len: usize, missing_frac: f32, major_frac: f32) -> Vec<Symbol> {
    let mut r = thread_rng();
    
    let major_sym: Symbol = Symbol::from_u8(r.gen_range(0..4)).unwrap();
    let mut minor_sym: Symbol;
    loop {
        minor_sym = Symbol::from_u8(r.gen_range(0..4)).unwrap();
        if minor_sym != major_sym { break }
    }
    
    (0..len).map(|_| {
        let num = r.gen_range(0f32..1f32);
        if num < missing_frac {
            Symbol::Missing
        } else if num < (missing_frac + major_frac) {
            major_sym
        } else {
            minor_sym
        }
    }).collect()
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("weighted_pair_ld");
    let mut r = thread_rng();

    for len in [10, 50, 100, 250, 500, 1000].iter() {
        let seq_a = make_sequence(*len, 0.1, 0.6);
        let seq_b = make_sequence(*len, 0.1, 0.6);
        let a_hist = SymbolHistogram::from_slice(&seq_a);
        let b_hist = SymbolHistogram::from_slice(&seq_b);

        let mut weights = vec![0f32; *len];
        // Defaults to uniforms in the range 0..1
        r.fill(&mut weights[..]);
        
        group.throughput(Throughput::Elements(*len as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(len),
            &(seq_a, seq_b, weights),
            |b, (seq_a, seq_b, weights)| {
                b.iter(|| single_weighted_ld_pair(seq_a, &a_hist, seq_b, &b_hist, weights));
            });
    }
    
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);