use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::prelude::*;

use weighted_ld::*;

fn make_sequence(len: usize, missing_frac: f32, major_frac: f32) -> Vec<MajMin> {
    let mut r = thread_rng();

    (0..len)
        .map(|_| {
            let num = r.gen_range(0f32..1f32);
            if num < missing_frac {
                MajMin::Other
            } else if num < (missing_frac + major_frac) {
                MajMin::Major
            } else {
                MajMin::Minor
            }
        })
        .collect()
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("weighted_pair_ld");
    let mut r = thread_rng();

    for len in [10, 50, 100, 250, 500, 1000].iter() {
        let seq_a = make_sequence(*len, 0.1, 0.6);
        let seq_b = make_sequence(*len, 0.1, 0.6);

        let mut weights = vec![0f32; *len];
        // Defaults to uniforms in the range 0..1
        r.fill(&mut weights[..]);

        group.throughput(Throughput::Elements(*len as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(len),
            &(seq_a, seq_b, weights),
            |b, (seq_a, seq_b, weights)| {
                b.iter(|| single_weighted_ld_pair(seq_a, seq_b, weights));
            },
        );
    }

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
