use bulletproofs::LinearProof;
use bulletproofs::{BulletproofGens, PedersenGens};
use core::iter;
use criterion::{criterion_group, criterion_main, Criterion};
use curve25519_dalek::traits::VartimeMultiscalarMul;
use curve25519_dalek::{ristretto::RistrettoPoint, scalar::Scalar};
use merlin::Transcript;

criterion_main!(bench);
criterion_group!(bench, benches);

fn u8_to_scalar(x: u8) -> Scalar {
    let mut bits = [0; 32];
    bits[0] = x;
    Scalar::from_bits(bits)
}

fn inner_product(a: &[Scalar], b: &[Scalar]) -> Scalar {
    let mut out = Scalar::zero();
    if a.len() != b.len() {
        panic!("inner_product(a,b): lengths of vectors do not match");
    }
    for i in 0..a.len() {
        out += a[i] * b[i];
    }
    out
}

fn prove<const LOG_N: usize>(criterion: &mut Criterion) {
    let vector_len = 1 << LOG_N;
    let mut rng = rand::thread_rng();
    let bp_gens = BulletproofGens::new(vector_len, 1);
    let G: Vec<RistrettoPoint> = bp_gens.share(0).G(vector_len).cloned().collect();

    let pedersen_gens = PedersenGens::default();
    let F = pedersen_gens.B;
    let B = pedersen_gens.B_blinding;

    let a: Vec<_> = (0..vector_len).map(|_| Scalar::random(&mut rng)).collect();
    let point = (0..LOG_N)
        .map(|_| Scalar::random(&mut rng))
        .collect::<Vec<_>>();
    let mut b = vec![u8_to_scalar(1)];
    for i in point {
        let len = b.len();
        for j in 0..len {
            b.push(b[j] * i);
        }
    }
    let mut transcript = Transcript::new(b"LinearProofBenchmark");

    // C = <a, G> + r * B + <a, b> * F
    let r = Scalar::random(&mut rng);
    let c = inner_product(&a, &b);
    let C = RistrettoPoint::vartime_multiscalar_mul(
        a.iter().chain(iter::once(&r)).chain(iter::once(&c)),
        G.iter().chain(iter::once(&B)).chain(iter::once(&F)),
    )
    .compress();

    let mut group = criterion.benchmark_group("prove");
    group.sample_size(10);
    group.bench_function("open point", move |bencher| {
        bencher.iter(|| {
            LinearProof::create(
                &mut transcript,
                &mut rng,
                &C,
                r,
                a.clone(),
                b.clone(),
                G.clone(),
                &F,
                &B,
            )
            .unwrap();
        })
    });
}

fn benches(c: &mut Criterion) {
    prove::<16>(c);
    verify::<16>(c);
}

fn verify<const LOG_N: usize>(criterion: &mut Criterion) {
    let vector_len = 1 << LOG_N;
    let bp_gens = BulletproofGens::new(vector_len, 1);
    let mut rng = rand::thread_rng();

    // Calls `.G()` on generators, which should be a pub(crate) function only.
    // For now, make that function public so it can be accessed from benches.
    // We can't simply use bp_gens directly because we don't need the H generators.
    let G: Vec<RistrettoPoint> = bp_gens.share(0).G(vector_len).cloned().collect();
    let pedersen_gens = PedersenGens::default();
    let F = pedersen_gens.B;
    let B = pedersen_gens.B_blinding;

    let point = (0..LOG_N)
        .map(|_| Scalar::random(&mut rng))
        .collect::<Vec<_>>();
    let mut b = vec![u8_to_scalar(1)];
    for i in &point {
        let len = b.len();
        for j in 0..len {
            b.push(b[j] * i);
        }
    }

    // Generate the proof in its own scope to prevent reuse of
    // prover variables by the verifier
    let (proof, C) = {
        // a and b are the vectors for which we want to prove c = <a,b>
        let a: Vec<_> = (0..vector_len).map(|_| Scalar::random(&mut rng)).collect();

        let mut transcript = Transcript::new(b"LinearProofBenchmark");

        // C = <a, G> + r * B + <a, b> * F
        let r = Scalar::random(&mut rng);
        let c = inner_product(&a, &b);
        let C = RistrettoPoint::vartime_multiscalar_mul(
            a.iter().chain(iter::once(&r)).chain(iter::once(&c)),
            G.iter().chain(iter::once(&B)).chain(iter::once(&F)),
        )
        .compress();

        let proof = LinearProof::create(
            &mut transcript,
            &mut rng,
            &C,
            r,
            a.clone(),
            b.clone(),
            G.clone(),
            &F,
            &B,
        )
        .unwrap();

        (proof, C)
    };

    let mut group = criterion.benchmark_group("prove");
    group.sample_size(10);
    group.bench_function("verify point value", move |bencher| {
        bencher.iter(|| {
            let mut verifier_transcript = Transcript::new(b"LinearProofBenchmark");
            let mut b = vec![u8_to_scalar(1)];
            for i in &point {
                let len = b.len();
                for j in 0..len {
                    b.push(b[j] * i);
                }
            }
            proof
                .verify(&mut verifier_transcript, &C, &G, &F, &B, b)
                .unwrap();
        })
    });
}
