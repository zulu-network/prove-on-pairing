use ark_bn254::{fq::Fq, Fq12, Fq2, Fq6, Fr, G1Affine};
use ark_ec::AffineRepr;
use ark_ff::{Field, PrimeField};
use ark_std::UniformRand;
use num_bigint::BigUint;
use num_traits::{Num, One};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::str::FromStr;

use crate::constant::MODULUS;
use crate::final_exp_verifier::tonelli_shanks_cubic;
use crate::utils::biguint_to_naf;
use crate::{
    fields::{fq12_to_frobenius, fq12_to_frobenius_p2, fq12_to_frobenius_p3},
    miller_loop_verify::line_evaluation,
    optimal_ate::mul_line_base,
    utils::to_naf,
};

// Inputs:
//
//
// verify c^lambda = f * wi, namely c_inv^lambda * f * wi = 1
pub fn verify_pairings(
    eval_points: Vec<G1Affine>,
    lines: &[Vec<(Fq2, Fq2)>; 2],
    e: BigUint,
    c: Fq12,
    c_inv: Fq12,
    wi: Fq12,
) -> Fq12 {
    assert_eq!(eval_points.len(), lines.len());
    assert_eq!(c * c_inv, Fq12::ONE);

    let mut lc = 0_usize;
    let mut f = c_inv;
    let mut naf_digits = biguint_to_naf(e);
    naf_digits.reverse();
    naf_digits.remove(0);

    let po1 = 9999;
    for (i, digit) in naf_digits.into_iter().enumerate() {
        if i == po1 {
            println!("before f.square, f = {}\n\n", f);
        }
        f = f.square();
        if i == po1 {
            println!("after f.square, f = {}\n\n", f);
        }
        // update c^lambda
        if digit.pow(2) == 1 {
            f = if digit == 1 { f * c_inv } else { f * c };
            if i == po1 {
                dbg!(digit);
                println!("update c^lambda f = {}\n\n", f);
            }
        } else {
            if i == po1 {
                println!("update c^lambda f = {}\n\n", f);
            }
        }
        for (&P, L) in eval_points.iter().zip(lines) {
            let (alpha, bias) = L[lc];
            if i == po1 {
                println!("alpha = {:?}\n", alpha.to_string());
                println!("bias = {:?}\n\n", bias.to_string());
            }
            let le = line_evaluation(alpha, bias, P);
            if i == po1 {
                println!(
                    "after line_evaluation le.x = {:?}\nle.y={:?}\nle.z={:?}\n\n",
                    le.0.to_string(),
                    le.1.to_string(),
                    le.2.to_string()
                );
            }
            f = mul_line_base(f, le.0, le.1, le.2);
            if i == po1 {
                println!("after mul_line_base1 f = {}\n\n", f);
            }

            if digit.pow(2) == 1 {
                let (alpha, bias) = L[lc + 1];
                let le = line_evaluation(alpha, bias, P);
                f = mul_line_base(f, le.0, le.1, le.2);
                if i == po1 {
                    println!("after mul_line_base2 f = {}\n\n", f);
                }
            } else {
                if i == po1 {
                    println!("after mul_line_base2 f = {}\n\n", f);
                }
            }
        }
        if i == po1 {
            println!("after for loop f = {}", f);
        }
        lc = if digit == 0 { lc + 1 } else { lc + 2 };
    }

    // update c^lambda
    println!("before fq12_to_frobenius f = {}\n\n", f);
    f = f * fq12_to_frobenius(c_inv) * fq12_to_frobenius_p2(c) * fq12_to_frobenius_p3(c_inv);
    println!(
        "fq12_to_frobenius_p2(c) = {}\n\n",
        fq12_to_frobenius_p2(c).to_string()
    );
    println!(
        "fq12_to_frobenius_p3(c_inv) = {}\n\n",
        fq12_to_frobenius_p3(c_inv).to_string()
    );
    println!("after fq12_to_frobenius f = {}\n\n", f);
    // update the scalar
    f = f * wi;
    println!("after update the scalar f = {}\n\n", f);

    let po2 = 9999;
    // frobenius map part, p - p^2 + p^3
    for (i, (P, L)) in eval_points.into_iter().zip(lines).enumerate() {
        println!("frobenius map part i = {}\n", i);
        for k in 0..3 {
            let (alpha, bias) = L[lc + k];
            println!("alpha = {:?}\n\n", alpha.to_string());
            println!("bias = {:?}\n\n", bias.to_string());
            if k == 2 {
                let eval = Fq12::new(
                    Fq6::new(
                        Fq2::new(P.x().expect("failed to read x").to_owned(), Fq::ZERO),
                        -bias,
                        Fq2::ZERO,
                    ),
                    Fq6::ZERO,
                );
                if i == po2 {
                    println!("k = {} eval = {}\n\n", k, eval.to_string());
                }
                f = f * eval;
                if i == po2 {
                    println!("k = {} f = {}\n\n", k, f);
                }
            } else {
                let le = line_evaluation(alpha, bias, P);
                if i == po2 {
                    println!("k == {} le.0 = {}\n\n", k, le.0.to_string());
                    println!("k == {} le.1 = {}\n\n", k, le.1.to_string());
                    println!("k == {} le.2 = {}\n\n", k, le.2.to_string());
                }
                f = mul_line_base(f, le.0, le.1, le.2);
                if i == po2 {
                    println!("k = {} f = {}\n\n", k, f);
                }
            }
        }
    }
    lc = lc + 3;
    assert_eq!(lc, lines[0].len());
    assert_eq!(f, Fq12::ONE);
    f
}

// Computing λ residues over BN curve
// Input:
//      f: output of a Miller loop
// Output:
//      c and wi,
//      satisfying c^lambda = f * wi
//
// Ref: Algorithm 5 of [On Proving Pairings](https://eprint.iacr.org/2024/640.pdf)
fn compute_lambda_residues(f: ark_bn254::Fq12) -> (ark_bn254::Fq12, ark_bn254::Fq12) {
    let p = MODULUS;
    let r = BigUint::from_str(
        "21888242871839275222246405745257275088548364400416034343698204186575808495617",
    )
    .unwrap();
    let lambda = BigUint::from_str(
            "10486551571378427818905133077457505975146652579011797175399169355881771981095211883813744499745558409789005132135496770941292989421431235276221147148858384772096778432243207188878598198850276842458913349817007302752534892127325269"
        ).unwrap();
    let s = 3_u32;
    let exp = p.pow(12_u32) - 1_u32;
    let h = &exp / &r;
    let t = &exp / 3_u32.pow(s);
    let k = (&t + 1_u32) / 3_u32;
    let m = &lambda / &r;
    let d = 3_u32;
    let mm = &m / d;

    let mut prng = ChaCha20Rng::seed_from_u64(0);
    let cofactor_cubic = 3_u32.pow(s - 1) * &t;

    // make f is r-th residue, but it's not cubic residue
    assert_eq!(f.pow(h.to_u64_digits()), ark_bn254::Fq12::ONE);
    assert_ne!(f.pow(cofactor_cubic.to_u64_digits()), ark_bn254::Fq12::ONE);

    // sample a proper scalar w which is cubic non-residue
    let w = {
        let (mut w, mut z) = (ark_bn254::Fq12::ONE, ark_bn254::Fq12::ONE);
        while w == ark_bn254::Fq12::ONE {
            // choose z which is 3-th non-residue
            let mut legendre = ark_bn254::Fq12::ONE;
            while legendre == ark_bn254::Fq12::ONE {
                z = ark_bn254::Fq12::rand(&mut prng);
                legendre = z.pow(cofactor_cubic.to_u64_digits());
            }
            // obtain w which is t-th power of z
            w = z.pow(t.to_u64_digits());
        }
        w
    };
    // make sure 27-th root w, is 3-th non-residue and r-th residue
    assert_ne!(w.pow(cofactor_cubic.to_u64_digits()), ark_bn254::Fq12::ONE);
    assert_eq!(w.pow(h.to_u64_digits()), ark_bn254::Fq12::ONE);

    // just two option, w and w^2, since w^3 must be cubic residue, leading f*w^3 must not be cubic residue
    let mut wi = w;
    if (f * wi).pow(cofactor_cubic.to_u64_digits()) != ark_bn254::Fq12::ONE {
        assert_eq!(
            (f * w * w).pow(cofactor_cubic.to_u64_digits()),
            ark_bn254::Fq12::ONE
        );
        wi = w * w;
    }
    assert_eq!(wi.pow(h.to_u64_digits()), ark_bn254::Fq12::ONE);

    assert_eq!(lambda, &d * &mm * &r);
    // f1 is scaled f
    let f1 = f * wi;

    // r-th root of f1, say f2
    let r_inv = r.modinv(&h).unwrap();
    assert_ne!(r_inv, BigUint::one());
    let f2 = f1.pow(r_inv.to_u64_digits());
    assert_ne!(f2, ark_bn254::Fq12::ONE);

    // m'-th root of f, say f3
    let mm_inv = mm.modinv(&(r * h)).unwrap();
    assert_ne!(mm_inv, BigUint::one());
    let f3 = f2.pow(mm_inv.to_u64_digits());
    assert_eq!(f3.pow(cofactor_cubic.to_u64_digits()), ark_bn254::Fq12::ONE);
    assert_ne!(f3, ark_bn254::Fq12::ONE);

    // d-th (cubic) root, say c
    let c = tonelli_shanks_cubic(f3, w, s, t, k);
    assert_ne!(c, ark_bn254::Fq12::ONE);
    assert_eq!(c.pow(lambda.to_u64_digits()), f * wi);

    (c, wi)
}

#[cfg(test)]
mod test {
    use super::*;
    use std::ops::{Mul, Neg};
    use std::{
        fs::File,
        io::{self, Read},
        str::FromStr,
    };

    use ark_bn254::{Bn254, Fq, Fq12, Fq2, Fq6, G1Affine, G1Projective, G2Affine, G2Projective};
    use ark_ec::pairing::Pairing;
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ff::{Field, MontFp, One};
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;

    use crate::constant;
    use crate::constant::LAMBDA;
    use crate::line_precompute::line_function;
    use crate::optimal_ate::{line_func_add, miller_loop};

    #[test]
    fn test_pairing_verify_1() {
        // 1. setup pairing: (p1, q1)=(p2,q2)
        //      To check (p1, q1)*(p2,-q2)=1
        let p1 = constant::g1
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();
        let p2 = constant::g1
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();

        let q1 = constant::g2
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();
        let q2 = constant::g2
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();

        // ====================================
        // ===== 2.Prover compute following data.
        // ====================================

        // 2.1 precompute lines of miller_loop
        let l1 = line_function(
            G2Projective::from(q1),
            constant::E.clone(),
            constant::LAMBDA.clone(),
        );
        let l2 = line_function(
            G2Projective::from(q2.neg()),
            constant::E.clone(),
            constant::LAMBDA.clone(),
        );

        // 2.2 precompute witness
        let cy = Fq6::new(
            Fq2::new(
                MontFp!(
                    "1398074605395234385742746982635607890999901394449920787952679263329618900423"
                ),
                MontFp!(
                    "14863623346637387528754574561712611762280560946730003365167624628783031806851"
                ),
            ),
            Fq2::new(
                MontFp!(
                    "3943322538927771188455999564370373075332486995274784857622626372657519173057"
                ),
                MontFp!(
                    "3097735610337843900368979344756038697740814458812421595304050335303356161701"
                ),
            ),
            Fq2::new(
                MontFp!(
                    "15535495333518547635388053375087086844764140428602229212987121796906694947701"
                ),
                MontFp!(
                    "17900230360671217726773770977461976340138958556757112953925160338294259755451"
                ),
            ),
        );
        let cx = Fq6::new(
            Fq2::new(
                MontFp!(
                    "4667318771591831111896913998599388465660409602488670805853346418592518763334"
                ),
                MontFp!(
                    "6713632336462274726021105432467423937237553861088593973094350524276003798886"
                ),
            ),
            Fq2::new(
                MontFp!(
                    "21081939668623356886869294872217712850424216414807721597624545100857744884557"
                ),
                MontFp!(
                    "16588937522064281852235800901297894214318954372559985137225140148647705064928"
                ),
            ),
            Fq2::new(
                MontFp!(
                    "19903885137215268899583999685812289845712323958279995809868941246364604339090"
                ),
                MontFp!(
                    "15864462405647635712161095062692932104398676326285281473670528335466546619372"
                ),
            ),
        );
        let c = Fq12::new(cx, cy);

        let c_inv = c.inverse().unwrap();

        let wi_x = Fq6::new(
            Fq2::ZERO,
            Fq2::new(
                MontFp!(
                    "1552599724427260165122619260710620891632201030296138444657693101145235169998"
                ),
                MontFp!(
                    "6872310070170355934577608915679474652968694353285292038893436960198570441414"
                ),
            ),
            Fq2::ZERO,
        );
        let wi = Fq12::new(wi_x, Fq6::ZERO);

        let verify_res =
            verify_pairings(vec![p1, p2], &[l1, l2], constant::E.clone(), c, c_inv, wi);
        assert_eq!(verify_res, Fq12::ONE);
        println!("========Successfully");
    }

    #[test]
    fn test_pairing_verify_native_full() {
        // 1. setup pairing: (p1, q1)=(p2,q2)
        //      To check (p1, q1)*(p2,-q2)=1
        let p1 = constant::g1
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();
        let p2 = constant::g1
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();

        let q1 = constant::g2
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();
        let q2 = constant::g2
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();

        // ====================================
        // ===== 2.Prover compute following data.
        // ====================================

        // 2.1 precompute lines of miller_loop
        let l1 = line_function(
            G2Projective::from(q1),
            constant::E.clone(),
            constant::LAMBDA.clone(),
        );
        let l2 = line_function(
            G2Projective::from(q2.neg()),
            constant::E.clone(),
            constant::LAMBDA.clone(),
        );

        // 2.2 precompute witness
        // TODO: meet error when replacing it with compute.
        let f1 = miller_loop(G1Projective::from(p1), G2Projective::from(q1));
        let f2 = miller_loop(G1Projective::from(p2), G2Projective::from(q2));

        // 2.3 precompute c,wi
        let (c, wi) = compute_lambda_residues(f1.mul(f2));
        let c_inv = c.inverse().unwrap();
        let verify_res = crate::pairing_verify::verify_pairings(
            vec![p1, p2],
            &[l1, l2],
            constant::E.clone(),
            c,
            c_inv,
            wi,
        );
        assert_eq!(verify_res, Fq12::ONE);
        println!("========Successfully");
    }

    #[ignore]
    #[test]
    fn test_pairing_verify_full() {
        // 1. setup pairing: (p1, q1)=(p2,q2)
        //      To check (p1, q1)*(p2,-q2)=1
        let p1 = constant::g1
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();
        let p2 = constant::g1
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();

        let q1 = constant::g2
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();
        let q2 = constant::g2
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();

        // ====================================
        // ===== 2.Prover compute following data.
        // ====================================

        // 2.1 precompute lines of miller_loop
        let l1 = line_function(
            G2Projective::from(q1),
            constant::E.clone(),
            constant::LAMBDA.clone(),
        );
        let l2 = line_function(
            G2Projective::from(q2.neg()),
            constant::E.clone(),
            constant::LAMBDA.clone(),
        );

        // 2.2 precompute witness
        // TODO: meet error when replacing it with compute.
        let f1 = Bn254::miller_loop(
            <G1Affine as Into<<Bn254 as ark_ec::pairing::Pairing>::G1Prepared>>::into(p1),
            <G2Affine as Into<<Bn254 as ark_ec::pairing::Pairing>::G2Prepared>>::into(q1),
        );
        let f2 = Bn254::miller_loop(
            <G1Affine as Into<<Bn254 as ark_ec::pairing::Pairing>::G1Prepared>>::into(p2),
            <G2Affine as Into<<Bn254 as ark_ec::pairing::Pairing>::G2Prepared>>::into(q2.neg()),
        );

        let (f1, f2) = (f1.0, f2.0);

        // 2.3 precompute c,wi
        let (c, wi) = compute_lambda_residues(f1.mul(f2));
        let c_inv = c.inverse().unwrap();
        let verify_res = crate::pairing_verify::verify_pairings(
            vec![p1, p2],
            &[l1, l2],
            constant::E.clone(),
            c,
            c_inv,
            wi,
        );
        assert_eq!(verify_res, Fq12::ONE);
        println!("========Successfully");
    }
}
