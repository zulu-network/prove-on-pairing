use ark_bn254::{fq::Fq, Fq12, Fq2, Fq6, G1Affine};
use ark_ec::AffineRepr;
use ark_ff::Field;
use ark_std::UniformRand;
use num_bigint::BigUint;
use num_traits::One;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::str::FromStr;

use crate::lambda_residues::LambdaResidues;
use crate::params::MODULUS;
use crate::utils::biguint_to_naf;
use crate::{
    miller_lines::MillerLines,
    utils::{fq12_to_frobenius, fq12_to_frobenius_p2, fq12_to_frobenius_p3},
};

// Inputs:
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
            let le = MillerLines::line_evaluation(alpha, bias, P);
            if i == po1 {
                println!(
                    "after line_evaluation le.x = {:?}\nle.y={:?}\nle.z={:?}\n\n",
                    le.0.to_string(),
                    le.1.to_string(),
                    le.2.to_string()
                );
            }
            f = MillerLines::mul_line_base(f, le.0, le.1, le.2);
            if i == po1 {
                println!("after mul_line_base1 f = {}\n\n", f);
            }

            if digit.pow(2) == 1 {
                let (alpha, bias) = L[lc + 1];
                let le = MillerLines::line_evaluation(alpha, bias, P);
                f = MillerLines::mul_line_base(f, le.0, le.1, le.2);
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
                let le = MillerLines::line_evaluation(alpha, bias, P);
                if i == po2 {
                    println!("k == {} le.0 = {}\n\n", k, le.0.to_string());
                    println!("k == {} le.1 = {}\n\n", k, le.1.to_string());
                    println!("k == {} le.2 = {}\n\n", k, le.2.to_string());
                }
                f = MillerLines::mul_line_base(f, le.0, le.1, le.2);
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

#[cfg(test)]
mod test {
    use super::*;
    use std::ops::{Mul, Neg};

    use ark_bn254::{Bn254, Fq12, G1Affine, G1Projective, G2Affine, G2Projective};
    use ark_ec::pairing::{Pairing, PairingOutput};
    use ark_ec::{AffineRepr, CurveGroup, Group};
    use ark_ff::{Field, One};
    use num_bigint::BigUint;
    use num_traits::{FromPrimitive, Zero};

    use crate::lambda_residues::LambdaResidues;
    use crate::optimal_ate::NativeMillerLoop;
    use crate::{dev, params};

    #[test]
    fn test_pairing_verify_native() {
        // 1. setup pairing: (p1, q1)=(p2,q2)
        //      To check (p1, q1)*(p2,-q2)=1
        let p1 = dev::g1
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();
        let p2 = dev::g1
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();

        let q1 = dev::g2
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();
        let q2 = dev::g2
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();

        // ====================================
        // ===== 2.Prover compute following data.
        // ====================================

        // 2.1 precompute lines of miller_loop
        let l1 = MillerLines::precompute_lines(
            G2Projective::from(q1),
            params::E.clone(),
            params::LAMBDA.clone(),
        );
        let l2 = MillerLines::precompute_lines(
            G2Projective::from(q2.neg()),
            params::E.clone(),
            params::LAMBDA.clone(),
        );

        // 2.2 precompute witness
        // TODO: meet error when replacing it with compute.
        let f1 = NativeMillerLoop::miller_loop(G1Projective::from(p1), G2Projective::from(q1));
        let f2 =
            NativeMillerLoop::miller_loop(G1Projective::from(p2), G2Projective::from(q2).neg());

        println!("f1: {:?}", f1.to_string());
        println!("f2: {:?}", f2.to_string());
        // 2.3 precompute c,wi
        let witness = LambdaResidues::finding_c(f1.mul(f2));
        let c_inv = witness.c.inverse().unwrap();
        println!("c_inv: {:?}", c_inv.to_string());
        // ====================================
        // ===== 2.Prover compute following data.
        // ====================================
        let verify_res = verify_pairings(
            vec![p1, p2],
            &[l1, l2],
            params::E.clone(),
            witness.c,
            c_inv,
            witness.wi,
        );
        assert_eq!(verify_res, Fq12::ONE);
        println!("========Successfully");
    }

    #[ignore]
    #[test]
    fn test_pairing_verify_full() {
        // 1. setup pairing: (p1, q1)=(p2,q2)
        //      To check (p1, q1)*(p2,-q2)=1
        let p1 = dev::g1
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();
        let p2 = dev::g1
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();

        let q1 = dev::g2
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();
        let q2 = dev::g2
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();

        // ====================================
        // ===== 2.Prover compute following data.
        // ====================================

        // 2.1 precompute lines of miller_loop
        let l1 = MillerLines::precompute_lines(
            G2Projective::from(q1),
            params::E.clone(),
            params::LAMBDA.clone(),
        );
        let l2 = MillerLines::precompute_lines(
            G2Projective::from(q2.neg()),
            params::E.clone(),
            params::LAMBDA.clone(),
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
        let witness = LambdaResidues::finding_c(f1.mul(f2));
        let c_inv = witness.c.inverse().unwrap();
        let verify_res = crate::pairing_verify::verify_pairings(
            vec![p1, p2],
            &[l1, l2],
            params::E.clone(),
            witness.c,
            c_inv,
            witness.wi,
        );
        assert_eq!(verify_res, Fq12::ONE);
        println!("========Successfully");
    }

    #[test]
    fn test_multi_pairing() {
        type E = Bn254;

        // 8 * 10 = 2 * 3 + 4 * 5 + 6 * 9
        let zero_g1: <E as Pairing>::G1 = <<E as Pairing>::G1Affine as AffineRepr>::zero().into();
        let one_g1 = <<E as Pairing>::G1 as Group>::generator();
        let one_g2 = <<E as Pairing>::G2 as Group>::generator();
        let a = zero_g1 - one_g1.mul_bigint([8]);
        let b = one_g2.mul_bigint([10]);
        let c = one_g1.mul_bigint([2]);
        let d = one_g2.mul_bigint([3]);
        let e = one_g1.mul_bigint([4]);
        let f = one_g2.mul_bigint([5]);
        let g = one_g1.mul_bigint([6]);
        let h = one_g2.mul_bigint([9]);
        let ans1 = E::multi_pairing(&[a, c, e, g], &[b, d, f, h]);
        let ans2 = PairingOutput::<E>::zero();
        assert_eq!(ans1, ans2);

        // let a_affine: G1Affine = a.into();
        // let b_affine: G2Affine = b.into();
        // let c_affine: G1Affine = c.into();
        // let d_affine: G2Affine = d.into();
        // let e_affine: G1Affine = e.into();
        // let f_affine: G2Affine = f.into();
        // let g_affine: G1Affine = g.into();
        // let h_affine: G2Affine = h.into();

        // dbg!(a_affine);
        // dbg!(b_affine);
        // dbg!(c_affine);
        // dbg!(d_affine);
        // dbg!(e_affine);
        // dbg!(f_affine);
        // dbg!(g_affine);
        // dbg!(h_affine);
    }
}
