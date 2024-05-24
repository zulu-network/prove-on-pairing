use ark_bn254::{fq::Fq, Fq12, Fq2, Fq6, G1Affine};
use ark_ec::AffineRepr;
use ark_ff::Field;
use num_bigint::BigUint;

use crate::poc::{
    miller_lines::MillerLines,
    utils::{fq12_to_frobenius, fq12_to_frobenius_p2, fq12_to_frobenius_p3},
};
use crate::utils::biguint_to_naf;

// To verify e(P1,Q1)=e(P2,Q2), which equal e(P1,Q1)*(P2,-Q2)=1
//
// params:
//  @eval_points: Pi points
//  @lines: precompute miller lines for Qi. Only support fixed Qi.
//  @c: c^lambda = f*w^i
//  @c_inv: inverse of c
//
// verify c^lambda = f * wi, namely c_inv^lambda * f * wi = 1
pub fn dual_miller_loop_with_c_wi(
    eval_points: Vec<G1Affine>,
    lines: &[Vec<(Fq2, Fq2)>],
    e: BigUint,
    c: Fq12,
    c_inv: Fq12,
    wi: Fq12,
) -> Fq12 {
    // all of  them are fixed-point.
    assert_eq!(eval_points.len(), lines.len());
    // 1. assert c · c^−1 = 1
    assert_eq!(c * c_inv, Fq12::ONE);

    // 2. f = c_inv
    let mut lc = 0_usize;
    let mut f = c_inv;

    let mut naf_digits = biguint_to_naf(e);
    naf_digits.reverse();
    naf_digits.remove(0);

    for (i, digit) in naf_digits.into_iter().enumerate() {
        f = f.square();

        // update c^lambda
        if digit.pow(2) == 1 {
            f = if digit == 1 { f * c_inv } else { f * c };
        }

        for (&P, L) in eval_points.iter().zip(lines) {
            let (alpha, bias) = L[lc];

            let le = MillerLines::line_evaluation(alpha, bias, P);

            f = MillerLines::mul_line_base(f, le.0, le.1, le.2);

            if digit.pow(2) == 1 {
                let (alpha, bias) = L[lc + 1];
                let le = MillerLines::line_evaluation(alpha, bias, P);
                f = MillerLines::mul_line_base(f, le.0, le.1, le.2);
            }
        }

        lc = if digit == 0 { lc + 1 } else { lc + 2 };
    }
    println!("2.f: {:?}", f.to_string());

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

    use ark_bn254::{Fq12, G1Projective, G2Projective};
    use ark_ec::pairing::Pairing;
    use ark_ec::{AffineRepr, CurveGroup, Group};
    use ark_ff::{Field, One};
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;

    use crate::lambda_residues::LambdaResidues;
    use crate::poc::optimal_ate::NativeMillerLoop;
    use crate::{params, poc::constants};

    #[test]
    fn test_dual_miller_loop_with_c_wi_fixed() {
        // 1. setup pairing: (p1, q1)=(p2,q2)
        //      To check (p1, q1)*(p2,-q2)=1
        let p1 = constants::g1
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();
        let p2 = constants::g1
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();

        let q1 = constants::g2
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();
        let q2 = constants::g2
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

        // println!("f1: {:?}", f1.to_string());
        // println!("f2: {:?}", f2.to_string());
        // 2.3 precompute c,wi
        let witness = LambdaResidues::finding_c(f1.mul(f2));
        let c_inv = witness.c.inverse().unwrap();
        // println!("c_inv: {:?}", c_inv.to_string());
        // ====================================
        // ===== 2.Prover compute following data.
        // ====================================
        let verify_res = dual_miller_loop_with_c_wi(
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
}
