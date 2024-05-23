use ark_bn254::{fq::Fq, Bn254, Fq12, Fq2, Fq6, G1Affine, G2Affine, G2Projective};
use ark_ec::bn::g2::{mul_by_char, G2HomProjective};
use ark_ec::bn::{BnConfig, G2Prepared};
use ark_ec::AffineRepr;
use ark_ff::Field;
use ark_std::UniformRand;
use num_bigint::BigUint;
use num_traits::One;
use rand::SeedableRng;
use std::ops::{Mul, Neg};
use std::str::FromStr;

use crate::lambda_residues::LambdaResidues;
use crate::params::MODULUS;
use crate::utils::biguint_to_naf;
use crate::{
    miller_lines::MillerLines,
    params,
    utils::{fq12_to_frobenius, fq12_to_frobenius_p2, fq12_to_frobenius_p3},
};

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

// To verify e(P1,Q1)*e(P2,Q2)*e(P3,Q3)*e(P4,Q4)=1
//
// Here is only support to verify groth16's pairing, which (Q1,Q2,Q3) are fixed, Q4 is non-fixed.
//
// params:
//  @eval_points: [P1,P2,P3]. which has fixed {Q1,Q2,Q3}
//  @P4: P4
//  @Q4: Q4
//  @lines: []precompute miller lines for Qi. Only support fixed Qi.
//  @c: c^lambda = f*w^i
//  @c_inv: inverse of c
//
// verify c^lambda = f * wi, namely c_inv^lambda * f * wi = 1
pub fn quad_miller_loop_with_c_wi(
    eval_points: Vec<G1Affine>,
    P4: G1Affine,
    Q4: G2Affine,
    // lines: &[Vec<(Fq2, Fq2)>],
    lines: &Vec<G2Prepared<ark_bn254::Config>>,
    c: Fq12,
    c_inv: Fq12,
    wi: Fq12,
    // TODO: What's B in stack
) -> Fq12 {
    assert_eq!(eval_points.len(), 3, "Should contains 4 G1Affine: P1,P2,P3");
    assert_eq!(lines.len(), 3, "Only precompute lines for Q1,Q2,Q3");
    assert_eq!(c * c_inv, Fq12::ONE, "Check if c·c^−1 = 1");

    // let P4 = eval_points[3].clone();
    // let Q4_projective: G2Projective = Q4.into_group();
    let mut T4 = G2HomProjective::<ark_bn254::Config> {
        x: Q4.x,
        y: Q4.y,
        z: ark_bn254::Fq2::one(),
    };

    // constants
    let two_inv = ark_bn254::Fq::one().double().inverse().unwrap();

    // 1. f = c_inv
    let mut f = c_inv;
    println!("1.f: {:?}", f.to_string());

    let mut lines_iters = lines
        .iter()
        .map(|item| item.ell_coeffs.iter())
        .collect::<Vec<_>>();

    // 2. miller loop part, 6x + 2
    for i in (1..ark_bn254::Config::ATE_LOOP_COUNT.len()).rev() {
        let bit = ark_bn254::Config::ATE_LOOP_COUNT[i - 1];

        // 2.1 double: f = f * f
        f = f.square();

        // 2.2 mul c
        //  f = f * c_inv, if digit == 1
        //  f = f * c, if digit == -1
        f = if 1 == bit {
            f * c_inv
        } else if bit == -1 {
            f * c
        } else if bit == 0 {
            f
        } else {
            panic!("bit is not in (-1,1), bit={bit}");
        };

        // 2.3 accumulate double lines (fixed and non-fixed)
        // 2.3.1(fixed) f = f * double_line_Q(P). fixed points: P1, P2, P3
        for (line_i, pi) in lines_iters.iter_mut().zip(eval_points.iter()) {
            let line_i_0 = line_i.next().unwrap();
            Bn254::ell(&mut f, line_i_0, pi);
        }

        // 2.3.2(non-fixed) double line with T4 (projective coordinates)
        let double_line = T4.double_in_place(&two_inv); // TODO: check if the param is 1/2

        // 2.3.3(non-fixed) evaluation double_line. non-fixed points: P4
        Bn254::ell(&mut f, &double_line, &P4);

        if bit == 1 || bit == -1 {
            // 2.4 accumulate add lines (fixed and non-fixed)
            // 2.4.1(fixed) f = f * add_line_eval. fixed points: P1, P2, P3
            for (line_i, pi) in lines_iters.iter_mut().zip(eval_points.iter()) {
                let line_i_1 = line_i.next().unwrap();
                Bn254::ell(&mut f, line_i_1, pi);
            }
            // 2.4.2(non-fixed) double line with T4 (projective coordinates)
            let add_line = T4.add_in_place(&Q4);

            // 2.4.3(non-fixed) evaluation double_line. non-fixed points: P4
            Bn254::ell(&mut f, &add_line, &P4);
        }
    }
    println!("2.f: {:?}", f.to_string());

    // 3. f = f * c_inv^p * c^{p^2}
    f = f
        * c_inv.pow(params::MODULUS.to_u64_digits())
        * c.pow(params::MODULUS.pow(2).to_u64_digits());
    println!("3.f: {:?}", f.to_string());

    // 4. f = f * wi . scale f
    f = f * wi;
    println!("4.f: {:?}", f.to_string());

    // 5 add lines (fixed and non-fixed)
    // 5.1(fixed) f = f * add_line_eval. fixed points: P1, P2, P3
    for (line_i, pi) in lines_iters.iter_mut().zip(eval_points.iter()) {
        let line_i_1 = line_i.next().unwrap();
        Bn254::ell(&mut f, line_i_1, pi);
    }
    // 5.2 one-time frobenius map to compute phi_Q
    //     compute phi(Q) with Q4
    let phi_Q = mul_by_char::<ark_bn254::Config>(Q4.clone());

    let add_line = T4.add_in_place(&phi_Q);

    // 5.4(non-fixed) evaluation add_lin. non-fixed points: P4
    Bn254::ell(&mut f, &add_line, &P4);
    println!("5.f: {:?}", f.to_string());

    // 6. add lines (fixed and non-fixed)
    // 6.1(fixed) f = f * add_line_eval. fixed points: P1, P2, P3
    for (line_i, pi) in lines_iters.iter_mut().zip(eval_points.iter()) {
        // TODO: where is f?? and where is double line?
        let line_i_1 = line_i.next().unwrap();
        Bn254::ell(&mut f, line_i_1, pi);
    }
    // 6.2 two-time frobenius map to compute phi_Q
    //     compute phi_Q_2 with phi_Q
    // mul_by_char: used to q's frob...map.
    let phi_Q_2 = mul_by_char::<ark_bn254::Config>(phi_Q.clone());

    let add_line = T4.add_in_place(&phi_Q_2);
    println!("6.2.f: {:?}", f.to_string());

    // 6.3(non-fixed) evaluation add_lin. non-fixed points: P4
    Bn254::ell(&mut f, &add_line, &P4);
    println!("6.3.f: {:?}", f.to_string());

    // return final_f
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
    fn test_dual_miller_loop_with_c_wi_fixed() {
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

    #[test]
    #[ignore]
    fn test_dual_miller_loop_with_c_wi_fixed_with_4_point() {
        // 1. setup pairing: (p1, q1)(p3, q3)=(p2,q2)(p4,q4)
        //      To check (p1, q1)*(p2,-q2)*(p3,q3)*(p4,-q4)=1
        let p1 = dev::g1
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();
        let p2 = dev::g1
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();

        let p3 = dev::g1
            .mul_bigint(BigUint::from_i8(5).unwrap().to_u64_digits())
            .into_affine();
        let p4 = dev::g1
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();

        let q1 = dev::g2
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();
        let q2 = dev::g2
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();

        let q3 = dev::g2
            .mul_bigint(BigUint::one().to_u64_digits())
            .into_affine();
        let q4 = dev::g2
            .mul_bigint(BigUint::from_i8(5).unwrap().to_u64_digits())
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
        let l3 = MillerLines::precompute_lines(
            G2Projective::from(q3),
            params::E.clone(),
            params::LAMBDA.clone(),
        );
        let l4 = MillerLines::precompute_lines(
            G2Projective::from(q4.neg()),
            params::E.clone(),
            params::LAMBDA.clone(),
        );

        // 2.2 precompute witness
        // TODO: meet error when replacing it with compute.
        let f1 = NativeMillerLoop::miller_loop(G1Projective::from(p1), G2Projective::from(q1));
        let f2 =
            NativeMillerLoop::miller_loop(G1Projective::from(p2), G2Projective::from(q2).neg());
        let f3 = NativeMillerLoop::miller_loop(G1Projective::from(p3), G2Projective::from(q3));
        let f4 =
            NativeMillerLoop::miller_loop(G1Projective::from(p4), G2Projective::from(q4).neg());

        println!("f1: {:?}", f1.to_string());
        println!("f2: {:?}", f2.to_string());
        // 2.3 precompute c,wi
        let witness = LambdaResidues::finding_c(f1.mul(f2).mul(f3).mul(f4));
        let c_inv = witness.c.inverse().unwrap();
        println!("c_inv: {:?}", c_inv.to_string());
        // ====================================
        // ===== 2.Prover compute following data.
        // ====================================
        let verify_res = dual_miller_loop_with_c_wi(
            vec![p1, p2, p3, p4],
            &[l1, l2, l3, l4],
            params::E.clone(),
            witness.c,
            c_inv,
            witness.wi,
        );
        assert_eq!(verify_res, Fq12::ONE);
        println!("========Successfully");
    }
}
