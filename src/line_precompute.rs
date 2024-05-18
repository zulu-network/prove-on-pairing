use crate::constant::{E, LAMBDA, MODULUS};
use crate::fields::Fq12Ext;
use crate::utils::{biguint_to_naf, to_naf};
use ark_bn254::{Fq, Fq12, Fq2, G1Affine, G2Affine, G2Projective};
use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_ff::{Field, PrimeField};
use num_bigint::BigUint;
use num_traits::Zero;
use std::ops::{Add, Div, Mul, Neg};

// stands for (alpha, bias)
type LiearRes = (Fq2, Fq2);

fn line_double(point: &G2Affine) -> LiearRes {
    // T = T.force_affine()
    // assert(T.z == T.one_element())
    // x, y = T.x, T.y
    let (x, y) = (point.x, point.y);

    // slope: alpha = 3 * x ^ 2 / (2 * y)
    let alpha = x.square().mul(Fq2::from(3)).div(y.mul(Fq2::from(2)));
    // bias = y - alpha * x
    let bias = y - alpha * x;

    (alpha, bias)
}

fn line_add(point: &G2Affine, other: &G2Affine) -> LiearRes {
    let (x1, y1) = (point.x, point.y);
    let (x2, y2) = (other.x, other.y);

    // slope: alpha = (y2 - y1) / (x2 - x1)
    let alpha = (y2 - y1) / (x2 - x1);
    // bias = y1 - alpha * x1
    let bias = y1 - alpha * x1;
    (alpha, bias)
}

// cache line line_function for [6x + 2 + p - p^2 + p^3]Q
// Input:
//      Q: The Fixed Point
fn line_function(Q: G2Projective, e: BigUint, lamb: BigUint) -> Vec<LiearRes> {
    let mut point_naf = biguint_to_naf(e.clone());
    point_naf.reserve(0);
    let naf_digits = point_naf[1..].to_vec();
    let mut line_vec = vec![];

    let mut T = Q.clone();

    // 1. double-add part, 6x + 2
    naf_digits.into_iter().enumerate().map(|(i, digit)| {
        let double_res = line_double(&T.into_affine());
        T = T.double();
        line_vec.push(double_res);
        if digit ^ 2 == 1 {
            let qt = if 1 == digit {
                Q.clone()
            } else {
                Q.clone().neg()
            };

            let qt_double_res = line_add(&T.into_affine(), &qt.into_affine());
            line_vec.push(qt_double_res);
        }
    });
    // assert_eq!(T, Q.mul_bigint(e));

    // 2. frobenius map part, p - p^2 + p^3
    // 2.1 Q1 = pi(Q)
    // x = x' * beta^(2 * (p - 1) / 6)
    // y = y' * beta^(3 * (p - 1) / 6))
    let (mut x, mut y) = (Q.x, Q.y);

    let pi_1_Q = G2Projective::new(
        x.conjugate_in_place().mul(Fq12Ext::beta_pi_1()[1]),
        y.conjugate_in_place().mul(Fq12Ext::beta_pi_1()[2]),
        Fq2::ONE,
    );
    // assert!(pi_1_Q.into_affine().is_on_curve());
    // assert_eq!(pi_1_Q, Q.into_affine().mul_bigint(e.clone()));

    // 2.2. Q2 = pi2(Q)
    // x = x * beta * (2 * (p^2 - 1) / 6)
    // y = y * beta * (3 * (p^2 - 1) / 6) = -y
    let (mut x, mut y) = (Q.x, Q.y);
    let pi_2_Q = G2Projective::new(
        x.conjugate_in_place().mul(Fq12Ext::beta_pi_2()[1]),
        y.conjugate_in_place().mul(Fq12Ext::beta_pi_2()[2]),
        Fq2::ONE,
    );
    // assert!(pi_2_Q.into_affine().is_on_curve());
    // assert_eq!(
    //     pi_2_Q,
    //     Q.mul_bigint(&BigUint::from(Fq::MODULUS).pow(2).to_u64_digits())
    // );

    // 2.3. Q3 = pi3(Q)
    // x = x' * beta * (2 * (p^3 - 1) / 6)
    // y = y' * beta * (3 * (p^3 - 1) / 6)
    let (mut x, mut y) = (Q.x, Q.y);
    let pi_3_Q = G2Projective::new(
        x.conjugate_in_place().mul(Fq12Ext::beta_pi_3()[1]),
        y.conjugate_in_place().mul(Fq12Ext::beta_pi_3()[2]),
        Fq2::ONE,
    );
    // assert!(pi_3_Q.into_affine().is_on_curve());
    // assert_eq!(pi_3_Q, Q.mul_bigint(BigUint::from(Fq::MODULUS).pow(3)));

    let line_pi_1 = line_add(&T.into_affine(), &pi_1_Q.into_affine());
    T = T.add(pi_1_Q);
    line_vec.push(line_pi_1);
    // assert_eq!(T, Q.mul_bigint(BigUint::from(Fq::MODULUS).add(E)));

    let line_pi_2 = line_add(&T.into_affine(), &pi_2_Q.into_affine());
    T = T.add(pi_2_Q.neg());
    line_vec.push(line_pi_2);

    // k = p - p^2 + e
    let k = BigUint::from(Fq::MODULUS) - BigUint::from(Fq::MODULUS).pow(2);
    let k = k + e;
    // assert_eq!(T, Q.mul_bigint(if k.gt(BigUint::ZERO){
    //     k
    // }else {
    // TODO
    //     // rx(x) - (-k % rx(x))
    // } ));

    let line_i = line_add(&T.into_affine(), &pi_3_Q.into_affine());
    let T = T.add(pi_3_Q);
    line_vec.push(line_i);

    // assert
    // assert_eq!(T, Q.mul_bigint(LAMBDA));
    // assert!(T.into_affine().is_on_curve());
    // assert!(line_i.0.is_zero());

    line_vec
}

#[cfg(test)]
mod test {
    use ark_bn254::{Fq, Fq12, Fq2, Fq6, Fr, G1Affine, G1Projective};
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ff::Field;
    use ark_std::{end_timer, start_timer};
    use num_bigint::BigUint;
    use num_traits::{FromPrimitive, One};
    use std::str::FromStr;

    use super::*;
    use crate::constant::{g1, g2, E};
    use crate::{
        fields::{fq12_to_frobenius, fq12_to_frobenius_p2, fq12_to_frobenius_p3},
        miller_loop_verify::line_evaluation,
        optimal_ate::mul_line_base,
        utils::to_naf,
    };

    #[test]
    fn test_line_precomputation() {
        // assume we want to prove e(P1, Q1) = e(P2, Q2), namely e(P1, Q1) * e(P2, -Q2) = 1
        // fixed point Q in G2, public known to verifier
        let P1 = g1
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();
        let P2 = g1.mul_bigint(BigUint::one().to_u64_digits()).into_affine();

        let (P1, P2) = (G1Projective::from(P1), G1Projective::from(P2));

        let Q1 = g2.mul_bigint(BigUint::one().to_u64_digits()).into_affine();
        let Q2 = g2
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();
        let (Q1, Q2) = (G2Projective::from(Q1), G2Projective::from(Q2));

        let lamb = BigUint::from_str(
            "10486551571378427818905133077457505975146652579011797175399169355881771981095211883813744499745558409789005132135496770941292989421431235276221147148858384772096778432243207188878598198850276842458913349817007302752534892127325269"
        ).unwrap();

        let start = start_timer!(|| "start compute line precomputation");
        // indexer (oracle) for fixed point Q
        let L1 = line_function(Q1, E.clone(), lamb.clone());
        let L2 = line_function(Q2.neg(), E.clone(), lamb);
        end_timer!(start);
    }
}
