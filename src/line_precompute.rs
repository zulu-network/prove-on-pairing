use crate::utils::to_naf;
use ark_bn254::{Fq, Fq12, Fq2, G1Affine, G2Affine};
use ark_ff::Field;
use std::ops::{Div, Mul, Neg};

// stands for (alpha, bias)
type LiearRes = (Fq2, Fq2);

fn line_double(point: G2Affine) -> LiearRes {
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

fn line_add(point: G2Affine, other: G2Affine) -> LiearRes {
    let (x1, y1) = (point.x, point.y);
    let (x2, y2) = (other.x, other.y);

    // slope: alpha = (y2 - y1) / (x2 - x1)
    let alpha = (y2 - y1) / (x2 - x1);
    // bias = y1 - alpha * x1
    let bias = y1 - alpha * x1;
    (alpha, bias)
}

// TODO: need redefine the type of params(e)
// cache line line_function for [6x + 2 + p - p^2 + p^3]Q
fn cache_line_function(Q: G2Affine, e: i128, lamb: i128) {
    let mut point_naf = to_naf(e);
    point_naf.reserve(0);
    let naf_digits = point_naf[1..].to_vec();

    // let mut line_vec = vec![];

    let mut T: G2Affine = Q.clone();
    // 1. double-add part, 6x + 2

    // naf_digits.into_iter().enumerate().map(|(i, digit)| {
    //     let double_res = line_double(T.clone());
    //     // T = G2Affine::from(T.mul(T.clone())); // TODO
    //     line_vec.push(double_res);
    //     if digit ^ 2 == 1 {
    //         let qt: G2Affine = if 1==digit {
    //             Q.clone()
    //         } else {
    //             Q.clone().neg()
    //         };
    //
    //         let qt_double_res = line_add(T.clone(), qt);
    //         line_vec.push(qt_double_res);
    //     }
    // });
    // assert_eq!(T, );

    // 2. frobenius map part, p - p^2 + p^3
    // Q1 = pi(Q)
    // x = x' * beta^(2 * (p - 1) // 6)
    // y = y' * beta^(3 * (p - 1) // 6))
    let (mut x, mut y) = (Q.x, Q.y);
    let x2 = *x.conjugate_in_place();
    let x2 = Fq2::frobenius_map(&x2, 1);
    let y2 = *y.conjugate_in_place();
    let y2 = Fq2::frobenius_map(&y2, 1);
    let pi_1_Q = G2Affine::new(x2, y2);
    assert!(pi_1_Q.is_on_curve());
    // assert_eq!(pi_1_Q, Q.mul_bigint(e));

    let pi_2_Q = G2Affine::new(Fq2::frobenius_map(&Q.x, 1), Fq2::frobenius_map(&Q.y, 1));
    assert!(pi_2_Q.is_on_curve());
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_line_precomputation() {
        // let
    }
}
