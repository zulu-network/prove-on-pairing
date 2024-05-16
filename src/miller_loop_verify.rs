use ark_bn254::{Fq, Fq12, Fq2, G1Affine};
use ark_ff::Field;
use std::ops::{Mul, Neg};

// we use affine coordinate to verify the line evaluation
// (-b) + y_P * w^3 + (-alpha * x_P) * w^2 where w \in Fp12
fn line_evaluation(alpha: Fq2, bias: Fq2, point: G1Affine) -> (Fq2, Fq2, Fq2) {
    let mut neg_alpha = alpha.neg();
    neg_alpha.mul_assign_by_basefield(&point.x);
    (bias.neg(), neg_alpha, Fq2::new(Fq::ZERO, point.y))
}

#[cfg(test)]
mod test {
    use super::*;
}
