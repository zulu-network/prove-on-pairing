use ark_bn254::{Fq, Fq12, Fq2, Fq6, Fr, G1Affine};
use ark_ec::AffineRepr;
use ark_ff::Field;

use crate::{
    fields::{fq12_to_frobenius, fq12_to_frobenius_p2, fq12_to_frobenius_p3},
    miller_loop_verify::line_evaluation,
    optimal_ate::mul_line_base,
    utils::scalar_to_naf,
};

pub struct Pairing;

impl Pairing {
    pub fn verify_pairings(
        eval_points: Vec<G1Affine>,
        lines: &[Vec<(Fq2, Fq2)>; 2],
        e: Fr,
        c: Fq12,
        c_inv: Fq12,
        wi: Fq12,
    ) -> Fq12 {
        assert_eq!(eval_points.len(), lines.len());
        assert_eq!(c * c_inv, Fq12::ONE);

        let mut lc = 0_usize;
        // let mut f = Fq12::ONE;
        let mut f = c_inv;
        let mut naf_digits = scalar_to_naf(e);
        naf_digits.reverse();
        naf_digits.remove(0);
        for digit in naf_digits.into_iter() {
            f = f.square();
            // update c^lambda
            if digit.pow(2) == 1 {
                f = if digit == 1 { f * c_inv } else { f * c };
            }
            for (&P, L) in eval_points.iter().zip(lines) {
                let (alpha, bias) = L[lc];
                let le = line_evaluation(alpha, bias, P);
                f = mul_line_base(f, le.0, le.1, le.2);

                if digit.pow(2) == 1 {
                    let (alpha, bias) = L[lc + 1];
                    let le = line_evaluation(alpha, bias, P);
                    f = mul_line_base(f, le.0, le.1, le.2);
                }
            }
            lc = if digit == 0 { lc + 1 } else { lc + 2 };
        }

        // update c^lambda
        f = f * fq12_to_frobenius(c_inv) * fq12_to_frobenius_p2(c) * fq12_to_frobenius_p3(c_inv);
        // update the scalar
        f = f * wi;

        // frobenius map part, p - p^2 + p^3
        for (P, L) in eval_points.into_iter().zip(lines) {
            for k in 0..3 {
                let (alpha, bias) = L[lc + k];
                if k == 2 {
                    let eval = Fq12::new(
                        Fq6::ZERO,
                        Fq6::new(
                            Fq2::ZERO,
                            -bias,
                            Fq2::new(Fq::ZERO, P.x().expect("failed to read x").to_owned()),
                        ),
                    );
                    f = f * eval;
                } else {
                    let le = line_evaluation(alpha, bias, P);
                    f = mul_line_base(f, le.0, le.1, le.2);
                }
            }
        }
        lc = lc + 3;
        assert_eq!(lc, lines[0].len());
        assert_eq!(f, Fq12::ONE);
        f
    }
}
