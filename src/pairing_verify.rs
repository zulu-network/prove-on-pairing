use ark_bn254::{Fq, Fq12, Fr, G1Affine};
use ark_ff::Field;

use crate::utils::scalar_to_naf;

// TODO
pub struct Pairing;

impl Pairing {
    pub fn verify_pairings(
        eval_points: Vec<G1Affine>,
        lines: Vec<(Fq, Fq)>,
        e: Fr,
        c: Fq12,
        c_inv: Fq12,
        wi: Fq12,
    ) -> bool {
        assert_eq!(eval_points.len(), lines.len());
        assert_eq!(c * c_inv, Fq12::ONE);

        // let mut f = Fq12::ONE;
        let mut f = c_inv;
        let mut naf_digits = scalar_to_naf(e);
        naf_digits.reverse();
        naf_digits.remove(0);
        for (i, &digit) in naf_digits.iter().enumerate() {}

        true
    }
}
