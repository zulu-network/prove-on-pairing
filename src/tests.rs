#[cfg(test)]
mod tests {
    use ark_bn254::{g1::G1Affine, Fq12};
    use ark_ff::Field;

    use crate::pairing_verify::Pairing;

    #[test]
    fn test_pairing_verify() {
        let P1 = G1Affine::new_unchecked(x, y)
        let verify_res = Pairing::verify_pairings(eval_points, lines, e, c, c_inv, wi);
        assert_eq!(verify_res, Fq12::ONE);
    }
}
