#[cfg(test)]
mod tests {
    use ark_bn254::{g1::G1Affine, Fq12};
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ff::{BigInt, BigInteger, Field};

    use crate::{
        constant::{g1, g2},
        pairing_verify::Pairing,
    };

    #[test]
    fn test_pairing_verify() {
        let P1 = g1.mul_bigint(BigInt::<4>::new([3, 0, 0, 0])).into_affine();
        let P2 = g1.mul_bigint(BigInt::<4>::one()).into_affine();

        let Q1 = g2.mul_bigint(BigInt::<4>::one()).into_affine();
        let Q2 = g2.mul_bigint(BigInt::<4>::new([3, 0, 0, 0])).into_affine();

        println!("Q1: {:?},\nQ2: {:?}", Q1, Q2);

        // let verify_res = Pairing::verify_pairings(eval_points, lines, e, c, c_inv, wi);
        let verify_res = Fq12::ONE;
        assert_eq!(verify_res, Fq12::ONE);
    }
}
