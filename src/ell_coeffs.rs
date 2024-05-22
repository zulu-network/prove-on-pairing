use ark_bn254::Fq2;
use ark_ec::CurveGroup;
use ark_ff::Field;

// Porting from ark_ec::bn::g2::EllCoeff
pub type EllCoeff = (Fq2, Fq2, Fq2);

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct G2Prepared {
    /// Stores the coefficients of the line evaluations as calculated in
    /// <https://eprint.iacr.org/2013/722.pdf>
    pub ell_coeffs: Vec<EllCoeff>,
}

impl Default for G2Prepared {
    fn default() -> Self {
        Self::from(ark_bn254::G2Affine::generator())
    }
}

impl From<ark_bn254::G2Affine> for G2Prepared {
    // equal with line_function.
    fn from(q: ark_bn254::G2Affine) -> Self {
        assert!(!q.infinity);
        let two_inv = ark_bn254::Fq::one().double().inverse().unwrap();
        let mut ell_coeffs = vec![];
        let mut r = G2HomProjective {
            x: q.x,
            y: q.y,
            z: ark_bn254::Fq2::one(),
        };

        let neg_q = -q;

        for bit in ark_bn254::Config::ATE_LOOP_COUNT.iter().rev().skip(1) {
            ell_coeffs.push(r.double_in_place(&two_inv));

            match bit {
                1 => ell_coeffs.push(r.add_in_place(&q)),
                -1 => ell_coeffs.push(r.add_in_place(&neg_q)),
                _ => continue,
            }
        }

        let q1 = mul_by_char(q);
        let mut q2 = mul_by_char(q1);

        q2.y = -q2.y;

        ell_coeffs.push(r.add_in_place(&q1));
        ell_coeffs.push(r.add_in_place(&q2));

        Self { ell_coeffs }
    }
}

impl From<ark_bn254::G2Projective> for G2Prepared {
    fn from(q: ark_bn254::G2Projective) -> Self {
        q.into_affine().into()
    }
}

impl<'a> From<&'a ark_bn254::G2Affine> for G2Prepared {
    fn from(other: &'a ark_bn254::G2Affine) -> Self {
        (*other).into()
    }
}

impl<'a> From<&'a ark_bn254::G2Projective> for G2Prepared {
    fn from(q: &'a ark_bn254::G2Projective) -> Self {
        q.into_affine().into()
    }
}
