// Rephrased from https://github.com/arkworks-rs/algebra/blob/master/ec/src/models/bn/g2.rs#L185
// Cannot directly obtain G2 because of visibility
//
// Note: BitVM::bn254::ell_coeffs also rephrased following structs

use crate::ell_coeffs::EllCoeff;
use ark_bn254::{Bn254, Fq12, Fr, G1Affine, G1Projective};
use ark_ec::bn::{BnConfig, TwistType};
use ark_ec::pairing::{MillerLoopOutput, Pairing};
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{CyclotomicMultSubgroup, Field, PrimeField};
use ark_groth16::{PreparedVerifyingKey, Proof};
use ark_relations::r1cs::{Result as R1CSResult, SynthesisError};
use ark_std::cfg_chunks_mut;
use itertools::Itertools;
use num_traits::One;
use std::ops::AddAssign;

// Only support Bn254 for now.
pub struct Groth16Verifier;

impl Groth16Verifier {
    // Porting from ark_groth16::Groth16::verify_proof
    pub fn verify_proof(
        pvk: &PreparedVerifyingKey<Bn254>,
        proof: &Proof<Bn254>,
        public_inputs: &[Fr],
    ) -> R1CSResult<bool> {
        let prepared_inputs = Self::prepare_inputs(pvk, public_inputs)?;
        Self::verify_proof_with_prepared_inputs(pvk, proof, &prepared_inputs)
    }

    // Porting from ark_groth16::Groth16::prepare_inputs
    pub fn prepare_inputs(
        pvk: &PreparedVerifyingKey<Bn254>,
        public_inputs: &[Fr],
    ) -> R1CSResult<G1Projective> {
        assert_eq!(public_inputs.len() + 1, pvk.vk.gamma_abc_g1.len());
        // g_ic = pvk.vk.gamma_abc_g1[0] + msm(pvk.vk.gamma_abc_g1[1..], public_inputs)
        let mut g_ic = pvk.vk.gamma_abc_g1[0].into_group();

        // msm(pvk.vk.gamma_abc_g1[1..], public_inputs)
        for (i, b) in public_inputs.iter().zip(pvk.vk.gamma_abc_g1.iter().skip(1)) {
            g_ic.add_assign(&b.mul_bigint(i.into_bigint()));
        }

        Ok(g_ic)
    }

    // Verifier by applying with new paper:
    //
    // Porting from groth16::verify_proof_with_prepared_inputs
    pub fn verify_proof_with_prepared_inputs(
        pvk: &PreparedVerifyingKey<Bn254>,
        proof: &Proof<Bn254>,
        prepared_inputs: &G1Projective,
    ) -> R1CSResult<bool> {
        let p_a: Vec<<Bn254 as Pairing>::G1Prepared> = vec![
            <G1Affine as Into<<Bn254 as Pairing>::G1Prepared>>::into(proof.a),
            prepared_inputs.into_affine().into(),
            proof.c.into(),
        ];

        let q_b: Vec<<Bn254 as Pairing>::G2Prepared> = vec![
            proof.b.into(),
            pvk.gamma_g2_neg_pc.clone(),
            pvk.delta_g2_neg_pc.clone(),
        ];

        let qap = Self::multi_miller_loop(p_a, q_b);

        let test = Bn254::final_exponentiation(qap).ok_or(SynthesisError::UnexpectedIdentity)?;

        Ok(test.0 == pvk.alpha_g1_beta_g2)
    }

    // porting from bn254::multi_miller_loop
    fn multi_miller_loop(
        a: Vec<<Bn254 as Pairing>::G1Prepared>,
        b: Vec<<Bn254 as Pairing>::G2Prepared>,
    ) -> MillerLoopOutput<Bn254> {
        let mut pairs = a
            .into_iter()
            .zip_eq(b)
            .filter_map(|(p, q)| {
                // let (p, q) = (p.into(), q.into());
                match !p.is_zero() && !q.is_zero() {
                    true => Some((p, q.ell_coeffs.into_iter())),
                    false => None,
                }
            })
            .collect::<Vec<_>>();

        let mut f = cfg_chunks_mut!(pairs, 4)
            .map(|pairs| {
                let mut f = Fq12::one();
                for i in (1..ark_bn254::Config::ATE_LOOP_COUNT.len()).rev() {
                    if i != ark_bn254::Config::ATE_LOOP_COUNT.len() - 1 {
                        f.square_in_place();
                    }

                    for (p, coeffs) in pairs.iter_mut() {
                        Self::ell(&mut f, &coeffs.next().unwrap(), &p.0);
                    }

                    let bit = ark_bn254::Config::ATE_LOOP_COUNT[i - 1];
                    if bit == 1 || bit == -1 {
                        for (p, coeffs) in pairs.iter_mut() {
                            Self::ell(&mut f, &coeffs.next().unwrap(), &p.0);
                        }
                    }
                }
                f
            })
            .product::<Fq12>();

        if ark_bn254::Config::X_IS_NEGATIVE {
            f.cyclotomic_inverse_in_place();
        }

        for (p, coeffs) in &mut pairs {
            Self::ell(&mut f, &coeffs.next().unwrap(), &p.0);
        }

        for (p, coeffs) in &mut pairs {
            Self::ell(&mut f, &coeffs.next().unwrap(), &p.0);
        }

        MillerLoopOutput(f)
    }

    // It's Fq12's ell
    //
    // Porting from ark_bn254::Bn254::ell
    fn ell(f: &mut Fq12, coeffs: &EllCoeff, p: &G1Affine) {
        let mut c0 = coeffs.0;
        let mut c1 = coeffs.1;
        let mut c2 = coeffs.2;

        match ark_bn254::Config::TWIST_TYPE {
            TwistType::M => {
                c2.mul_assign_by_fp(&p.y);
                c1.mul_assign_by_fp(&p.x);
                f.mul_by_014(&c0, &c1, &c2);
            }
            TwistType::D => {
                c0.mul_assign_by_fp(&p.y);
                c1.mul_assign_by_fp(&p.x);
                f.mul_by_034(&c0, &c1, &c2);
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::groth16::gen_dummy_groth16_proof;
    use ark_bn254::Bn254;
    use ark_ec::pairing::Pairing;
    use ark_groth16::{prepare_verifying_key, Groth16};

    #[test]
    fn test_groth16_verifier() {
        type E = Bn254;

        let (proof, pvk, pi) = gen_dummy_groth16_proof::<E>();
        assert!(Groth16Verifier::verify_proof(&pvk, &proof, &pi).unwrap());
    }
}
