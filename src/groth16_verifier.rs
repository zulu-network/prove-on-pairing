use crate::lambda_residues::LambdaResidues;
use crate::pairing_verifier::quad_miller_loop_with_c_wi;
use crate::params;
use ark_bn254::{Bn254, Fq12, Fr, G1Affine, G1Projective};
use ark_ec::bn::{BnConfig, G2Prepared};
use ark_ec::pairing::{MillerLoopOutput, Pairing};
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::{CyclotomicMultSubgroup, Field, PrimeField};
use ark_groth16::{Groth16, PreparedVerifyingKey, Proof};
use ark_relations::r1cs::{Result as R1CSResult, SynthesisError};
use ark_std::cfg_chunks_mut;
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::One;
use std::ops::{AddAssign, Neg};

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
        Groth16::<Bn254>::verify_proof_with_prepared_inputs(pvk, proof, &prepared_inputs)
    }
    pub fn verify_proof_with_c_wi(
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
        // for (i, b) in public_inputs.iter().zip(pvk.vk.gamma_abc_g1.iter().skip(1)) {
        //     g_ic.add_assign(&b.mul_bigint(i.into_bigint()));
        // }
        let g_ic = g_ic + G1Projective::msm(&pvk.vk.gamma_abc_g1[1..], &public_inputs).unwrap();

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
        let p_pow3 = params::MODULUS.pow(3_u32);
        let lambda = params::LAMBDA.clone();
        let (exp, sign) = if lambda > p_pow3 {
            (lambda - p_pow3, true)
        } else {
            (p_pow3 - lambda, false)
        };

        let beta_prepared: G2Prepared<ark_bn254::Config> = (pvk.vk.beta_g2.clone().neg()).into();
        let gamma_g2_neg_pc = pvk.gamma_g2_neg_pc.clone();
        let delta_g2_neg_pc = pvk.delta_g2_neg_pc.clone();
        let q_prepared_lines = [gamma_g2_neg_pc, delta_g2_neg_pc, beta_prepared.clone()].to_vec();
        let sum_ai_abc_gamma = prepared_inputs.into_affine();

        let a = vec![
            sum_ai_abc_gamma.into(),
            proof.c.into(),
            pvk.vk.alpha_g1.into(),
            <G1Affine as Into<<Bn254 as Pairing>::G1Prepared>>::into(proof.a),
        ];
        let b = vec![
            pvk.gamma_g2_neg_pc.clone(),
            pvk.delta_g2_neg_pc.clone(),
            beta_prepared,
            proof.b.into(),
        ];

        // compute f. base line
        let qap = Bn254::multi_miller_loop(a.clone(), b.clone());
        let f = qap.0;
        // finding_c
        let witness = LambdaResidues::finding_c(f);
        let (c, wi) = (witness.c, witness.wi);
        let c_inv = c.inverse().unwrap();

        let eval_points = vec![sum_ai_abc_gamma, proof.c, pvk.vk.alpha_g1];

        let final_f = quad_miller_loop_with_c_wi(
            eval_points,
            proof.a,
            proof.b,
            &q_prepared_lines,
            c,
            c_inv,
            wi,
        );

        // check
        let hint = if sign {
            f * wi * (c_inv.pow(exp.to_u64_digits()))
        } else {
            f * wi * (c_inv.pow(exp.to_u64_digits()).inverse().unwrap())
        };

        let p_pow3 = params::MODULUS.pow(3_u32);
        assert_eq!(hint, c.pow(p_pow3.to_u64_digits()));
        assert_eq!(final_f, hint);

        Ok(true)
    }

    // porting from Bn254::multi_miller_loop
    pub fn multi_miller_loop(
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
                        Bn254::ell(&mut f, &coeffs.next().unwrap(), &p.0);
                    }

                    let bit = ark_bn254::Config::ATE_LOOP_COUNT[i - 1];
                    if bit == 1 || bit == -1 {
                        for (p, coeffs) in pairs.iter_mut() {
                            Bn254::ell(&mut f, &coeffs.next().unwrap(), &p.0);
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
            Bn254::ell(&mut f, &coeffs.next().unwrap(), &p.0);
        }

        for (p, coeffs) in &mut pairs {
            Bn254::ell(&mut f, &coeffs.next().unwrap(), &p.0);
        }

        MillerLoopOutput(f)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::circuit::gen_dummy_groth16_proof;
    use ark_bn254::Bn254;
    use ark_ec::pairing::Pairing;

    #[test]
    fn test_groth16_verifier_native() {
        type E = Bn254;

        let (proof, pvk, pi) = gen_dummy_groth16_proof::<E>();
        assert!(Groth16Verifier::verify_proof(&pvk, &proof, &pi).unwrap());
    }

    #[test]
    fn test_groth16_verifier_with_c_wi() {
        type E = Bn254;

        let (proof, pvk, pi) = gen_dummy_groth16_proof::<E>();
        assert!(Groth16Verifier::verify_proof_with_c_wi(&pvk, &proof, &pi).unwrap());
    }
}
