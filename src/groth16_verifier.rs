use crate::lambda_residues::LambdaResidues;
use crate::pairing_verifier::quad_miller_loop_with_c_wi;
use crate::params;
use ark_bn254::{Bn254, Fr, G1Affine, G1Projective};
use ark_ec::bn::{BnConfig, G2Prepared};
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::Field;
use ark_groth16::{PreparedVerifyingKey, Proof};
use ark_relations::r1cs::Result as R1CSResult;
use std::ops::Neg;

// Only support Bn254 for now.
pub struct Groth16Verifier;

impl Groth16Verifier {
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
        let g_ic = pvk.vk.gamma_abc_g1[0].into_group();

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
        assert_eq!(hint, c.pow(p_pow3.to_u64_digits()), "hint is wrong");
        assert_eq!(final_f, hint, "final_f not equal hint");

        Ok(true)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::dummy_circuit::gen_groth16_dummy_circuit_proof;
    use ark_bn254::Bn254;
    use ark_groth16::Groth16;

    #[test]
    fn test_groth16_verifier() {
        type E = Bn254;

        let k = 6;

        // 1. gen proof
        let (proof, pvk, pi) = gen_groth16_dummy_circuit_proof::<E>(k);

        // 2. verify with native verifier
        assert!(
            Groth16::<E>::verify_proof(&pvk, &proof, &pi).unwrap(),
            "native verifier can't pass"
        );

        // 3. verifier with new one
        assert!(Groth16Verifier::verify_proof_with_c_wi(&pvk, &proof, &pi).unwrap());
    }
}
