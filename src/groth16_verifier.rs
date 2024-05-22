// use ark_ec::pairing::Pairing;
// use ark_ec::CurveGroup;
// use ark_groth16::r1cs_to_qap::{LibsnarkReduction, R1CSToQAP};
// use ark_groth16::{Groth16, PreparedVerifyingKey, Proof};
// use ark_relations::r1cs::{Result as R1CSResult, SynthesisError};
// use std::marker::PhantomData;
//
// pub struct Groth16Verifier<E: Pairing>(PhantomData<E>);
//
// impl<E: Pairing> Groth16Verifier<E> {
//     // Porting from ark_groth16::Groth16::verify_proof
//     pub fn verify_proof(
//         pvk: &PreparedVerifyingKey<E>,
//         proof: &Proof<E>,
//         public_inputs: &[E::ScalarField],
//     ) -> R1CSResult<bool> {
//         let prepared_inputs = Groth16::prepare
//         Self::verify_proof_with_prepared_inputs(pvk, proof, &prepared_inputs)
//     }
//
//     /// Verify a Groth16 proof `proof` against the prepared verification key `pvk` and prepared public
//     /// inputs. This should be preferred over [`verify_proof`] if the instance's public inputs are
//     /// known in advance.
//     ///
//     /// Verifier by applying with new paper:
//     pub fn verify_proof_with_prepared_inputs(
//         pvk: &PreparedVerifyingKey<E>,
//         proof: &Proof<E>,
//         prepared_inputs: &E::G1,
//     ) -> R1CSResult<bool> {
//         let qap = E::multi_miller_loop(
//             [
//                 <E::G1Affine as Into<E::G1Prepared>>::into(proof.a),
//                 prepared_inputs.into_affine().into(),
//                 proof.c.into(),
//             ],
//             [
//                 proof.b.into(),
//                 pvk.gamma_g2_neg_pc.clone(),
//                 pvk.delta_g2_neg_pc.clone(),
//             ],
//         );
//
//         let test = E::final_exponentiation(qap).ok_or(SynthesisError::UnexpectedIdentity)?;
//
//         Ok(test.0 == pvk.alpha_g1_beta_g2)
//     }
// }
