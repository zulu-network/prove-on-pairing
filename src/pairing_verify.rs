use ark_bn254::{Fq, Fq12, G1Affine};

// TODO
pub struct Pairing;

impl Pairing {
    pub fn verify_pairings(
        eval_points: Vec<G1Affine>,
        lines: Vec<(Fq, Fq)>,
        e: Fq,
        c: Fq12,
        c_inv: Fq12,
        wi: Fq12,
    ) -> bool {
        true
    }
}
