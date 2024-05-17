// def conjugate_of(self):
// """
//         For gamma = A + iB \in gfp2
//         gamma^p = A - iB
//         """
// return Fp2(self.x.additive_inverse(), self.y)
// pub fn fq2_conjugate_of()

use ark_bn254::{Fq, Fq12, Fq2};
use ark_ff::Field;

// const beta: Fq2 = Fq2::new(Fq::ONE, Fq::from(9));

pub fn fq12_to_frobenius(field: Fq12) -> Fq12 {
    todo!()
}

pub fn fq12_to_frobenius_p2(field: Fq12) -> Fq12 {
    todo!()
}

pub fn fq12_to_frobenius_p3(field: Fq12) -> Fq12 {
    todo!()
}
