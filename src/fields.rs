// def conjugate_of(self):
// """
//         For gamma = A + iB \in gfp2
//         gamma^p = A - iB
//         """
// return Fp2(self.x.additive_inverse(), self.y)
// pub fn fq2_conjugate_of()

use ark_bn254::{Fq, Fq12, Fq2, Fq6};
use ark_ff::Field;

// const beta: Fq2 = Fq2::new(Fq::ONE, Fq::from(9));

// TODO: check if origin q12 is changed @Payne
pub fn fq12_to_frobenius(mut q12: Fq12) -> Fq12 {
    // let e1_x = q12.c0.c0.conjugate_in_place() * beta_pi_1[4];
    // let e1_y = q12.c0.c1.conjugate_in_place() * beta_pi_1[2];
    // let e1_z = q12.c0.c2.conjugate_in_place() * beta_pi_1[0];

    // let e2_x = q12.c1.c0.conjugate_in_place() * beta_pi_1[3];
    // let e2_y = q12.c1.c1.conjugate_in_place() * beta_pi_1[1];
    // let e2_z = q12.c1.c2.conjugate_in_place().to_owned();

    // Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e1_x, e1_y, e1_z))
    todo!()
}

pub fn fq12_to_frobenius_p2(q12: Fq12) -> Fq12 {
    // let e1_x = q12.c0.c0 * beta_pi_2[4];
    // let e1_y = q12.c0.c1 * beta_pi_2[2];
    // let e1_z = q12.c0.c2 * beta_pi_2[0];

    // let e2_x = q12.c1.c0 * beta_pi_2[3];
    // let e2_y = q12.c1.c1 * beta_pi_2[1];
    // let e2_z = q12.c1.c2;

    // Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e1_x, e1_y, e1_z))
    todo!()
}

pub fn fq12_to_frobenius_p3(mut q12: Fq12) -> Fq12 {
    // let e1_x = q12.c0.c0.conjugate_in_place() * beta_pi_3[4];
    // let e1_y = q12.c0.c1.conjugate_in_place() * beta_pi_3[2];
    // let e1_z = q12.c0.c2.conjugate_in_place() * beta_pi_3[0];

    // let e2_x = q12.c1.c0.conjugate_in_place() * beta_pi_3[3];
    // let e2_y = q12.c1.c1.conjugate_in_place() * beta_pi_3[1];
    // let e2_z = q12.c1.c2.conjugate_in_place().to_owned();

    // Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e1_x, e1_y, e1_z))
    todo!()
}
