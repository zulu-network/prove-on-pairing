// def conjugate_of(self):
// """
//         For gamma = A + iB \in gfp2
//         gamma^p = A - iB
//         """
// return Fp2(self.x.additive_inverse(), self.y)
// pub fn fq2_conjugate_of()

use crate::constant;
use ark_bn254::{Fq12, Fq6};
use ark_ff::Field;
use num_traits::{FromPrimitive, One, Pow};
use std::ops::{DivAssign, Mul, MulAssign, SubAssign};

pub fn fq12_to_frobenius(mut q12: Fq12) -> Fq12 {
    let e2_z = q12.c1.c2.conjugate_in_place().mul(constant::BETA_PI_1[4]);
    let e2_y = q12.c1.c1.conjugate_in_place().mul(constant::BETA_PI_1[2]);
    let e2_x = q12.c1.c0.conjugate_in_place().mul(constant::BETA_PI_1[0]);

    let e1_z = q12.c0.c2.conjugate_in_place().mul(constant::BETA_PI_1[3]);
    let e1_y = q12.c0.c1.conjugate_in_place().mul(constant::BETA_PI_1[1]);
    let e1_x = q12.c0.c0.conjugate_in_place().to_owned();

    Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e2_x, e2_y, e2_z))
}

pub fn fq12_to_frobenius_p2(q12: Fq12) -> Fq12 {
    let e2_z = q12.c1.c2 * constant::BETA_PI_2[4];
    let e2_y = q12.c1.c1 * constant::BETA_PI_2[2];
    let e2_x = q12.c1.c0 * constant::BETA_PI_2[0];

    let e1_z = q12.c0.c2 * constant::BETA_PI_2[3];
    let e1_y = q12.c0.c1 * constant::BETA_PI_2[1];
    let e1_x = q12.c0.c0;

    Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e2_x, e2_y, e2_z))
}

pub fn fq12_to_frobenius_p3(mut q12: Fq12) -> Fq12 {
    let e2_z = q12.c1.c2.conjugate_in_place().mul(constant::BETA_PI_3[4]);
    let e2_y = q12.c1.c1.conjugate_in_place().mul(constant::BETA_PI_3[2]);
    let e2_x = q12.c1.c0.conjugate_in_place().mul(constant::BETA_PI_3[0]);
    let e1_z = q12.c0.c2.conjugate_in_place().mul(constant::BETA_PI_3[3]);
    let e1_y = q12.c0.c1.conjugate_in_place().mul(constant::BETA_PI_3[1]);
    let e1_x = q12.c0.c0.conjugate_in_place().to_owned();

    Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e2_x, e2_y, e2_z))
}

#[cfg(test)]
mod test {
    use super::*;
    use std::ops::Deref;

    #[test]
    fn test_beta_pi() {
        for x in constant::BETA_PI_1.deref() {
            println!("beta_pi_1: {:?}", x.to_string());
        }
        println!("");
        for x in constant::BETA_PI_1.deref() {
            println!("beta_pi_2: {:?}", x.to_string());
        }
        println!("");
        for x in constant::BETA_PI_3.deref() {
            println!("beta_pi_3: {:?}", x.to_string());
        }
    }
}
