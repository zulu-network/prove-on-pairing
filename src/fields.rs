// def conjugate_of(self):
// """
//         For gamma = A + iB \in gfp2
//         gamma^p = A - iB
//         """
// return Fp2(self.x.additive_inverse(), self.y)
// pub fn fq2_conjugate_of()

use crate::constant;
use crate::constant::BETA;
use ark_bn254::{Fq12, Fq2, Fq6};
use ark_ff::Field;
use num_bigint::BigUint;
use num_traits::{FromPrimitive, One, Pow};
use std::ops::{DivAssign, Mul, MulAssign, SubAssign};

// const beta: Fq2 = Fq2::new(Fq::ONE, Fq::from(9));
pub struct Fq12Ext;

impl Fq12Ext {
    pub fn beta_pi_1() -> Vec<Fq2> {
        let mut res = vec![];
        for i in 1..6 {
            // exp = i * ((module-1)/6)
            let mut t = constant::MODULUS.clone();
            t.sub_assign(BigUint::one());
            t.div_assign(BigUint::from_i32(6).unwrap());
            t.mul_assign(BigUint::from_i32(i).unwrap());
            let exp = t;
            let pi = BETA.pow(exp.to_u64_digits());
            res.push(pi);
        }
        res
    }

    pub fn beta_pi_2() -> Vec<Fq2> {
        let mut res = vec![];
        for i in 1..6 {
            // exp = i * ((module^2 -1)/6)
            let mut t = constant::MODULUS.clone();
            t = t.pow(2_u32);
            t.sub_assign(BigUint::one());
            t.div_assign(BigUint::from_i32(6).unwrap());
            t.mul_assign(BigUint::from_i32(i).unwrap());
            let exp = t;
            let pi = BETA.pow(exp.to_u64_digits());
            res.push(pi);
        }
        res
    }
    pub fn beta_pi_3() -> Vec<Fq2> {
        let mut res = vec![];
        for i in 1..6 {
            // exp = i * ((module^3 -1)/6)
            let mut t = constant::MODULUS.clone();
            t = t.pow(3_u32);
            t.sub_assign(BigUint::one());
            t.div_assign(BigUint::from_i32(6).unwrap());
            t.mul_assign(BigUint::from_i32(i).unwrap());
            let exp = t;
            let pi = BETA.pow(exp.to_u64_digits());
            res.push(pi);
        }
        res
    }
}

pub fn fq12_to_frobenius(mut q12: Fq12) -> Fq12 {
    // for i in 0..5 {
    //     println!(
    //         "in frobenius, beta_pi_1[{}] = {:?}\n\n",
    //         i,
    //         Fq12Ext::beta_pi_1()[i].to_string()
    //     );
    // let r: BigInt<4> = FqConfig::MODULUS;
    // println!("FqConfig::MODULUS = {:?}", r.to_string());
    // }

    // println!("frobenius Self = {:?}\n\n", q12.to_string());

    let e2_z = q12.c1.c2.conjugate_in_place().mul(Fq12Ext::beta_pi_1()[4]);
    let e2_y = q12.c1.c1.conjugate_in_place().mul(Fq12Ext::beta_pi_1()[2]);
    let e2_x = q12.c1.c0.conjugate_in_place().mul(Fq12Ext::beta_pi_1()[0]);

    let e1_z = q12.c0.c2.conjugate_in_place().mul(Fq12Ext::beta_pi_1()[3]);
    let e1_y = q12.c0.c1.conjugate_in_place().mul(Fq12Ext::beta_pi_1()[1]);
    let e1_x = q12.c0.c0.conjugate_in_place().to_owned();

    // println!("frobenius e1_x = {:?}\n\n", e1_x.to_string());
    // println!("frobenius e1_y = {:?}\n\n", e1_y.to_string());
    // println!("frobenius e1_z = {:?}\n\n", e1_z.to_string());
    // println!("frobenius e2_x = {:?}\n\n", e2_x.to_string());
    // println!("frobenius e2_y = {:?}\n\n", e2_y.to_string());
    // println!("frobenius e2_z = {:?}\n\n", e2_z.to_string());

    Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e2_x, e2_y, e2_z))
}

pub fn fq12_to_frobenius_p2(q12: Fq12) -> Fq12 {
    // for i in 0..5 {
    //     println!(
    //         "in frobenius, beta_pi_2[{}] = {:?}\n\n",
    //         i,
    //         Fq12Ext::beta_pi_2()[i].to_string()
    //     );
    //     // let r: BigInt<4> = FqConfig::MODULUS;
    //     // println!("FqConfig::MODULUS = {:?}", r.to_string());
    // }

    // println!("frobenius_p2 Self = {:?}\n\n", q12.to_string());

    let e2_z = q12.c1.c2 * Fq12Ext::beta_pi_2()[4];
    let e2_y = q12.c1.c1 * Fq12Ext::beta_pi_2()[2];
    let e2_x = q12.c1.c0 * Fq12Ext::beta_pi_2()[0];

    let e1_z = q12.c0.c2 * Fq12Ext::beta_pi_2()[3];
    let e1_y = q12.c0.c1 * Fq12Ext::beta_pi_2()[1];
    let e1_x = q12.c0.c0;

    // println!("frobenius_p2 e1_x = {:?}\n\n", e1_x.to_string());
    // println!("frobenius_p2 e1_y = {:?}\n\n", e1_y.to_string());
    // println!("frobenius_p2 e1_z = {:?}\n\n", e1_z.to_string());
    // println!("frobenius_p2 e2_x = {:?}\n\n", e2_x.to_string());
    // println!("frobenius_p2 e2_y = {:?}\n\n", e2_y.to_string());
    // println!("frobenius_p2 e2_z = {:?}\n\n", e2_z.to_string());

    Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e2_x, e2_y, e2_z))
}

pub fn fq12_to_frobenius_p3(mut q12: Fq12) -> Fq12 {
    // for i in 0..5 {
    //     println!(
    //         "in frobenius, beta_pi_3[{}] = {:?}\n\n",
    //         i,
    //         Fq12Ext::beta_pi_3()[i].to_string()
    //     );
    //     // let r: BigInt<4> = FqConfig::MODULUS;
    //     // println!("FqConfig::MODULUS = {:?}", r.to_string());
    // }

    let e2_z = q12.c1.c2.conjugate_in_place().mul(Fq12Ext::beta_pi_3()[4]);
    let e2_y = q12.c1.c1.conjugate_in_place().mul(Fq12Ext::beta_pi_3()[2]);
    let e2_x = q12.c1.c0.conjugate_in_place().mul(Fq12Ext::beta_pi_3()[0]);
    let e1_z = q12.c0.c2.conjugate_in_place().mul(Fq12Ext::beta_pi_3()[3]);
    let e1_y = q12.c0.c1.conjugate_in_place().mul(Fq12Ext::beta_pi_3()[1]);
    let e1_x = q12.c0.c0.conjugate_in_place().to_owned();

    // println!("frobenius_p3 e1_x = {:?}\n\n", e1_x.to_string());
    // println!("frobenius_p3 e1_y = {:?}\n\n", e1_y.to_string());
    // println!("frobenius_p3 e1_z = {:?}\n\n", e1_z.to_string());
    // println!("frobenius_p3 e2_x = {:?}\n\n", e2_x.to_string());
    // println!("frobenius_p3 e2_y = {:?}\n\n", e2_y.to_string());
    // println!("frobenius_p3 e2_z = {:?}\n\n", e2_z.to_string());

    Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e2_x, e2_y, e2_z))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_beta_pi() {
        for x in Fq12Ext::beta_pi_1() {
            println!("beta_pi_1: {:?}", x.to_string());
        }
        println!("");
        for x in Fq12Ext::beta_pi_2() {
            println!("beta_pi_2: {:?}", x.to_string());
        }
        println!("");
        for x in Fq12Ext::beta_pi_3() {
            println!("beta_pi_3: {:?}", x.to_string());
        }
    }
}
