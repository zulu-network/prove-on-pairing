// def conjugate_of(self):
// """
//         For gamma = A + iB \in gfp2
//         gamma^p = A - iB
//         """
// return Fp2(self.x.additive_inverse(), self.y)
// pub fn fq2_conjugate_of()

use crate::constant::BETA;
use ark_bn254::{Fq, Fq12, Fq2, Fq6, FqConfig};
use ark_ff::{BigInt, BigInteger, Field, MontConfig};
use num_bigint::BigUint;
use std::ops::Mul;

// const beta: Fq2 = Fq2::new(Fq::ONE, Fq::from(9));
pub struct Fq12Ext;

impl Fq12Ext {
    pub fn beta_pi_1() -> Vec<Fq2> {
        let mut res = vec![];
        for i in 1..6 {
            // exp = i * ((module-1)/6)
            let mut t = FqConfig::MODULUS;
            t.sub_with_borrow(&BigInt::one());
            t.divn(6);
            t.muln(i);
            let exp = t;
            let pi = BETA.pow(exp);
            res.push(pi);
        }
        res
    }

    pub fn beta_pi_2() -> Vec<Fq2> {
        let mut res = vec![];
        for i in 1..6 {
            // exp = i * ((module^2 -1)/6)
            let mut t = FqConfig::MODULUS;
            // todo: add module's pow t= t^2;
            t.sub_with_borrow(&BigInt::one());
            t.divn(6);
            t.muln(i);
            let exp = t;
            let pi = BETA.pow(exp);
            res.push(pi);
        }
        res
    }
    pub fn beta_pi_3() -> Vec<Fq2> {
        let mut res = vec![];
        for i in 1..6 {
            // exp = i * ((module^2 -1)/6)
            let mut a = FqConfig::MODULUS;
            let mut t = FqConfig::MODULUS;
            // todo: add module's pow t=t^2;
            t.sub_with_borrow(&BigInt::one());
            t.divn(6);
            t.muln(i);
            let exp = t;
            let pi = BETA.pow(exp);
            res.push(pi);
        }
        res
    }
}

// const beta: Fq2 = Fq2::new(Fq::ONE, Fq::from(9));

// TODO: check if origin q12 is changed @Payne
pub fn fq12_to_frobenius(mut q12: Fq12) -> Fq12 {
    let e1_x = q12.c0.c0.conjugate_in_place().mul(Fq12Ext::beta_pi_1()[4]);
    let e1_y = q12.c0.c1.conjugate_in_place().mul(Fq12Ext::beta_pi_1()[2]);
    let e1_z = q12.c0.c2.conjugate_in_place().mul(Fq12Ext::beta_pi_1()[0]);

    let e2_x = q12.c1.c0.conjugate_in_place().mul(Fq12Ext::beta_pi_1()[3]);
    let e2_y = q12.c1.c1.conjugate_in_place().mul(Fq12Ext::beta_pi_1()[1]);
    let e2_z = q12.c1.c2.conjugate_in_place().to_owned();

    Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e2_x, e2_y, e2_z))
}

pub fn fq12_to_frobenius_p2(q12: Fq12) -> Fq12 {
    let e1_x = q12.c0.c0 * Fq12Ext::beta_pi_2()[4];
    let e1_y = q12.c0.c1 * Fq12Ext::beta_pi_2()[2];
    let e1_z = q12.c0.c2 * Fq12Ext::beta_pi_2()[0];

    let e2_x = q12.c1.c0 * Fq12Ext::beta_pi_2()[3];
    let e2_y = q12.c1.c1 * Fq12Ext::beta_pi_2()[1];
    let e2_z = q12.c1.c2;

    Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e2_x, e2_y, e2_z))
}

pub fn fq12_to_frobenius_p3(mut q12: Fq12) -> Fq12 {
    let e1_x = q12.c0.c0.conjugate_in_place().mul(Fq12Ext::beta_pi_3()[4]);
    let e1_y = q12.c0.c1.conjugate_in_place().mul(Fq12Ext::beta_pi_3()[2]);
    let e1_z = q12.c0.c2.conjugate_in_place().mul(Fq12Ext::beta_pi_3()[0]);

    let e2_x = q12.c1.c0.conjugate_in_place().mul(Fq12Ext::beta_pi_3()[3]);
    let e2_y = q12.c1.c1.conjugate_in_place().mul(Fq12Ext::beta_pi_3()[1]);
    let e2_z = q12.c1.c2.conjugate_in_place().to_owned();

    Fq12::new(Fq6::new(e1_x, e1_y, e1_z), Fq6::new(e2_x, e2_y, e2_z))
}
