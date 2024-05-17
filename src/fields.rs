// def conjugate_of(self):
// """
//         For gamma = A + iB \in gfp2
//         gamma^p = A - iB
//         """
// return Fp2(self.x.additive_inverse(), self.y)
// pub fn fq2_conjugate_of()

use crate::constant::get_BETA;
use ark_bn254::{Fq, Fq12, Fq2, FqConfig};
use ark_ff::{BigInt, BigInteger, Field, MontConfig};
use num_bigint::BigUint;

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

pub struct Fq12Ext;

impl Fq12Ext {
    fn beta_pi_1() -> Vec<Fq2> {
        let mut res = vec![];
        for i in 1..6 {
            // exp = i * ((module-1)/6)
            let mut t = FqConfig::MODULUS;
            t.sub_with_borrow(&BigInt::one());
            t.divn(6);
            t.muln(i);
            let exp = t;
            let pi = get_BETA().pow(exp);
            res.push(pi);
        }
        res
    }

    fn beta_pi_2() -> Vec<Fq2> {
        let mut res = vec![];
        for i in 1..6 {
            // exp = i * ((module^2 -1)/6)
            let mut t = FqConfig::MODULUS;
            // todo: add module's pow t= t^2;
            t.sub_with_borrow(&BigInt::one());
            t.divn(6);
            t.muln(i);
            let exp = t;
            let pi = get_BETA().pow(exp);
            res.push(pi);
        }
        res
    }
    fn beta_pi_3() -> Vec<Fq2> {
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
            let pi = get_BETA().pow(exp);
            res.push(pi);
        }
        res
    }
}
