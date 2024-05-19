use crate::constant;
use ark_bn254::{Fq, Fq2, G1Affine, G1Projective, G2Affine, G2Projective};
use ark_ec::{AffineRepr, CurveGroup, Group};
use ark_ff::{Field, MontFp, PrimeField};
use num_bigint::BigUint;
use num_traits::{FromPrimitive, Num, One, Pow};
use once_cell::sync::Lazy;
use std::ops::{Add, Deref, DivAssign, Mul, MulAssign, Sub, SubAssign};

// constant modulus of Fq
pub const MODULUS_STR: &str = "30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47";
pub const MODULUS: Lazy<BigUint> = Lazy::new(|| BigUint::from_str_radix(MODULUS_STR, 16).unwrap());

// const X: &'static [u64] = &[4965661367192848881]. See more on: Config::X
pub static X: Lazy<BigUint> = Lazy::new(|| BigUint::from_i128(4965661367192848881).unwrap());

// e = 6X + 2;
pub static E: Lazy<BigUint> = Lazy::new(|| {
    let x6 = BigUint::from_i8(6).unwrap() * BigUint::from_i128(4965661367192848881).unwrap();
    // 6x + 2
    x6 + BigUint::from_i8(2).unwrap()
});

pub static BETA: Lazy<Fq2> = Lazy::new(|| Fq2::new(Fq::from(9), Fq::ONE));

// optimal lambda in miller loop, lambda
pub const LAMBDA: Lazy<BigUint> = Lazy::new(|| {
    // lambdax = 6X + 2 + p - p^2 + p^3
    let p = MODULUS.clone();
    let p_pow2 = p.clone().pow(2_u32);
    let p_pow3 = p.clone().pow(3_u32);
    let x6 = BigUint::from_i8(6).unwrap().mul(X.deref());
    let two = BigUint::from_i8(2).unwrap();
    p_pow3.sub(p_pow2).add(p).add(two).add(x6)
});

pub static BETA_PI_1: Lazy<Vec<Fq2>> = Lazy::new(|| {
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
});

pub static BETA_PI_2: Lazy<Vec<Fq2>> = Lazy::new(|| {
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
});

pub static BETA_PI_3: Lazy<Vec<Fq2>> = Lazy::new(|| {
    let mut res = vec![];
    for i in 1..6 {
        // exp = i * ((module^3 -1)/6)
        let mut t = MODULUS.clone();
        t = t.pow(3_u32);
        t.sub_assign(BigUint::one());
        t.div_assign(BigUint::from_i32(6).unwrap());
        t.mul_assign(BigUint::from_i32(i).unwrap());
        let exp = t;
        let pi = BETA.pow(exp.to_u64_digits());
        res.push(pi);
    }
    res
});

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::G1Projective;
    use ark_ec::{CurveGroup, Group};
    use std::str::FromStr;

    #[test]
    fn test_lamb() {
        let actual = LAMBDA.clone();
        let expect = BigUint::from_str(
            "10486551571378427818905133077457505975146652579011797175399169355881771981095211883813744499745558409789005132135496770941292989421431235276221147148858384772096778432243207188878598198850276842458913349817007302752534892127325269"
        ).unwrap();

        assert_eq!(actual, expect);
    }

    #[test]
    fn test_modulus() {
        let actual = BigUint::from(Fq::MODULUS);
        let expect = BigUint::from_str_radix(MODULUS_STR, 16).unwrap();
        assert_eq!(actual, expect);
    }

    #[test]
    fn test_e() {
        let actual = E.clone();
        let expect = BigUint::from_str("29793968203157093288").unwrap();
        assert_eq!(actual, expect);
    }

    #[test]
    fn test_beta() {
        // passed
        println!("Beta: {:?}", BETA.clone());
    }
}
