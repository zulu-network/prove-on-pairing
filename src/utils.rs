use ark_ff::One;
use num_bigint::BigUint;
use num_traits::{FromPrimitive, ToPrimitive};

pub fn biguint_to_naf(num: BigUint) -> Vec<i8> {
    to_naf(num.to_i128().unwrap())
}

// TODO: This should not be public.
pub fn to_naf(mut x: i128) -> Vec<i8> {
    let mut z = vec![];
    while x > 0 {
        if x % 2 == 0 {
            z.push(0);
        } else {
            let zi: i8 = 2 - (x % 4) as i8;
            x -= zi as i128;
            z.push(zi)
        }

        x = x / 2
    }
    return z;
}

pub fn px(x: BigUint) -> BigUint {
    let p1 = BigUint::from_i8(36).unwrap();
    let p2 = BigUint::from_i8(24).unwrap();
    p1.clone() * x.pow(4) + p1 * x.pow(3) + p2 * x.pow(2) + BigUint::one()
}

pub fn rx(x: BigUint) -> BigUint {
    let p1 = BigUint::from_i8(36).unwrap();
    let p2 = BigUint::from_i8(18).unwrap();
    p1.clone() * x.pow(4) + p1 * x.pow(3) + p2 * x.pow(2) + BigUint::one()
}

pub fn tx(x: BigUint) -> BigUint {
    let p1 = BigUint::from_i8(6).unwrap();
    p1 * x.pow(2) + BigUint::one()
}

pub fn hx(x: BigUint) -> BigUint {
    let p1 = BigUint::from_i8(18).unwrap();
    let numerator = p1 * x.pow(12) - BigUint::one();
    let denominator = rx(x);
    numerator / denominator
}

pub fn lambdax(x: BigUint) -> BigUint {
    let p1 = BigUint::from_i8(36).unwrap();
    let p2 = BigUint::from_i8(18).unwrap();
    p1.clone() * x.pow(4) + p1 * x.pow(3) + p2 * x.pow(2) + BigUint::one()
}

pub fn mx(x: BigUint) -> BigUint {
    let p1 = BigUint::from_i8(36).unwrap();
    let p2 = BigUint::from_i8(18).unwrap();
    p1.clone() * x.pow(4) + p1 * x.pow(3) + p2 * x.pow(2) + BigUint::one()
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::constant;
    use std::ops::Deref;

    #[test]
    fn test_biguint_naf() {
        let mut expect = to_naf(29793968203157093288);
        expect.reverse();
        expect.remove(0);
        println!("res: {:?}", expect);

        println!("E :{:?}", constant::E.deref());
        let mut actual = biguint_to_naf(constant::E.clone());
        actual.reverse();
        actual.remove(0);
        println!("res: {:?}", actual);

        assert_eq!(expect, actual);
    }
}
