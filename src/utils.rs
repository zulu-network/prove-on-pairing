use ark_bn254::Fr;
use ark_ff::One;
use num_bigint::BigUint;
use num_traits::{FromPrimitive, ToPrimitive};

pub fn biguint_to_naf(mut num: BigUint) -> Vec<i32> {
    let mut naf = Vec::new();
    let mut last_bit = 0;

    while num > BigUint::ZERO {
        let least_significant_bit = num.clone() & BigUint::from(1u8);
        let digit = least_significant_bit.pow(last_bit).to_i32().unwrap();
        naf.push(digit);
        last_bit = least_significant_bit.to_u32().unwrap();
        num >>= 1;
    }

    naf
}

pub fn naf_to_biguint(naf: &[i32]) -> BigUint {
    let mut val = BigUint::ZERO;

    for (i, digit) in naf.iter().enumerate() {
        // Shift left by the current digit's position
        val <<= 1;
        if *digit != 0 {
            // Add 1 for positive digit, 2 for negative
            val += BigUint::from_i32(if *digit > 0 { 1 } else { 2 }).unwrap();
        }
    }

    val
}

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
    use num_traits::One;

    #[test]
    fn test_naf_to_biguint() {
        let one = BigUint::one();

        let naf = biguint_to_naf(one);

        let actual = naf_to_biguint(&naf);

        assert_eq!(BigUint::one(), actual);
    }
}
