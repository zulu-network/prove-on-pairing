use ark_bn254::{Config, Fq, Fq2, FqConfig, G1Affine, G2Affine};
use ark_ec::bn::BnConfig;
use ark_ff::{MontFp, PrimeField};
use num_bigint::{BigInt, BigUint};
use num_traits::{FromPrimitive, Num};
use once_cell::sync::Lazy;
use std::clone::Clone;
use std::ops::{Add, Mul, Sub};
use std::str::FromStr;

pub static BETA: Lazy<Fq2> = Lazy::new(|| Fq2::new(1.into(), Fq::from(9)));

pub const G1_GENERATOR_X: Fq =
    MontFp!("19491323635986486980056165026003970884581302300479364565163758691834883767296");
pub const G1_GENERATOR_Y: Fq =
    MontFp!("2503817206389621232991390790939417031444960302945150474681637705185779211401");

pub const G2_GENERATOR_X: Fq2 = Fq2::new(
    MontFp!("10022529265301880767558967801827554994678953177337994173174782310334418209951"),
    MontFp!("11403269112307582471523194844678173363615200121780745962919905543513926078845"),
);
pub const G2_GENERATOR_Y: Fq2 = Fq2::new(
    MontFp!("14403937293889182757621054345090826401263455856569175761852807173588543872656"),
    MontFp!("7417909083002664933410862546938954664060641619680344911439335935535164894254"),
);

pub const g1: G1Affine = G1Affine::new_unchecked(G1_GENERATOR_X, G1_GENERATOR_Y);
pub const g2: G2Affine = G2Affine::new_unchecked(G2_GENERATOR_X, G2_GENERATOR_Y);

// constant modulus of Fq
pub const MODULUS: &str = "30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47";

// const X: &'static [u64] = &[4965661367192848881]. See more on: Config::X
pub const X: Lazy<BigUint> = Lazy::new(|| BigUint::from_i128(4965661367192848881).unwrap());

// e = 6X + 2;
pub const E: Lazy<BigUint> = Lazy::new(|| {
    let x6 = BigUint::from_i8(6).unwrap() * BigUint::from_i128(4965661367192848881).unwrap();
    // 6x + 2
    x6 + BigUint::from_i8(2).unwrap()
});

// optimal lambda in miller loop, lambda

pub const LAMBDA: Lazy<BigUint> = Lazy::new(|| lambda(&X));

// lambdax = 6x + 2 + p - p^2 + p^3
fn lambda(x: &BigUint) -> BigUint {
    println!("x: {:?}", x);
    let p = BigUint::from(Fq::MODULUS);
    let p_pow3 = p.pow(3_u32);
    let p_pow2 = p.pow(2_u32);
    let x6 = BigUint::from_i8(6).unwrap().mul(x);
    let two = BigUint::from_i8(2).unwrap();
    p_pow3.sub(p_pow2).add(p).add(two).add(x6)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_lamb() {
        let actual = lambda(&X);
        let expect = BigUint::from_str(
            "10486551571378427818905133077457505975146652579011797175399169355881771981095211883813744499745558409789005132135496770941292989421431235276221147148858384772096778432243207188878598198850276842458913349817007302752534892127325269"
        ).unwrap();

        assert_eq!(actual, expect);
    }

    #[test]
    fn test_modulus() {
        let actual = BigUint::from(Fq::MODULUS);
        let expect = BigUint::from_str_radix(MODULUS, 16).unwrap();
        assert_eq!(actual, expect);
    }

    #[test]
    fn test_e() {
        let actual = E.clone();
        let expect = BigUint::from_str("29793968203157093288").unwrap();
        assert_eq!(actual, expect);
    }
}
