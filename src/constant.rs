use ark_bn254::{Config, Fq, Fq2, FqConfig, G1Affine, G2Affine};
use ark_ec::bn::BnConfig;
use ark_ff::MontFp;
use num_bigint::{BigInt, BigUint};
use once_cell::sync::Lazy;
use std::clone::Clone;
use std::ops::{Add, Mul};
use std::str::FromStr;

// pub static  BETA: OnceCell<Fq2> = {
//         OnceCell::from(Fq2::new(Fq::from(1), Fq::from(9)))
// };

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


// const X: &'static [u64] = &[4965661367192848881]. See more on: Config::X
// pub const X: BigUint = BigUint::from(4965661367192848881);
// // e = 6x + 2;
// pub const E: Lazy<BigUint> = Lazy::new(||  BigUint::from(6).mul(X.clone()).add(2));
