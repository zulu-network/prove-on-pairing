use std::ops::Mul;
use ark_bn254::{Fq12, Fq2, Fq6};
use ark_ff::Field;


pub fn mul_line_base(r: Fq12, a:Fq2, b:Fq2 , c:Fq2 )-> Fq12{
    let fq6_1 = Fq6::new(Fq2::ZERO, a, b);
    let fq6_2 = Fq6::new(Fq2::ZERO, Fq2::ZERO, c);
    let fl = Fq12::new(fq6_1, fq6_2);
    r.mul(fl)
}
