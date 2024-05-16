use ark_bn254::{Fq12, Fq2, Fq6};
use ark_ff::Field;
use std::ops::{Add, Mul, Sub};

pub fn mul_line(r: Fq12, a: Fq2, b: Fq2, c: Fq2) -> Fq12 {
    // See function fp12e_mul_line in dclxvi
    let t1 = Fq6::new(Fq2::ZERO, a, b);
    let t2 = Fq6::new(Fq2::ZERO, a, b.add(c));

    let t1 = t1.mul(r.c0);
    let mut t3 = r.c1;
    t3.mul_by_fp2(&c);

    let y = t3;
    let mut x = r.c0.add(r.c1);
    x = x.mul(t2);
    x = x.sub(t1);
    let x = x.sub(y);
    // TODO: y = y.add(t1.mul_tau())
    // let y = y.add()

    Fq12::new(x, y)
}

pub fn mul_line_base(r: Fq12, a: Fq2, b: Fq2, c: Fq2) -> Fq12 {
    let fq6_1 = Fq6::new(Fq2::ZERO, a, b);
    let fq6_2 = Fq6::new(Fq2::ZERO, Fq2::ZERO, c);
    let fl = Fq12::new(fq6_1, fq6_2);
    r.mul(fl)
}
