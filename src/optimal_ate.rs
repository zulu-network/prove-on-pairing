use crate::utils::to_naf;
use ark_bn254::g2::Config;
use ark_bn254::{Bn254, Fq12, Fq2, Fq2Config, Fq6, G1Affine, G1Projective, G2Affine, G2Projective};
use ark_ec::pairing::Pairing;
use ark_ec::short_weierstrass::Projective;
use ark_ec::AffineRepr;
use ark_ff::{Field, Fp2ConfigWrapper, QuadExtField};
use std::ops::{Add, Mul, Neg, Sub};

pub fn line_func_add(
    r: G2Projective,
    p: G2Projective,
    q: G2Projective,
    r2: Fq2,
) -> (Fq2, Fq2, Fq2, G2Projective) {
    let r_t = r.z.square();
    let B = p.x * r_t;
    let mut D = p.y + r.z;
    D = D.square();
    D -= r2;
    D -= r_t;
    D *= r_t;

    let H = B - r.x;
    let I = H.square();

    let E = I.double().double();

    let J = H * E;
    let mut L1 = D - r.y;
    L1 -= r.y;

    let V = r.x * E;

    let mut r_x = L1.square();
    r_x -= J;
    r_x -= V.double();

    let mut r_z = r.z + H;
    r_z = r_z.square();
    r_z -= r_t;
    r_z -= I;

    let mut t = V - r_x;
    t *= L1;
    let mut t2 = r.y * J;
    t2 = t2.double();
    let r_y = t - t2;

    let r_out = G2Projective::new(r_x, r_y, r_z);

    let mut t = p.y + r_z;
    t = t.square();
    t = t - r2;
    t = t - (r_z.square());

    let mut t2 = L1 * p.x;
    t2 = t2.double();

    let mut a = t2 - t;

    let x = q.y;
    let mut c = r_z.mul(q.y).double();

    let mut b = L1.neg();
    b = b.mul(q.x).double();

    // abandon the convenience of projective coordinate, be consistent with verifier
    // 2 * z_r
    let aux_inv = r_z.double().inverse().unwrap();

    let a = a.mul(aux_inv);
    let b = b.mul(aux_inv);
    let c = c.mul(aux_inv);

    (a, b, c, r_out)
}

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

pub fn miller(Q: G2Affine, P: G1Affine) {
    let f = Fq12::new(Fq6::ZERO, Fq6::ONE);

    let T = Q.clone();

    let Qp = Q.y.square();

    // 6x + 2 in NAF
    // let mut naf_6xp2 = to_naf();
    // naf_6xp2.reverse();
    // naf_6xp2.remove(0);
    //
    // let mut f_list = vec![];
    //
    // naf_6xp2.iter().enumerate().map(|a| {
    //     let fi = f.square();
    // });
}
