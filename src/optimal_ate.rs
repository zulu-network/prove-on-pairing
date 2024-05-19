use ark_bn254::{Bn254, Fq12, Fq2, Fq2Config, Fq6, G1Affine, G1Projective, G2Affine, G2Projective};
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{Field, Fp2ConfigWrapper, QuadExtField};
use num_bigint::BigUint;
use std::ops::{Add, Mul, Neg, Sub};

pub fn line_func_add(
    r: G2Projective,
    p: G2Projective,
    q: G1Projective,
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
    let mut c = r_z;
    c.mul_assign_by_basefield(&q.y);
    c = c.double();

    let mut b = L1.neg();
    b.mul_assign_by_basefield(&q.x);
    b = b.double();

    // abandon the convenience of projective coordinate, be consistent with verifier
    // 2 * z_r
    let aux_inv = r_z.double().inverse().unwrap();

    let a = a.mul(aux_inv);
    let b = b.mul(aux_inv);
    let c = c.mul(aux_inv);

    (a, b, c, r_out)
}

fn line_func_double(r: G2Projective, q: G1Projective) -> (Fq2, Fq2, Fq2, G2Projective) {
    let r_t = r.z.square();

    let A = r.x.square();
    let B = r.y.square();
    let C = B.square();

    let mut D = r.x + B;
    D = D.square();
    D -= A;
    D -= C;
    D = D.double();

    let E = A.double() + A;
    let F = E.square();

    // C^8
    let C8 = C.double().double().double();

    let r_x = F - D.double();
    let r_y = E * (D - r_x) - C8;

    // (y+z)*(y+z) - (y*y) - (z*z) = 2*y*z
    let r_z = (r.y + r.z).square() - B - r_t;

    assert_eq!(r_z, r.y * r.z.double());

    let r_out = G2Projective::new(r_x, r_y, r_z);
    assert!(r_out.into_affine().is_on_curve());

    let mut a = r.x + E;
    a = a.square();
    a -= A + F + B.double().double();

    let mut t = E * r_t;
    t = t.double();
    let mut b = t.neg();
    b.mul_assign_by_basefield(&q.x);

    let mut c = r_z * r_t;
    c = c.double();
    c.mul_assign_by_basefield(&q.y);

    // abandon the convenience of projective coordinate, be consistant with verifier
    // 2 * z_r * z_t^2
    let mut aux_inv = r_t.mul(r_z).double().inverse().unwrap();
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
    let fq6_y = Fq6::new(b, a, Fq2::ZERO);
    let fq6_x = Fq6::new(c, Fq2::ZERO, Fq2::ZERO);
    let fl = Fq12::new(fq6_x, fq6_y);
    r.mul(fl)
}

#[cfg(test)]
mod test {
    use crate::constant::{g1, g2};
    use ark_bn254::{Bn254, G1Affine, G2Affine};
    use ark_ec::pairing::Pairing;
    use ark_ec::{AffineRepr, CurveGroup};
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;

    #[test]
    fn test_miller_loop() {
        let Q = g2
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();
        let P = g1
            .mul_bigint(BigUint::from_i8(4).unwrap().to_u64_digits())
            .into_affine();

        println!("P:{:?}", P);
        println!("\n Q:{:?}", Q);
        let actual = Bn254::miller_loop(
            <G1Affine as Into<<Bn254 as Pairing>::G1Prepared>>::into(P),
            <G2Affine as Into<<Bn254 as Pairing>::G2Prepared>>::into(Q),
        );

        println!("\n actual: {:?}", actual.0.to_string());
    }
}
