use ark_bn254::Fr;
use num_bigint::BigUint;

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

pub fn px(x: BigUint) {}
