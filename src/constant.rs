use ark_bn254::{Fq, Fq2};
use std::cell::OnceCell;

// pub static  BETA: OnceCell<Fq2> = {
//         OnceCell::from(Fq2::new(Fq::from(1), Fq::from(9)))
// };

pub fn get_BETA() -> Fq2 {
    Fq2::new(Fq::from(1), Fq::from(9))
}
