#[cfg(test)]
mod tests {
    use std::{
        fs::File,
        io::{self, Read},
        str::FromStr,
    };

    use ark_bn254::{Fq, Fq12, Fq2, Fq6};
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ff::{Field, MontFp, One};
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;

    use crate::{
        constant::{g1, g2},
        pairing_verify::Pairing,
    };

    #[test]
    fn test_pairing_verify() {
        let p1 = g1
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();
        let p2 = g1.mul_bigint(BigUint::one().to_u64_digits()).into_affine();

        let q1 = g2.mul_bigint(BigUint::one().to_u64_digits()).into_affine();
        let q2 = g2
            .mul_bigint(BigUint::from_i8(3).unwrap().to_u64_digits())
            .into_affine();

        println!("Q1: {:?},\nQ2: {:?}", q1, q2);
        let l1 = generate_lines("l1.txt");
        let l2 = generate_lines("l2.txt");
        // println!("l1 = {:?}", l1);
        // println!("l2 = {:?}", l2);

        let e: i128 = 29793968203157093288;

        let cx = Fq6::new(
            Fq2::new(
                MontFp!(
                    "20201738709783866217331415010360981929515908053985546990877702762237532098619"
                ),
                MontFp!(
                    "15456279988053324345541467341779085106204653073005659040613522482047331531215"
                ),
            ),
            Fq2::new(
                MontFp!(
                    "13489608147144759151428672340723831207483960909985311771672648530728538547919"
                ),
                MontFp!(
                    "5839925900105875329292204483554224423437555932138073078712709683459700700695"
                ),
            ),
            Fq2::new(
                MontFp!(
                    "15126867401573936295269375558559019240807443361889083511177200289849932772721"
                ),
                MontFp!(
                    "17764970281802625026130810332747444901517554817814187548834331281542778642427"
                ),
            ),
        );
        let cy = Fq6::new(
            Fq2::new(
                MontFp!(
                    "4032083302069892693408840942613137942537568963287332711747870178718866954397"
                ),
                MontFp!(
                    "11300710317364955496310923778035288238945437277150320360850689389938305374038"
                ),
            ),
            Fq2::new(
                MontFp!(
                    "2164084145544483012866601008567729233848223519042058191386097738904720442710"
                ),
                MontFp!(
                    "9171884817913122022865873159746318373988003274101841511199162153768410371384"
                ),
            ),
            Fq2::new(
                MontFp!(
                    "13448531161494720028237301515094431586317095588653693447905281697278723034814"
                ),
                MontFp!(
                    "13906303506703881772293742616476616035816177707476145711555641424944353982125"
                ),
            ),
        );
        let c = Fq12::new(cx, cy);

        let c_inv = c.inverse().unwrap();

        let wy = Fq6::new(
            Fq2::ZERO,
            Fq2::new(
                MontFp!(
                    "14396857495183920858389897556068566964512197819080094036164390672375860233648"
                ),
                MontFp!(
                    "20412823004746074851973712366746204934010010499129554554879858850258371975807"
                ),
            ),
            Fq2::ZERO,
        );
        let wi = Fq12::new(Fq6::ZERO, wy);

        let verify_res = Pairing::verify_pairings(vec![p1, p2], &[l1, l2], e, c, c_inv, wi);
        println!("verify_res = {}", verify_res);
        assert_eq!(verify_res, Fq12::ONE);
    }

    fn generate_lines(file_path: &str) -> Vec<(Fq2, Fq2)> {
        let input = read_from_file(file_path).unwrap();
        let result = parse_to_vec(input.as_str());
        let mut res = vec![];
        for point in result {
            let x1 = Fq::from_str(point.0 .0.as_str()).unwrap();
            let x2 = Fq::from_str(point.0 .1.as_str()).unwrap();
            let y1 = Fq::from_str(point.1 .0.as_str()).unwrap();
            let y2 = Fq::from_str(point.1 .1.as_str()).unwrap();
            let x = Fq2::new(x1, x2);
            let y = Fq2::new(y1, y2);
            res.push((x, y));
        }
        res
    }

    fn parse_to_vec(input: &str) -> Vec<((String, String), (String, String))> {
        let cleaned_input = input.replace(&['[', ']', '(', ')', ' '][..], "");
        let numbers: Vec<&str> = cleaned_input.split(',').collect();

        numbers
            .chunks(4)
            .map(|chunk| {
                (
                    (chunk[0].to_string(), chunk[1].to_string()),
                    (chunk[2].to_string(), chunk[3].to_string()),
                )
            })
            .collect()
    }

    fn read_from_file(file_path: &str) -> io::Result<String> {
        let mut file = File::open(file_path)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;
        Ok(contents)
    }
}
