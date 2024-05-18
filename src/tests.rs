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

        println!("P1: {:?},\nP2: {:?}", p1, p2);
        // println!("Q1: {:?},\nQ2: {:?}", q1, q2);
        let l1 = generate_lines("l1.txt");
        let l2 = generate_lines("l2.txt");
        // println!("l1 = {:?}", l1);
        // println!("l2 = {:?}", l2);

        let e: i128 = 29793968203157093288;

        let cy = Fq6::new(
            Fq2::new(
                MontFp!(
                    "1398074605395234385742746982635607890999901394449920787952679263329618900423"
                ),
                MontFp!(
                    "14863623346637387528754574561712611762280560946730003365167624628783031806851"
                ),
            ),
            Fq2::new(
                MontFp!(
                    "3943322538927771188455999564370373075332486995274784857622626372657519173057"
                ),
                MontFp!(
                    "3097735610337843900368979344756038697740814458812421595304050335303356161701"
                ),
            ),
            Fq2::new(
                MontFp!(
                    "15535495333518547635388053375087086844764140428602229212987121796906694947701"
                ),
                MontFp!(
                    "17900230360671217726773770977461976340138958556757112953925160338294259755451"
                ),
            ),
        );
        let cx = Fq6::new(
            Fq2::new(
                MontFp!(
                    "4667318771591831111896913998599388465660409602488670805853346418592518763334"
                ),
                MontFp!(
                    "6713632336462274726021105432467423937237553861088593973094350524276003798886"
                ),
            ),
            Fq2::new(
                MontFp!(
                    "21081939668623356886869294872217712850424216414807721597624545100857744884557"
                ),
                MontFp!(
                    "16588937522064281852235800901297894214318954372559985137225140148647705064928"
                ),
            ),
            Fq2::new(
                MontFp!(
                    "19903885137215268899583999685812289845712323958279995809868941246364604339090"
                ),
                MontFp!(
                    "15864462405647635712161095062692932104398676326285281473670528335466546619372"
                ),
            ),
        );
        let c = Fq12::new(cx, cy);

        let c_inv = c.inverse().unwrap();
        // println!("c_inv = {}", c_inv);

        // let cx = Fq6::new(
        //     Fq2::new(
        //         MontFp!(
        //             "8410820985486828179323384661117119350186451376490980649824012215447143383070"
        //         ),
        //         MontFp!(
        //             "14026349898266845133139930106306509595218002555547363018683961997134399146594"
        //         ),
        //     ),
        //     Fq2::new(
        //         MontFp!(
        //             "4010845789757773467267306813104728443780980948449336471499146318221939612736"
        //         ),
        //         MontFp!(
        //             "8981943063447879822320961732267243641159099802682360217961181068020761062833"
        //         ),
        //     ),
        //     Fq2::new(
        //         MontFp!(
        //             "21846041236494056335818375613098650413853674650084417248122340331839228071908"
        //         ),
        //         MontFp!(
        //             "12587353521189678719456917396579902327589368988358444919247389819839883246651"
        //         ),
        //     ),
        // );
        // let cy = Fq6::new(
        //     Fq2::new(
        //         MontFp!(
        //             "3793897967655062532148090808831468241380877204407587825177490572927068583760"
        //         ),
        //         MontFp!(
        //             "13795480458105548408116770738811294846253384391086666718941865371697734616207"
        //         ),
        //     ),
        //     Fq2::new(
        //         MontFp!(
        //             "1637095471437938380234191854670070198758411158651802935529746171229487189119"
        //         ),
        //         MontFp!(
        //             "4673921058059591577865555722224213869676984985353750255526409269282781480930"
        //         ),
        //     ),
        //     Fq2::new(
        //         MontFp!(
        //             "1777319391200362946206099411836734672046605317764214380952126362400716640765"
        //         ),
        //         MontFp!(
        //             "12075052841639945456894285024169683189507774633491517635567025828799169291349"
        //         ),
        //     ),
        // );
        // let c_inv = Fq12::new(cy, cx);

        // println!("c = {}", c);
        // println!("c_inv = {}", c_inv);

        assert_eq!(c * c_inv, Fq12::ONE);

        let wi_x = Fq6::new(
            Fq2::ZERO,
            Fq2::new(
                MontFp!(
                    "1552599724427260165122619260710620891632201030296138444657693101145235169998"
                ),
                MontFp!(
                    "6872310070170355934577608915679474652968694353285292038893436960198570441414"
                ),
            ),
            Fq2::ZERO,
        );
        let wi = Fq12::new(wi_x, Fq6::ZERO);

        // println!("wi = {}", wi);

        let verify_res = Pairing::verify_pairings(vec![p1, p2], &[l1, l2], e, c, c_inv, wi);
        // println!("verify_res = {}", verify_res);
        // assert_eq!(verify_res, Fq12::ONE);
    }

    fn generate_lines(file_path: &str) -> Vec<(Fq2, Fq2)> {
        let input = read_from_file(file_path).unwrap();
        let result = parse_to_vec(input.as_str());
        let mut res = vec![];
        for point in result {
            let x2 = Fq::from_str(point.0 .0.as_str()).unwrap();
            let x1 = Fq::from_str(point.0 .1.as_str()).unwrap();
            let y2 = Fq::from_str(point.1 .0.as_str()).unwrap();
            let y1 = Fq::from_str(point.1 .1.as_str()).unwrap();
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
