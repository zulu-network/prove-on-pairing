pub fn to_naf(mut x: usize) -> Vec<usize> {
    let mut z = vec![];
    while x > 0 {
        if x % 2 == 0 {
            z.push(0);
        } else {
            let zi = 2 - (x % 4);
            x -= zi;
            z.push(zi)
        }

        x = x / 2
    }
    return z;
}
