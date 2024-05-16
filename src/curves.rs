use ark_bn254::G1Projective;

impl G1Projective {
    pub fn force_affine(&self) {
        point_force_affine()
    }
}

pub fn point_force_affine() {}