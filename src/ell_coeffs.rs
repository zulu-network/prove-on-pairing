use ark_bn254::Fq2;

// Porting from ark_ec::bn::g2::EllCoeff
pub type EllCoeff = (Fq2, Fq2, Fq2);
