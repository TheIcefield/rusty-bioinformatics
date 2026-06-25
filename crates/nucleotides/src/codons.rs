use std::{fmt::Display, str::FromStr};

use crate::nucleotide::Nucleotide;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Codon {
    A,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    K,
    L,
    M,
    N,
    P,
    Q,
    R,
    S,
    T,
    V,
    W,
    Y,
    Stop,
}

pub struct CodonSequence(pub Vec<Codon>);

pub struct ProteinString(pub CodonSequence);

impl Display for Codon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::F => "F",
                Self::L => "L",
                Self::S => "S",
                Self::Y => "Y",
                Self::C => "C",
                Self::W => "W",
                Self::I => "I",
                Self::M => "M",
                Self::T => "T",
                Self::P => "P",
                Self::H => "H",
                Self::Q => "Q",
                Self::R => "R",
                Self::N => "N",
                Self::K => "K",
                Self::V => "V",
                Self::A => "A",
                Self::D => "D",
                Self::E => "E",
                Self::G => "G",
                Self::Stop => "Stop",
            }
        )
    }
}

impl TryFrom<char> for Codon {
    type Error = String;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value.to_ascii_lowercase() {
            'a' => Ok(Self::A),
            'c' => Ok(Self::C),
            'd' => Ok(Self::D),
            'e' => Ok(Self::E),
            'f' => Ok(Self::F),
            'g' => Ok(Self::G),
            'h' => Ok(Self::H),
            'i' => Ok(Self::I),
            'k' => Ok(Self::K),
            'l' => Ok(Self::L),
            'm' => Ok(Self::M),
            'n' => Ok(Self::N),
            'p' => Ok(Self::P),
            'q' => Ok(Self::Q),
            'r' => Ok(Self::R),
            's' => Ok(Self::S),
            't' => Ok(Self::T),
            'v' => Ok(Self::V),
            'w' => Ok(Self::W),
            'y' => Ok(Self::Y),
            _ => Err(format!("\"{value}\" is not a codon")),
        }
    }
}

impl Codon {
    pub fn try_from_triplet(triplet: [Nucleotide; 3]) -> Option<Self> {
        for (entry_triplet, entry_codon) in RNA_CODON_TABLE.iter() {
            if *entry_triplet == triplet {
                return Some(*entry_codon);
            }
        }

        None
    }

    pub fn get_monoisotopic_mass(&self) -> Option<f64> {
        CODON_MONOISOTOPIC_MASS_TABLE
            .iter()
            .find(|cand| cand.0 == *self)
            .map(|cand| cand.1)
    }
}

impl FromStr for CodonSequence {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut codons = Vec::<Codon>::with_capacity(s.len());
        for ch in s.chars() {
            codons.push(Codon::try_from(ch)?);
        }

        Ok(Self(codons))
    }
}

impl From<&[Nucleotide]> for CodonSequence {
    fn from(seq: &[Nucleotide]) -> Self {
        let mut codons = Vec::new();

        for start in (0..seq.len()).step_by(3) {
            let Some(triplet) = seq
                .get(start..(start + 3))
                .map(|slice| [slice[0], slice[1], slice[2]])
            else {
                continue;
            };

            let Some(codon) = Codon::try_from_triplet(triplet) else {
                continue;
            };

            if codon == Codon::Stop {
                continue;
            }

            codons.push(codon);
        }

        Self(codons)
    }
}

impl CodonSequence {
    pub fn to_string(seq: &[Codon]) -> String {
        seq.iter().map(|c| c.to_string()).collect()
    }

    pub fn get_monoisotopic_mass(&self) -> f64 {
        Self::get_monoisotopic_mass_of_seq(&self.0)
    }

    fn get_monoisotopic_mass_of_seq(seq: &[Codon]) -> f64 {
        seq.iter().filter_map(|c| c.get_monoisotopic_mass()).sum()
    }
}

const RNA_CODON_TABLE: [([Nucleotide; 3], Codon); 16 * 4] = [
    // Uracil combination
    ([Nucleotide::U, Nucleotide::U, Nucleotide::U], Codon::F),
    ([Nucleotide::U, Nucleotide::U, Nucleotide::C], Codon::F),
    ([Nucleotide::U, Nucleotide::U, Nucleotide::A], Codon::L),
    ([Nucleotide::U, Nucleotide::U, Nucleotide::G], Codon::L),
    ([Nucleotide::U, Nucleotide::C, Nucleotide::U], Codon::S),
    ([Nucleotide::U, Nucleotide::C, Nucleotide::C], Codon::S),
    ([Nucleotide::U, Nucleotide::C, Nucleotide::A], Codon::S),
    ([Nucleotide::U, Nucleotide::C, Nucleotide::G], Codon::S),
    ([Nucleotide::U, Nucleotide::A, Nucleotide::U], Codon::Y),
    ([Nucleotide::U, Nucleotide::A, Nucleotide::C], Codon::Y),
    ([Nucleotide::U, Nucleotide::A, Nucleotide::A], Codon::Stop),
    ([Nucleotide::U, Nucleotide::A, Nucleotide::G], Codon::Stop),
    ([Nucleotide::U, Nucleotide::G, Nucleotide::U], Codon::C),
    ([Nucleotide::U, Nucleotide::G, Nucleotide::C], Codon::C),
    ([Nucleotide::U, Nucleotide::G, Nucleotide::A], Codon::Stop),
    ([Nucleotide::U, Nucleotide::G, Nucleotide::G], Codon::W),
    // Cytosine combination
    ([Nucleotide::C, Nucleotide::U, Nucleotide::U], Codon::L),
    ([Nucleotide::C, Nucleotide::U, Nucleotide::C], Codon::L),
    ([Nucleotide::C, Nucleotide::U, Nucleotide::A], Codon::L),
    ([Nucleotide::C, Nucleotide::U, Nucleotide::G], Codon::L),
    ([Nucleotide::C, Nucleotide::C, Nucleotide::U], Codon::P),
    ([Nucleotide::C, Nucleotide::C, Nucleotide::C], Codon::P),
    ([Nucleotide::C, Nucleotide::C, Nucleotide::A], Codon::P),
    ([Nucleotide::C, Nucleotide::C, Nucleotide::G], Codon::P),
    ([Nucleotide::C, Nucleotide::A, Nucleotide::U], Codon::H),
    ([Nucleotide::C, Nucleotide::A, Nucleotide::C], Codon::H),
    ([Nucleotide::C, Nucleotide::A, Nucleotide::A], Codon::Q),
    ([Nucleotide::C, Nucleotide::A, Nucleotide::G], Codon::Q),
    ([Nucleotide::C, Nucleotide::G, Nucleotide::U], Codon::R),
    ([Nucleotide::C, Nucleotide::G, Nucleotide::C], Codon::R),
    ([Nucleotide::C, Nucleotide::G, Nucleotide::A], Codon::R),
    ([Nucleotide::C, Nucleotide::G, Nucleotide::G], Codon::R),
    // Adenine combination
    ([Nucleotide::A, Nucleotide::U, Nucleotide::U], Codon::I),
    ([Nucleotide::A, Nucleotide::U, Nucleotide::C], Codon::I),
    ([Nucleotide::A, Nucleotide::U, Nucleotide::A], Codon::I),
    ([Nucleotide::A, Nucleotide::U, Nucleotide::G], Codon::M),
    ([Nucleotide::A, Nucleotide::C, Nucleotide::U], Codon::T),
    ([Nucleotide::A, Nucleotide::C, Nucleotide::C], Codon::T),
    ([Nucleotide::A, Nucleotide::C, Nucleotide::A], Codon::T),
    ([Nucleotide::A, Nucleotide::C, Nucleotide::G], Codon::T),
    ([Nucleotide::A, Nucleotide::A, Nucleotide::U], Codon::N),
    ([Nucleotide::A, Nucleotide::A, Nucleotide::C], Codon::N),
    ([Nucleotide::A, Nucleotide::A, Nucleotide::A], Codon::K),
    ([Nucleotide::A, Nucleotide::A, Nucleotide::G], Codon::K),
    ([Nucleotide::A, Nucleotide::G, Nucleotide::U], Codon::S),
    ([Nucleotide::A, Nucleotide::G, Nucleotide::C], Codon::S),
    ([Nucleotide::A, Nucleotide::G, Nucleotide::A], Codon::R),
    ([Nucleotide::A, Nucleotide::G, Nucleotide::G], Codon::R),
    // Guanine combination
    ([Nucleotide::G, Nucleotide::U, Nucleotide::U], Codon::V),
    ([Nucleotide::G, Nucleotide::U, Nucleotide::A], Codon::V),
    ([Nucleotide::G, Nucleotide::U, Nucleotide::C], Codon::V),
    ([Nucleotide::G, Nucleotide::U, Nucleotide::G], Codon::V),
    ([Nucleotide::G, Nucleotide::C, Nucleotide::U], Codon::A),
    ([Nucleotide::G, Nucleotide::C, Nucleotide::C], Codon::A),
    ([Nucleotide::G, Nucleotide::C, Nucleotide::A], Codon::A),
    ([Nucleotide::G, Nucleotide::C, Nucleotide::G], Codon::A),
    ([Nucleotide::G, Nucleotide::A, Nucleotide::U], Codon::D),
    ([Nucleotide::G, Nucleotide::A, Nucleotide::C], Codon::D),
    ([Nucleotide::G, Nucleotide::A, Nucleotide::A], Codon::E),
    ([Nucleotide::G, Nucleotide::A, Nucleotide::G], Codon::E),
    ([Nucleotide::G, Nucleotide::G, Nucleotide::U], Codon::G),
    ([Nucleotide::G, Nucleotide::G, Nucleotide::C], Codon::G),
    ([Nucleotide::G, Nucleotide::G, Nucleotide::A], Codon::G),
    ([Nucleotide::G, Nucleotide::G, Nucleotide::G], Codon::G),
];

const CODON_MONOISOTOPIC_MASS_TABLE: [(Codon, f64); 20] = [
    (Codon::A, 71.03711),
    (Codon::C, 103.00919),
    (Codon::D, 115.02694),
    (Codon::E, 129.04259),
    (Codon::F, 147.06841),
    (Codon::G, 57.02146),
    (Codon::H, 137.05891),
    (Codon::I, 113.08406),
    (Codon::K, 128.09496),
    (Codon::L, 113.08406),
    (Codon::M, 131.04049),
    (Codon::N, 114.04293),
    (Codon::P, 97.05276),
    (Codon::Q, 128.05858),
    (Codon::R, 156.10111),
    (Codon::S, 87.03203),
    (Codon::T, 101.04768),
    (Codon::V, 99.06841),
    (Codon::W, 186.07931),
    (Codon::Y, 163.06333),
];

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use crate::sequence::Sequence;

    use super::*;

    #[test]
    fn rna_translation_test() {
        // Given
        const RAW_DATA: &str = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA";

        // When
        let seq = Sequence::from_str(RAW_DATA).unwrap();
        let translated = seq.translate().unwrap();

        // Then
        const EXPECTED: &str = "MAMAPRTEINSTRING";
        assert_eq!(CodonSequence::to_string(&translated.0), EXPECTED);
    }

    #[test]
    fn monoisotopic_mass_test() {
        // Given
        const RAW_DATA: &str = "SKADYEK";

        // When
        let seq = CodonSequence::from_str(RAW_DATA).unwrap();
        let mass = seq.get_monoisotopic_mass();

        // Then
        assert_eq!(mass.round(), 821.0);
    }
}
