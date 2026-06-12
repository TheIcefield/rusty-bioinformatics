use std::str::FromStr;

pub enum Nucleotyde {
    Adenine,
    Guanine,
    Cytosine,
    Thymine,
}

#[derive(Default)]
pub struct Nucleotides(pub Vec<Nucleotyde>);

pub type DNA = Nucleotides;
pub type RNA = Nucleotides;

impl TryFrom<char> for Nucleotyde {
    type Error = String;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value.to_ascii_lowercase() {
            'a' => Ok(Self::Adenine),
            'g' => Ok(Self::Guanine),
            'c' => Ok(Self::Cytosine),
            't' => Ok(Self::Thymine),
            ch => Err(format!("\"{ch}\" is not a nucleotide")),
        }
    }
}

impl FromStr for Nucleotides {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut nucleotides = Vec::<Nucleotyde>::with_capacity(s.len());
        for ch in s.chars() {
            let nucleotyde = Nucleotyde::try_from(ch)?;
            nucleotides.push(nucleotyde);
        }

        Ok(Nucleotides(nucleotides))
    }
}

impl Nucleotides {
    pub fn get_gc_content(&self) -> f64 {
        let gc = self
            .0
            .iter()
            .filter(|n| matches!(n, Nucleotyde::Guanine | Nucleotyde::Cytosine)) // yeild only G and C
            .count();

        (gc as f64) / (self.0.len() as f64) * 100.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_content_test() {
        // Given
        let raw_data = "ATGGCCTAA";

        // When
        let dna = DNA::from_str(raw_data).unwrap();

        // Then
        assert_eq!(dna.get_gc_content(), 44.44444444444444);
    }
}
