use std::{fmt::Display, str::FromStr};

use crate::nucleotide::Nucleotide;

#[derive(Default, Clone)]
pub struct Sequence(pub Vec<Nucleotide>);

impl Display for Sequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for n in &self.0 {
            write!(f, "{n}")?;
        }
        Ok(())
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum SequenceKind {
    Dna,
    Rna,
}

impl Display for SequenceKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Dna => write!(f, "DNA"),
            Self::Rna => write!(f, "RNA"),
        }
    }
}

impl FromStr for Sequence {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut nucleotides = Vec::<Nucleotide>::with_capacity(s.len());
        for ch in s.chars() {
            nucleotides.push(Nucleotide::try_from(ch)?);
        }

        Ok(Sequence(nucleotides))
    }
}

impl Sequence {
    pub fn is_dna(&self) -> bool {
        self.get_kind() == SequenceKind::Dna
    }

    pub fn get_kind(&self) -> SequenceKind {
        if self.count(Nucleotide::Uracil) == 0 {
            SequenceKind::Dna
        } else {
            SequenceKind::Rna
        }
    }

    /// Counts Guanine and Cytosine nucleotides
    pub fn get_gc_content(&self) -> f64 {
        let gc = self
            .0
            .iter()
            .filter(|n| matches!(n, Nucleotide::Guanine | Nucleotide::Cytosine)) // yeild only G and C
            .count();

        (gc as f64) / (self.0.len() as f64) * 100.0
    }

    /// Returns count of nucleotides in sequence
    pub fn length(&self) -> usize {
        self.0.len()
    }

    /// Returns count of concrete nucleotyde in sequence
    pub fn count(&self, kind: Nucleotide) -> usize {
        self.0.iter().filter(|found| **found == kind).count()
    }

    /// Replaces Thymine to Uracil. Kind of sequnce is changed to RNA
    /// Applicable to: DNA sequence kind
    pub fn transcribe(&self) -> Result<Self, String> {
        if !self.is_dna() {
            return Err(format!(
                "Only DNA can be transcribed. Sequence kind is {}",
                self.get_kind()
            ));
        }

        let mut transcribed = self.clone();
        transcribed.0.iter_mut().for_each(|n| {
            if *n == Nucleotide::Thymine {
                *n = Nucleotide::Uracil
            }
        });

        Ok(transcribed)
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
        let dna = Sequence::from_str(raw_data).unwrap();

        // Then
        assert_eq!(dna.get_gc_content(), 44.44444444444444);
    }

    #[test]
    fn dna_transcribe_test() {
        // Given
        let raw_data = "ATGGCCTAA";

        // When
        let dna = Sequence::from_str(raw_data).unwrap();
        let rna = dna.transcribe().unwrap();

        // Then
        assert_eq!(rna.to_string(), "AUGGCCUAA");
    }
}
