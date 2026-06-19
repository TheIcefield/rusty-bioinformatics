use std::{fmt::Display, str::FromStr};

use crate::nucleotide::Nucleotide;

#[derive(Default, Clone)]
pub struct Sequence(pub Vec<Nucleotide>);

#[derive(Default, Clone)]
pub struct SubSequence {
    pub start: usize,
    pub end: usize,
    pub seq: Sequence,
}

pub type CpgIsland = SubSequence;
pub type Orf = SubSequence;

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
    pub fn to_string(seq: &[Nucleotide]) -> String {
        seq.iter().map(|n| n.to_char()).collect()
    }

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
        Self::gc_percent(&self.0)
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

    /// Find ORFs
    pub fn find_orfs(&self, min_len: usize) -> Result<Vec<Orf>, Box<dyn std::error::Error>> {
        Self::find_orfs_in_seq(&self.0, min_len)
    }

    fn find_orfs_in_seq(
        seq: &[Nucleotide],
        min_len: usize,
    ) -> Result<Vec<Orf>, Box<dyn std::error::Error>> {
        /// ATG
        const START_CODON: [Nucleotide; 3] = [
            Nucleotide::Adenine,
            Nucleotide::Thymine,
            Nucleotide::Guanine,
        ];

        // TAA, TAG, TGA
        const STOP_CODONS: [[Nucleotide; 3]; 3] = [
            [
                Nucleotide::Thymine,
                Nucleotide::Adenine,
                Nucleotide::Adenine,
            ],
            [
                Nucleotide::Thymine,
                Nucleotide::Adenine,
                Nucleotide::Guanine,
            ],
            [
                Nucleotide::Thymine,
                Nucleotide::Guanine,
                Nucleotide::Adenine,
            ],
        ];

        let mut orfs = Vec::new();

        for frame_id in 0..3 {
            // Find start-codon
            for i in (frame_id..seq.len()).step_by(3) {
                let Some(start_codon) = seq.get(i..(i + 3)) else {
                    continue;
                };

                if start_codon != START_CODON {
                    continue;
                }

                // Find stop-codon
                for j in (i..seq.len()).step_by(3) {
                    let Some(stop_codon) = seq.get(j..(j + 3)) else {
                        continue;
                    };

                    if !STOP_CODONS.contains(stop_codon.try_into()?) {
                        continue;
                    }

                    let Some(orf_candidate) = seq.get(i..(j + 3)) else {
                        continue;
                    };

                    if orf_candidate.len() < min_len {
                        continue;
                    }

                    let mut orf_seq = Sequence::default();
                    orf_seq.0.extend_from_slice(orf_candidate);

                    orfs.push(Orf {
                        start: i,
                        end: j + 3,
                        seq: orf_seq,
                    });

                    break;
                }
            }
        }

        Ok(orfs)
    }

    /// Find CpG islands
    pub fn find_cpg_islands(
        &self,
        min_len: usize,
        min_gc: f64,
        min_oe: f64,
        step: usize,
    ) -> Vec<CpgIsland> {
        let mut islands = Vec::<(usize, usize)>::new();
        let seq_len = self.length();

        let mut i = 0;
        while i <= (seq_len - min_len) {
            let window = &self.0[i..(i + min_len)];
            let gc = Self::gc_percent(window);
            let observed_expected = Self::observed_expected_cpg(window);

            if gc >= min_gc && observed_expected >= min_oe {
                let mut start = i;
                let mut end = i + min_len;

                // Expand to the left side
                while start > 0 {
                    let prev = &self.0[(start - 1)..(start + min_len - 1)];

                    if Self::gc_percent(prev) >= min_gc
                        && Self::observed_expected_cpg(prev) >= min_oe
                    {
                        start -= 1;
                    } else {
                        break;
                    }
                }

                // Expand to the right side
                while end < seq_len {
                    let next = &self.0[(end - min_len + 1)..(end + 1)];

                    if Self::gc_percent(next) >= min_gc
                        && Self::observed_expected_cpg(next) >= min_oe
                    {
                        end += 1;
                    } else {
                        break;
                    }
                }

                islands.push((start, end));
                i = end;
            } else {
                i += step;
            }
        }

        // Sort islands
        islands.sort_by_key(|island| island.0);

        // Merge intersected CpG islands
        let merged = Self::merge_cpg_islands(islands);

        // Create slices
        merged
            .into_iter()
            .map(|(start, end)| {
                let mut seq = Sequence::default();
                seq.0.extend_from_slice(&self.0[start..end]);

                CpgIsland { start, end, seq }
            })
            .collect()
    }

    pub fn cpg_islands_metrics(
        &self,
        window: usize,
        step: usize,
    ) -> (Vec<usize>, Vec<f64>, Vec<f64>) {
        let mut positions = Vec::new();
        let mut gc_values = Vec::new();
        let mut oe_values = Vec::new();

        let mut start = 0;
        while start + window <= self.length() {
            let slice = &self.0[start..start + window];

            let gc = Self::gc_percent(slice);
            let oe = Self::observed_expected_cpg(slice);

            positions.push(start + window / 2);
            gc_values.push(gc);
            oe_values.push(oe);

            start += step;
        }

        (positions, gc_values, oe_values)
    }

    fn observed_expected_cpg(seq: &[Nucleotide]) -> f64 {
        let observed_gc_count = Self::cg_pair_count(seq);

        let c_count = Self::nucleotide_count(seq, Nucleotide::Cytosine);
        let g_count = Self::nucleotide_count(seq, Nucleotide::Guanine);
        let length = seq.len();

        if c_count == 0 || g_count == 0 || length == 0 {
            return 0.0;
        }

        let expected_gc_count = c_count * g_count / seq.len();
        (observed_gc_count as f64) / (expected_gc_count as f64)
    }

    fn merge_cpg_islands(islands: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
        if islands.is_empty() {
            return islands;
        }

        let mut merged = Vec::with_capacity(islands.len());

        let mut current = islands[0];
        for next in islands.iter().skip(1) {
            // Check if two islands intersected
            if next.0 <= current.1 {
                // merge intersected islands
                current.1 = current.1.max(next.1);
            } else {
                // store island and move to the next island
                merged.push(current);
                current = *next;
            }
        }

        // Store last island
        merged.push(current);

        merged.shrink_to_fit();
        merged
    }

    fn gc_percent(seq: &[Nucleotide]) -> f64 {
        let gc = seq
            .iter()
            .filter(|n| matches!(n, Nucleotide::Guanine | Nucleotide::Cytosine)) // yeild only G and C
            .count();

        (gc as f64) / (seq.len() as f64) * 100.0
    }

    fn nucleotide_count(seq: &[Nucleotide], expected: Nucleotide) -> usize {
        seq.iter().filter(|n| **n == expected).count()
    }

    fn cg_pair_count(seq: &[Nucleotide]) -> usize {
        Self::pair_count(seq, [Nucleotide::Cytosine, Nucleotide::Guanine])
    }

    fn pair_count(seq: &[Nucleotide], pair: [Nucleotide; 2]) -> usize {
        seq.windows(2)
            .filter(|w| w[0] == pair[0] && w[1] == pair[1])
            .count()
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

    #[test]
    fn find_cpg_islands_test() {
        // Given
        const WINDOW_LEN: usize = 150;
        const WINDOW_STEP: usize = 10;
        const MIN_GC: f64 = 50.0;
        const MIN_OE: f64 = 0.6;

        let mut raw_data = String::new();
        raw_data.push_str(&"A".repeat(100));
        raw_data.push_str(&"CGCGCGCGCGCGCGCGCG".repeat(20));
        raw_data.push_str(&"T".repeat(100));

        let seq = Sequence::from_str(&raw_data).unwrap();

        // When
        let islands = seq.find_cpg_islands(WINDOW_LEN, MIN_GC, MIN_OE, WINDOW_STEP);

        assert_eq!(islands.len(), 1);
        assert_eq!(islands[0].start, 25);
        assert_eq!(islands[0].end, 535);
    }
}
