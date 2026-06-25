use std::{fmt::Display, str::FromStr};

use crate::{codons::CodonSequence, nucleotide::Nucleotide};

#[derive(Default, Clone)]
pub struct Sequence(pub Vec<Nucleotide>);

#[derive(Default, Clone)]
pub struct SubSequence {
    pub start: usize,
    pub end: usize,
    pub seq: Sequence,
}

impl SubSequence {
    pub fn new_in(seq: &[Nucleotide], start: usize, end: usize) -> Option<Self> {
        let mut sub_seq = Sequence::default();

        sub_seq.0.extend_from_slice(seq.get(start..end)?);

        Some(Self {
            start,
            end,
            seq: sub_seq,
        })
    }
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
        Self::is_dna_seq(&self.0)
    }

    pub fn is_rna(&self) -> bool {
        !self.is_dna()
    }

    fn is_dna_seq(seq: &[Nucleotide]) -> bool {
        Self::get_seq_kind(seq) == SequenceKind::Dna
    }

    pub fn get_kind(&self) -> SequenceKind {
        Self::get_seq_kind(&self.0)
    }

    fn get_seq_kind(seq: &[Nucleotide]) -> SequenceKind {
        if seq.contains(&Nucleotide::U) {
            SequenceKind::Rna
        } else {
            SequenceKind::Dna
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
            if *n == Nucleotide::T {
                *n = Nucleotide::U
            }
        });

        Ok(transcribed)
    }

    /// Applicable to: RNA sequence kind
    pub fn translate(&self) -> Result<CodonSequence, String> {
        if !self.is_rna() {
            return Err(format!(
                "Only DNA can be transcribed. Sequence kind is {}",
                self.get_kind()
            ));
        }

        let translated = CodonSequence::from(self.0.as_ref());

        Ok(translated)
    }

    pub fn get_hamming_distance(&self, other: &Self) -> Result<usize, String> {
        Self::get_hamming_distance_in_seq(&self.0, &other.0)
    }

    fn get_hamming_distance_in_seq(
        first: &[Nucleotide],
        second: &[Nucleotide],
    ) -> Result<usize, String> {
        let first_kind = Self::get_seq_kind(first);
        let second_kind = Self::get_seq_kind(second);

        if first_kind != second_kind {
            return Err(format!(
                "Wrong comparison. Sequence kind should be same, actually: {first_kind} and {second_kind}"
            ));
        }

        if first.len() != second.len() {
            return Err(format!(
                "Wrong comparison. Sequence length should be equal, actually: {} and {}",
                first.len(),
                second.len()
            ));
        }

        Ok(first
            .iter()
            .zip(second.iter())
            .filter(|(a, b)| **a != **b)
            .count())
    }

    pub fn splice(&self, introns: &[Sequence]) -> Self {
        let introns = introns.iter().map(|a| a.0.as_ref()).collect::<Vec<_>>();

        Self(Self::splice_seq(&self.0, &introns))
    }

    fn splice_seq(seq: &[Nucleotide], introns: &[&[Nucleotide]]) -> Vec<Nucleotide> {
        let mut new_seq = Vec::from(seq);

        for intron in introns {
            new_seq = Self::splice_in_seq_one(&new_seq, intron);
        }

        new_seq
    }

    fn splice_in_seq_one(seq: &[Nucleotide], intron: &[Nucleotide]) -> Vec<Nucleotide> {
        let mut seq = Vec::from(seq);

        if let Some(idx) = seq
            .windows(intron.len())
            .position(|window| window == intron)
        {
            // Удаляем найденную подпоследовательность
            seq.drain(idx..idx + intron.len());
        }

        seq
    }

    /// Find intersections of two given sequences
    pub fn find_subsequence_intersections(
        &self,
        first: &[SubSequence],
        second: &[SubSequence],
    ) -> Vec<(SubSequence, usize, usize)> {
        Self::find_subsequence_intersections_in(&self.0, first, second)
    }

    fn find_subsequence_intersections_in(
        seq: &[Nucleotide],
        first: &[SubSequence],
        second: &[SubSequence],
    ) -> Vec<(SubSequence, usize, usize)> {
        let mut intersections = Vec::new();

        let mut i = 0;
        let mut j = 0;

        while i < first.len() && j < second.len() {
            let a = &first[i];
            let b = &second[j];

            if a.start <= b.end && b.start <= a.end {
                let start = std::cmp::max(a.start, b.start);
                let end = std::cmp::max(a.end, b.end);

                let intersection = SubSequence::new_in(seq, start, end).unwrap();

                intersections.push((intersection, i, j));
            }

            if a.end < b.end {
                i += 1;
            } else {
                j += 1;
            }
        }

        intersections
    }

    /// Find ORFs
    pub fn find_orfs(&self, min_len: usize) -> Result<Vec<Orf>, Box<dyn std::error::Error>> {
        Self::find_orfs_in_seq(&self.0, min_len)
    }

    fn find_orfs_in_seq(
        seq: &[Nucleotide],
        min_len: usize,
    ) -> Result<Vec<Orf>, Box<dyn std::error::Error>> {
        const DNA_START_CODON: [Nucleotide; 3] = [Nucleotide::A, Nucleotide::T, Nucleotide::G];
        const RNA_START_CODON: [Nucleotide; 3] = [Nucleotide::A, Nucleotide::U, Nucleotide::G];

        const DNA_STOP_CODONS: [[Nucleotide; 3]; 3] = [
            [Nucleotide::T, Nucleotide::A, Nucleotide::A],
            [Nucleotide::T, Nucleotide::A, Nucleotide::G],
            [Nucleotide::T, Nucleotide::G, Nucleotide::A],
        ];

        const RNA_STOP_CODONS: [[Nucleotide; 3]; 3] = [
            [Nucleotide::U, Nucleotide::A, Nucleotide::A],
            [Nucleotide::U, Nucleotide::A, Nucleotide::G],
            [Nucleotide::U, Nucleotide::G, Nucleotide::A],
        ];

        let (start_codon, stop_codons) = if Self::is_dna_seq(seq) {
            (DNA_START_CODON, DNA_STOP_CODONS)
        } else {
            (RNA_START_CODON, RNA_STOP_CODONS)
        };

        let mut orfs = Vec::new();

        for start in (0..seq.len()).step_by(3) {
            let Some(codon) = seq.get(start..(start + 3)) else {
                continue;
            };

            if codon != start_codon {
                continue;
            }

            // Find stop-codon
            for end in ((start + 3)..seq.len()).step_by(3) {
                let Some(codon) = seq.get(end..(end + 3)) else {
                    continue;
                };

                if !stop_codons.contains(codon.try_into()?) {
                    continue;
                }

                let Some(candidate) = SubSequence::new_in(seq, start, end + 3) else {
                    continue;
                };

                if candidate.seq.length() >= min_len {
                    orfs.push(candidate);
                    break;
                }
            }
        }

        Ok(orfs)
    }

    pub fn best_orfs(orfs: Vec<SubSequence>) -> Vec<SubSequence> {
        if orfs.is_empty() {
            return orfs;
        }

        let mut merged = Vec::with_capacity(orfs.len());

        let mut current = orfs[0].clone();
        for next in orfs.iter().skip(1) {
            // Check if two orfs intersected
            if !(next.start >= current.start && next.end <= current.end) {
                // store ORF and move to the next
                merged.push(current.clone());
                current = next.clone();
            }
        }

        // Store last ORF
        merged.push(current.clone());

        merged.shrink_to_fit();
        merged
    }

    pub fn find_motifs(&self, motif: &[Nucleotide]) -> Vec<usize> {
        Self::find_motifs_in_seq(&self.0, motif)
    }

    fn find_motifs_in_seq(seq: &[Nucleotide], motif: &[Nucleotide]) -> Vec<usize> {
        seq.windows(motif.len())
            .enumerate()
            .filter_map(|(pos, sub_seq)| if sub_seq == motif { Some(pos) } else { None })
            .collect()
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

    pub fn get_window_gc_oe_metrics(
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

        let c_count = Self::nucleotide_count(seq, Nucleotide::C);
        let g_count = Self::nucleotide_count(seq, Nucleotide::G);
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
            .filter(|n| matches!(n, Nucleotide::G | Nucleotide::C)) // yeild only G and C
            .count();

        (gc as f64) / (seq.len() as f64) * 100.0
    }

    fn nucleotide_count(seq: &[Nucleotide], expected: Nucleotide) -> usize {
        seq.iter().filter(|n| **n == expected).count()
    }

    fn cg_pair_count(seq: &[Nucleotide]) -> usize {
        Self::pair_count(seq, [Nucleotide::C, Nucleotide::G])
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

    #[test]
    fn find_motif_test() {
        // Given
        const RAW_DATA: &str = "GATATATGCATATACTT";
        const MOTIF_STR: &str = "ATAT";

        let seq = Sequence::from_str(RAW_DATA).unwrap();
        let motif = Sequence::from_str(MOTIF_STR).unwrap();

        // When
        let positions = seq.find_motifs(&motif.0);
        assert_eq!(positions.len(), 3);
        assert_eq!(positions[0], 1);
        assert_eq!(positions[1], 3);
        assert_eq!(positions[2], 9);
    }

    #[test]
    fn get_hamming_distance_test() {
        // Given
        const FIRST: &str = "GAGCCTACTAACGGGAT";
        const SECOND: &str = "CATCGTAATGACGGCCT";

        let first = Sequence::from_str(FIRST).unwrap();
        let second = Sequence::from_str(SECOND).unwrap();

        // When
        let res = first.get_hamming_distance(&second).unwrap();

        // Then
        const EXPECTED: usize = 7;

        assert_eq!(res, EXPECTED);
    }

    #[test]
    fn splice_test() {
        // Given
        const RAW_DATA: &str = "ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG";
        const INTRON_1: &str = "ATCGGTCGAA";
        const INTRON_2: &str = "ATCGGTCGAGCGTGT";

        let seq = Sequence::from_str(RAW_DATA).unwrap();
        let introns = vec![
            Sequence::from_str(INTRON_1).unwrap(),
            Sequence::from_str(INTRON_2).unwrap(),
        ];

        // When
        let spliced = seq.splice(&introns);
        let transribed = spliced.transcribe().unwrap();
        let translated = transribed.translate().unwrap();

        // Then
        const EXPECTED: &str = "MVYIADKQHVASREAYGHMFKVCA";
        assert_eq!(translated.0.len(), 24);
        assert_eq!(CodonSequence::to_string(&translated.0), EXPECTED);
    }
}
