use std::{collections::HashMap, io::Read, path::Path, str::FromStr};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GffFeatureType {
    // Areas
    Gene,
    PseudoGene,
    Transcript, // mRNA
    IncRna,     // miRNA, rRNA, tRNA
    MobileGeneticElement,

    // Structures
    Exon,
    Intron,
    Cds,
    CdsExon,
    FivePrime,  // 5' UTR
    ThreePrime, // 3' UTR
    StemLoop,
    MatureProteinRegionOfCds,

    // Signals
    StartCodon,
    StopCodon,
    SelenoCystein,

    // Chromosomes
    Chromosome,
    RepeatRegion,
    Contig, // scaffold

    Unknown(String),
}

pub struct GffFeature {
    pub seq_id: String,
    pub source: String,
    pub feature_type: GffFeatureType,
    pub start: usize,
    pub end: usize,
    pub score: String,
    pub strand: String,
    pub phase: String,
    pub attributes: HashMap<String, String>,
}

pub struct GffFile {
    pub features: Vec<GffFeature>,
}

impl GffFile {
    pub fn open(p: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        let mut input = std::fs::File::open(p)?;
        let mut raw_data = String::new();
        input.read_to_string(&mut raw_data)?;

        Self::from_text(&raw_data)
    }

    pub fn from_text(s: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let mut features = Vec::new();

        for line in s.lines() {
            let trimmed = line.trim();

            if line.starts_with("#") || trimmed.is_empty() {
                continue;
            }

            let parts: Vec<&str> = trimmed.split('\t').collect();

            let Some(seq_id) = parts.first().map(|s| s.to_string()) else {
                return Err("Can't get seq_id".into());
            };

            let Some(source) = parts.get(1).map(|s| s.to_string()) else {
                return Err("Can't get source".into());
            };

            let Some(feature_type) = parts.get(2).map(|s| s.to_string()) else {
                return Err("Can't get feature_type".into());
            };

            let feature_type = GffFeatureType::from_str(&feature_type)
                .unwrap_or(GffFeatureType::Unknown(feature_type));

            let Some(Ok(start)) = parts.get(3).map(|s| s.parse::<usize>()) else {
                return Err("Can't get start".into());
            };

            let Some(Ok(end)) = parts.get(4).map(|s| s.parse::<usize>()) else {
                return Err("Can't get end".into());
            };

            let Some(score) = parts.get(5).map(|s| s.to_string()) else {
                return Err("Can't get score".into());
            };

            let Some(strand) = parts.get(6).map(|s| s.to_string()) else {
                return Err("Can't get strand".into());
            };

            let Some(phase) = parts.get(7).map(|s| s.to_string()) else {
                return Err("Can't get phase".into());
            };

            let Some(attributes) = parts.get(8) else {
                return Err("Can't get attributes".into());
            };

            let attributes = attributes
                .split(';')
                .filter_map(|attribute| {
                    if attribute.contains('=') {
                        let splitted_attribute = attribute.splitn(2, '=').collect::<Vec<&str>>();
                        let key = splitted_attribute[0].to_string();
                        let value = splitted_attribute[1].to_string();

                        Some((key, value))
                    } else {
                        None
                    }
                })
                .collect::<HashMap<String, String>>();

            features.push(GffFeature {
                seq_id,
                source,
                feature_type,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
            });
        }

        Ok(Self { features })
    }
}

impl FromStr for GffFeatureType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let lowercase_s = s.to_lowercase();

        match (s, &lowercase_s[..]) {
            (_, "gene") => Ok(Self::Gene),
            (_, "pseudogene") => Ok(Self::PseudoGene),
            (_, "transcript") | ("mRNA", _) => Ok(Self::Transcript),
            ("Inc_RNA", _) | ("miRNA", _) | ("rRNA", _) | ("tRNA", _) => Ok(Self::IncRna),
            (_, "mobile_genetic_element") => Ok(Self::MobileGeneticElement),

            (_, "exon") => Ok(Self::Exon),
            (_, "intron") => Ok(Self::Intron),
            (_, "cds") => Ok(Self::Cds),
            (_, "cds_exon") => Ok(Self::CdsExon),
            (_, "five_prime_utr") => Ok(Self::FivePrime),
            (_, "three_prime_utr") => Ok(Self::ThreePrime),
            (_, "stem_loop") | (_, "stem-loop") => Ok(Self::StemLoop),
            (_, "mature_protein_region_of_cds") => Ok(Self::MatureProteinRegionOfCds),

            (_, "start_codon") => Ok(Self::StartCodon),
            (_, "stop_codon") => Ok(Self::StopCodon),
            (_, "selenocystein") => Ok(Self::SelenoCystein),

            (_, "chromosome") => Ok(Self::Chromosome),
            (_, "repeat_region") => Ok(Self::RepeatRegion),
            (_, "contig") | (_, "scaffold") => Ok(Self::Contig),
            _ => Err(format!("Unknown feature: {s}")),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_one_feature_test() {
        // Given
        const SEQ_ID: &str = "SEQ_ID_0.0";
        const SOURCE: &str = "RefSeq";
        const FEATURE_TYPE: &str = "gene";
        const START: usize = 1;
        const END: usize = 100;
        const SCORE: &str = ".";
        const STRAND: &str = "+";
        const PHASE: &str = ".";

        const ATTRIBUTES: [(&str, &str); 5] = [
            ("ID", "SEQ_ID_0.0:1..100"),
            ("Dbxref", "taxon:2222"),
            ("collection-date", "June-2026"),
            ("country", "A"),
            ("gb-acronym", "TEST-COV"),
        ];

        let raw_data = format!(
            "#SomeHeader\
          \n{SEQ_ID}\t{SOURCE}\t{FEATURE_TYPE}\t{START}\t{END}\t{SCORE}\t{STRAND}\t{PHASE}\
          \t{}={};{}={};{}={};{}={};{}={};",
            ATTRIBUTES[0].0,
            ATTRIBUTES[0].1,
            ATTRIBUTES[1].0,
            ATTRIBUTES[1].1,
            ATTRIBUTES[2].0,
            ATTRIBUTES[2].1,
            ATTRIBUTES[3].0,
            ATTRIBUTES[3].1,
            ATTRIBUTES[4].0,
            ATTRIBUTES[4].1,
        );

        // When
        let gff = GffFile::from_text(&raw_data).unwrap();

        // Then
        assert_eq!(gff.features.len(), 1);
        assert_eq!(gff.features[0].seq_id, SEQ_ID);
        assert_eq!(gff.features[0].source, SOURCE);
        assert_eq!(gff.features[0].feature_type, GffFeatureType::Gene);
        assert_eq!(gff.features[0].start, START);
        assert_eq!(gff.features[0].end, END);
        assert_eq!(gff.features[0].score, SCORE);
        assert_eq!(gff.features[0].strand, STRAND);
        assert_eq!(gff.features[0].phase, PHASE);

        assert_eq!(
            gff.features[0]
                .attributes
                .get(ATTRIBUTES[0].0)
                .map(|x| x.to_string()),
            Some(ATTRIBUTES[0].1).map(|x| x.to_string())
        );

        for (key, value) in ATTRIBUTES {
            assert_eq!(
                gff.features[0].attributes.get(key).map(|s| s.as_str()),
                Some(value)
            );
        }
    }
}
