use std::str::FromStr;

use nucleotides::sequence::Sequence;

#[derive(Default)]
pub struct FnaRecord {
    pub header: String,
    pub content: Sequence,
}

pub struct FnaFile {
    pub records: Vec<FnaRecord>,
}

impl FnaFile {
    pub fn read(s: &str) -> Result<Self, String> {
        let mut records = Vec::<FnaRecord>::new();
        let mut current_record = FnaRecord::default();

        // Read by lines
        for line in s.split('\n').map(|line| line.trim()) {
            if line.starts_with('>') {
                // If already has something - store it
                if !current_record.header.is_empty() {
                    // Store current record
                    records.push(current_record);

                    // Reset current record
                    current_record = FnaRecord::default();
                }

                current_record
                    .header
                    .push_str(line.get(1..).ok_or("failed to get slice \"1..\"")?);
            } else {
                let mut nucleotides = Sequence::from_str(line)?;

                current_record.content.0.append(&mut nucleotides.0);
            }
        }

        // Store the rest
        if !current_record.header.is_empty() {
            records.push(current_record);
        }

        Ok(FnaFile { records })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_one_record_test() {
        // Given
        let raw_data = ">some_header\
                      \nAATTGGCCTAGCAGAC\
                      \nACGTCGTAGTACGTGC";

        // When
        let fna = FnaFile::read(raw_data).unwrap();

        // Then
        assert_eq!(fna.records.len(), 1);
        assert_eq!(fna.records[0].header, "some_header");
        assert_eq!(fna.records[0].content.0.len(), 32);
        assert_eq!(fna.records[0].content.get_gc_content(), 53.125);
    }
}
