use std::fmt::Display;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Nucleotide {
    Adenine,
    Guanine,
    Cytosine,
    Thymine,
    Uracil,
}

impl Display for Nucleotide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

impl TryFrom<char> for Nucleotide {
    type Error = String;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value.to_ascii_lowercase() {
            'a' => Ok(Self::Adenine),
            'g' => Ok(Self::Guanine),
            'c' => Ok(Self::Cytosine),
            't' => Ok(Self::Thymine),
            'u' => Ok(Self::Uracil),
            _ => Err(format!("\"{value}\" is not a nucleotide")),
        }
    }
}

impl Nucleotide {
    pub fn to_char(&self) -> char {
        match self {
            Self::Adenine => 'A',
            Self::Guanine => 'G',
            Self::Cytosine => 'C',
            Self::Thymine => 'T',
            Self::Uracil => 'U',
        }
    }
}
