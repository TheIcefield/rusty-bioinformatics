use std::fmt::Display;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Nucleotide {
    A, // Adenine
    G, // Guanine
    C, // Cytosine
    T, // Thymine
    U, // Uracil
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
            'a' => Ok(Self::A),
            'g' => Ok(Self::G),
            'c' => Ok(Self::C),
            't' => Ok(Self::T),
            'u' => Ok(Self::U),
            _ => Err(format!("\"{value}\" is not a nucleotide")),
        }
    }
}

impl Nucleotide {
    pub fn to_char(&self) -> char {
        match self {
            Self::A => 'A',
            Self::G => 'G',
            Self::C => 'C',
            Self::T => 'T',
            Self::U => 'U',
        }
    }
}
