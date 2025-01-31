# Introduction to the Bioconductor Project

**Bioconductor** is an open-source software ecosystem for **bioinformatics** and **genomic analysis** in R.  
It provides specialized tools for working with DNA sequences, gene expression data, and biological annotations.

In this guide, we explore **BSgenome** and **Biostrings**, Bioconductor packages that provide whole-genome sequences and efficient biological string manipulation in R.

---

## **Biological Background (For Non-Biologists)**

### **What is a Genome?**
A **genome** is the complete set of DNA in an organism, encoded in a four-letter **alphabet**:
- **A (Adenine)**  
- **T (Thymine)**  
- **G (Guanine)**  
- **C (Cytosine)**  

### **Genome Structure**

A genome consists of:
- **Chromosomes** – Large DNA structures containing genes.
- **Genes** – Segments of DNA that code for proteins (some are non-coding but still have functions).
- **Proteins** – Molecules responsible for cellular functions, produced through:
  1. **Transcription**: DNA → RNA
  2. **Translation**: RNA → Protein

![DNA-Split](images/DNA-Transcription-Translation.png)

 --- 

## **Installation and Setup**

To use Bioconductor, first install the `BiocManager` package (recommended by R for managing Bioconductor packages).  
Then, install the **GenomicRanges**, **BSgenome**, and **Biostrings** packages, which are required for working with genome sequences and biological strings.

```r
# Install BiocManager and required packages
install.packages("BiocManager")
BiocManager::install("GenomicRanges") 
BiocManager::install("BSgenome")
BiocManager::install("Biostrings")

# Load the packages and check their versions
library(BiocManager)
library(GenomicRanges)
library(BSgenome)
library(Biostrings)

# Check R session details (useful for reproducibility)
sessionInfo()
```

---

## **Understanding the `BSgenome` Class**

`BSgenome` is an **S4 class**, meaning it has a structured object-oriented representation in R.

```r
# Suppose we have an object `a_genome` from class BSgenome
class(a_genome)       # Returns: "BSgenome"
is(a_genome)          # Lists inherited classes
isS4(a_genome)        # Returns: TRUE
slotNames(a_genome)   # Returns: c("organism", "provider", "seqinfo", ...)
```

### **Key Accessor Functions**
These functions help extract information from a `BSgenome` object.

```r
show(a_genome)       # Prints an overview of the genome object
organism(a_genome)   # Returns: "Saccharomyces cerevisiae"
provider(a_genome)   # Returns: "UCSC"
seqinfo(a_genome)    # Displays sequence metadata (chromosome lengths, circularity, etc.)
```

---

## **Exploring Available Genomes**
Bioconductor provides **pre-built genomes** for various species.  
You can check available genomes using:

```r
available.genomes()   # Returns: c("BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Mmusculus.UCSC.mm10", ...)
```

### **Example: Working with the Yeast Genome**
Let's load a yeast (`Saccharomyces cerevisiae`) genome from Bioconductor.

```r
# Load the yeast genome package
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Create a yeast genome object
yeast <- BSgenome.Scerevisiae.UCSC.sacCer3

# Inspect the number of chromosomes
length(yeast)    # Returns: 17

# Get the chromosome names (names of stored sequences)
names(yeast)     # Returns: c("chrI", "chrII", "chrIII", ..., "chrM")

# Retrieve chromosome lengths (number of DNA base pairs)
seqlengths(yeast)    # Returns named vector: c(chrI = 230218,
                     #          chrII = 813184, chrIII = 316620, ..., chrM = 85779)
```

### **Extracting DNA Sequences**
We can extract genomic sequences using `getSeq()`.

```r
# Retrieve the sequence of chromosome M
getSeq(yeast, "chrM")    # Returns: 85779-letter DNAString object
                         # seq: TTCATAATTAATTTTTTATATATATATTATATTATA...TACAGAAATATGCTTAATTATAATATAATATCCATA

# Retrieve the first 10 base pairs of the genome
getSeq(yeast, end = 10)  # Returns: DNAStringSet object of length 17:
                         # width seq        names               
                         # [1]  10 CCACACCACA chrI
                         # [2]  10 AAATAGCCCT chrII
                         # [3]  10 CCCACACACC chrIII
                         # [4]  10 ACACCACACC chrIV
                         # [5]  10 CGTCTCCTCC chrV
                         # ...
```

---

## **Introduction to Biostrings**

The **Biostrings** package implements algorithms for **fast manipulation of large biological sequences**. It is widely used, with more than 200 Bioconductor packages depending on it.

### **Biological String Containers**
Biological sequences in Biostrings are stored in **memory-efficient containers**, which allow efficient subsetting and pattern matching. The containers include:
- `BString` – Stores a generic big string.
- `DNAString` – Stores DNA sequences.
- `RNAString` – Stores RNA sequences.
- `AAString` – Stores amino acid sequences.

To store **multiple sequences**, we use **StringSet containers**:
- `BStringSet`
- `DNAStringSet`
- `RNAStringSet`
- `AAStringSet`

### **Checking Class Structures**
We can inspect how these classes are related using `showClass()`:

```r
showClass("XString")   # Returns:
                      # Virtual Class "XString" [package "Biostrings"]
                      # Slots:
                      # Name:      shared      offset      length   elementMetadata  metadata
                      # Class:  SharedRaw    integer     integer  DataFrame_OR_NULL  list
                      # Extends: "XRaw", "XVector", "Vector", "Annotated"
                      # Known Subclasses: "BString", "DNAString", "RNAString", "AAString"
```

### **Biostring Alphabets**
Biological sequences follow **predefined alphabets**:
- `DNA_BASES` – A, C, G, T
- `RNA_BASES` – A, C, G, U (U replaces T in RNA)
- `AA_STANDARD` – 20 standard amino acid letters

Additional predefined alphabets include:
- `DNA_ALPHABET` and `RNA_ALPHABET` – Includes **IUPAC_CODE_MAP** (allows ambiguity codes like N for any base).
- `AA_ALPHABET` – Based on **AMINO_ACID_CODE**.

For more details on IUPAC DNA codes, visit: [UCSC Genome IUPAC Codes](http://genome.ucsc.edu/goldenPath/help/iupac.html).

---

## **Transcription and Translation**

### **Transcription: DNA to RNA**
In transcription, DNA is converted into RNA by replacing **T** with **U**.

```r
# Create a DNA string
dna_seq <- DNAString("ATGATCTCGTAA")
dna_seq              # Returns: 12-letter DNAString object
                     # seq: ATGATCTCGTAA

# Transcription: Convert DNA to RNA
rna_seq <- RNAString(dna_seq)
rna_seq              # Returns: 12-letter RNAString object
                     # seq: AUGAUCUCGUAA
```

### **Translation: RNA to Amino Acids**
To translate RNA into **Amino Acids**, we use `translate()`.

```r
# Translate RNA to Amino Acid sequence
aa_seq <- translate(rna_seq)
aa_seq               # Returns: 4-letter AAString object
                     # seq: MIS*

# Shortcut: DNA to Amino Acids
translate(dna_seq)
```

---

## **Sequence Manipulation Functions**

### **Understanding Single Sequences vs. Sets**

Biostrings provides two main ways to work with sequences:
- **Single Sequences** (`XString` objects): Hold one sequence of a predefined alphabet
- **Sets** (`StringSet` objects): Store multiple sequences that can have varying lengths

The choice between single sequences and sets affects which functions are available and how sequence lengths are reported.

### **Converting Between Single Sequences and Sets**

#### From Set to Single Sequence
```r
# Read a sequence file into a StringSet
zikaVirus <- readDNAStringSet("data/zika.fa")
length(zikaVirus)    # Returns: 1 (number of sequences in set)
width(zikaVirus)     # Returns: 10794 (number of bases in sequence)

# Convert StringSet to single sequence
zikaVirus_seq <- unlist(zikaVirus)    # Creates a DNAString
length(zikaVirus_seq)  # Returns: 10794 (length of single sequence)
```

Note the difference between `length()` and `width()`:
- For sets: `length()` gives number of sequences, `width()` gives sequence lengths
- For single sequences: only `length()` is applicable, giving sequence length

#### From Single Sequence to Set
```r
# Create a set with multiple subsequences
zikaSet <- DNAStringSet(zikaVirus_seq, 
                       start = c(1, 101, 201),    # Starting positions
                       end = c(100, 200, 300))    # Ending positions
# Creates a set with three 100-base sequences
length(zikaSet)     # Returns: 3 (number of sequences)
width(zikaSet)      # Returns: c(100, 100, 100) (length of each sequence)
```

### **Sequence Manipulation Operations**

#### DNA Strand Operations
DNA naturally exists as two complementary strands (A pairs with T, G pairs with C). Biostrings efficiently stores just one strand and can computationally derive the other when needed:

```r
# Create a DNA sequence
dna_seq <- DNAString("ATGATCTCGTAA")

# Get complementary sequence (A↔T, G↔C)
complement(dna_seq)    # Returns: TACTAGAGCATT

# Reverse a sequence
reverse(dna_seq)       # Returns: AATGCTCTAGTA

# Reverse complement (common operation in molecular biology)
reverseComplement(dna_seq)    # Returns: TTACGAGATCAT
```

#### Working with Sets
```r
# Create a short set for demonstration
shortSet <- DNAStringSet(c(
  seq1 = "AGCTCGTAGCTAGCTAG",
  seq2 = "CGATCGATCGATCGATC"
))

# Reorder sequences in set (top to bottom)
rev(shortSet)    # Returns set with seq2 first, then seq1

# Reverse each sequence (right to left)
reverse(shortSet)    # Returns set with reversed sequences

# Get complement of all sequences in set
complement(shortSet)  # Returns complementary sequences for entire set
```

### **Function Summary Table**

| Function | Description | Single Sequence | Set |
|----------|-------------|-----------------|-----|
| `length()` | Returns sequence length or number of sequences in set | Yes (sequence length) | Yes (number of sequences) |
| `width()` | Number of characters per sequence | No | Yes |
| `unlist()` | Collates set elements into single sequence | No | Yes |
| `complement()` | Returns paired DNA strand | Yes | Yes |
| `rev()` | Reorders sequences (bottom to top) for sets; same as reverse() for single | Yes | Yes |
| `reverse()` | Changes sequence order (right to left) | Yes | Yes |
| `reverseComplement()` | Combines reverse and complement operations efficiently | Yes | Yes |

### **Performance and Best Practices**

1. **Use `reverseComplement()` for Efficiency**
   - When you need both operations, use `reverseComplement()` instead of separate `reverse()` and `complement()`
   - This is particularly important for large sequences as it's optimized for memory and speed

2. **Choosing Between Single Sequences and Sets**
   - Use single sequences when working with one continuous sequence
   - Use sets when you need to:
     - Store multiple related sequences
     - Break a sequence into chunks for parallel processing
     - Compare multiple sequence variants

3. **Memory Considerations**
   - StringSets can hold sequences of different lengths
   - For very large sequences, work with subsequences using `subseq()` to manage memory
   - Use `width()` to check sequence lengths before operations that might be memory-intensive

### **Practical Examples**

```r
# Extract and analyze a subsequence
short_seq <- subseq(zikaVirus_seq, end = 21)

# Transcribe DNA to RNA
rna_seq <- RNAString(short_seq)

# Translate RNA to amino acids
aa_seq <- translate(rna_seq)

# Create a set of overlapping sequences
overlap_set <- DNAStringSet(zikaVirus_seq, 
                          start = seq(1, 100, by = 10),
                          width = 20)
```

---

## **Conclusion**
This document introduces Bioconductor, **BSgenome**, and **Biostrings** for working with genome sequences. 
For further exploration, check out:

🔗 [Bioconductor Website](https://bioconductor.org/)