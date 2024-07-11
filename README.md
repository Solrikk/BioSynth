![Logo](https://github.com/Solrikk/BioSynth/blob/main/images/results/uHfj7EgN6S.png)

<div align="center">
  <h3>
    <a href="https://github.com/Solrikk/BioSynth/blob/main/README.md">⭐ English ⭐</a> |
    <a href="https://github.com/Solrikk/BioSynth/blob/main/README_RU.md">🇷🇺 Russian</a> |
    <a href="https://github.com/Solrikk/BioSynth/blob/main/README_GE.md">🇩🇪 German</a> |
    <a href="https://github.com/Solrikk/BioSynth/blob/main/README_JP.md">🇯🇵 Japanese</a> |
    <a href="https://github.com/Solrikk/BioSynth/blob/main/README_KR.md">🇰🇷 Korean</a> |
    <a href="https://github.com/Solrikk/BioSynth/blob/main/README_CN.md">🇨🇳 Chinese</a>
  </h3>
</div>

-----------------

# BioSynth 🔬

_**BioSynth**_ is an innovative tool designed to simulate and visualize protein synthesis from DNA sequences. It integrates modern bioinformatics methods to provide a detailed understanding of the relationship between DNA sequences and their corresponding protein structures.

## Features ✨

- **🧬 DNA Transcription and Translation**: Converts DNA sequences into RNA sequences and subsequently translates them into protein sequences following the central dogma of molecular biology. [Learn more about Transcription and Translation](https://www.ncbi.nlm.nih.gov/books/NBK26887/).
- **🧪 Mutation Simulation**: Introduces random mutations in the DNA sequence to study their effects on the resulting protein. [Read about Mutation Types](https://en.wikipedia.org/wiki/Mutation).
- **🖼️ 3D Visualization**: Generates and visualizes 3D structures of proteins and DNA-protein complexes using py3Dmol. [Py3Dmol Documentation](https://3dmol.csb.pitt.edu/doc/).
- **⚡ Asynchronous Processing**: Utilizes Python's asyncio for efficient handling of multiple synthesis tasks concurrently. [Asyncio Official Documentation](https://docs.python.org/3/library/asyncio.html).

## How it Works

BioSynth performs several key steps to transform a DNA sequence into a visualized 3D modeled protein. In this section, we'll delve into each of these steps in detail:

### 1. 📝 Transcription
Transcription is the first step in protein synthesis, where the DNA sequence is transcribed into a messenger RNA (mRNA) sequence. This process involves replacing thymine (T) in the DNA with uracil (U) in the RNA.
- **Process**: The DNA sequence is unwound, and enzymes copy one of the strands to form mRNA.
- **Example**:
  ```python
  dna_sequence = "ATGCGTACGTTAG"
  rna_sequence = dna_sequence.replace("T", "U")
  print("RNA Sequence:", rna_sequence)
  ```
2. **🔄 Translation**: The RNA sequence is then translated into a protein sequence based on codon tables. [A Guide to Translation](https://www.ncbi.nlm.nih.gov/books/NBK9849/).
3. **🧬 Protein Mutation**: A random base in the DNA sequence is mutated to simulate real-world genetic variations. [Understanding Genetic Mutations](https://www.genome.gov/Genetic-Disorders/What-are-Genetic-Mutations).
4. **🏗️ Structure Generation**: 3D structures of the synthesized protein and DNA-protein complexes are generated and saved as PDB files. [PDB File Format](https://www.rcsb.org/pages/help/format/pdbformat).
5. **🔍 Visualization**: The generated structures are visualized in an interactive 3D viewer. [3D Molecular Visualization Tools](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6168215/).
