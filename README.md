![Logo](https://github.com/Solrikk/BioSynth/blob/main/assets/OpenCV%20-%20result/bee.jpg)

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

- **🧬 DNA Transcription and Translation**: Converts DNA sequences into RNA sequences and subsequently translates them into protein sequences following the central dogma of molecular biology.
- **🧪 Mutation Simulation**: Introduces random mutations in the DNA sequence to study their effects on the resulting protein.
- **🖼️ 3D Visualization**: Generates and visualizes 3D structures of proteins and DNA-protein complexes using py3Dmol.
- **⚡ Asynchronous Processing**: Utilizes Python's asyncio for efficient handling of multiple synthesis tasks concurrently.

## How it Works

1. **📝 Transcription**: The DNA sequence is transcribed into an RNA sequence by replacing thymine (T) with uracil (U).
2. **🔄 Translation**: The RNA sequence is then translated into a protein sequence based on codon tables.
3. **🧬 Protein Mutation**: A random base in the DNA sequence is mutated to simulate real-world genetic variations.
4. **🏗️ Structure Generation**: 3D structures of the synthesized protein and DNA-protein complexes are generated and saved as PDB files.
5. **🔍 Visualization**: The generated structures are visualized in an interactive 3D viewer.
