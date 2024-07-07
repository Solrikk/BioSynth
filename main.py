import asyncio
from synthesis import synthesize_until_five_successful, visualize_protein

if __name__ == "__main__":
  dna_sequence = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTCGTCCGGGCAGGGCAGGGCGATGCCACCCACCAGGTACCGAAGCGCGA"
  asyncio.run(synthesize_until_five_successful(dna_sequence))

  pdb_files = [f"results/protein_dna_structure_{i}.pdb" for i in range(1, 6)]
  visualize_protein(pdb_files)
