import time
import asyncio
import random
import math
import os
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue, Atom
import py3Dmol

CODON_TABLE = {
    'AUG': 'M',  # Methionine
    'UUU': 'F',
    'UUC': 'F',  # Phenylalanine
    'UUA': 'L',
    'UUG': 'L',  # Leucine
    'UCU': 'S',
    'UCC': 'S',
    'UCA': 'S',
    'UCG': 'S',  # Serine
    'UAU': 'Y',
    'UAC': 'Y',  # Tyrosine
    'UGU': 'C',
    'UGC': 'C',  # Cysteine
    'UGG': 'W',  # Tryptophan
    'UAA': 'STOP',
    'UAG': 'STOP',
    'UGA': 'STOP'  # Stop codons
}

DNA_TO_RNA = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
RNA_TO_DNA = {v: k for k, v in DNA_TO_RNA.items()}


def log(message):
  print(f"[LOG] {message}")


def generate_dna_coords(length,
                        radius=10,
                        rise_per_base=3.4,
                        twist_per_base=36):
  coords_a = []
  coords_b = []
  for i in range(length):
    angle = math.radians(i * twist_per_base)
    x = radius * math.cos(angle)
    y = radius * math.sin(angle)
    z = i * rise_per_base
    coords_a.append((x, y, z))
    coords_b.append((-x, -y, z))
  return coords_a, coords_b


async def transcribe(dna_sequence):
  return ''.join([DNA_TO_RNA[nuc] for nuc in dna_sequence])


async def translate(rna_sequence):
  protein_sequence = []
  for i in range(0, len(rna_sequence) - 2, 3):
    codon = rna_sequence[i:i + 3]
    amino_acid = CODON_TABLE.get(codon, 'STOP')
    if amino_acid == 'STOP':
      break
    protein_sequence.append(amino_acid)
  return protein_sequence


def mutate_dna(dna_sequence):
  dna_bases = list(dna_sequence)
  mutation_index = random.randint(0, len(dna_bases) - 1)
  original_base = dna_bases[mutation_index]
  bases = ['A', 'T', 'C', 'G']
  bases.remove(original_base)
  dna_bases[mutation_index] = random.choice(bases)
  mutated_sequence = ''.join(dna_bases)
  return mutated_sequence


async def synthesize_protein(dna_sequence):
  log("Starting protein synthesis")
  start_time = time.time()
  try:
    rna_sequence = await transcribe(dna_sequence)
    log(f"Transcribed RNA sequence: {rna_sequence}")
    protein_sequence = await translate(rna_sequence)
    log(f"Translated protein sequence: {protein_sequence}")
  except (TypeError, ValueError) as e:
    return f"Error during synthesis: {e}"
  end_time = time.time()
  synthesis_time = end_time - start_time
  mutated_dna = mutate_dna(dna_sequence)
  return {
      "DNA Sequence": dna_sequence,
      "RNA Sequence": rna_sequence,
      "Protein Sequence": protein_sequence,
      "Synthesis Time (seconds)": synthesis_time,
      "Mutated DNA Sequence": mutated_dna
  }


async def synthesize_until_five_successful(dna_sequence):
  results = []
  iteration = 0
  while len(results) < 5:
    iteration += 1
    log(f"Iteration: {iteration}")
    result = await synthesize_protein(dna_sequence)
    if isinstance(result, dict):
      results.append(result)
      log(f"Successful synthesis #{len(results)}")
      save_protein_structure(result, iteration)
      dna_sequence = result["Mutated DNA Sequence"]
    else:
      log(f"Error encountered: {result}")
  log(f"Total successful results collected: {len(results)}")
  return results


def save_protein_structure(result, iteration):
  protein_sequence = result['Protein Sequence']
  dna_sequence = result['DNA Sequence']

  structure = Structure.Structure(f"Protein_{iteration}")
  model = Model.Model(iteration)
  chain = Chain.Chain("A")
  atom_serial = 1

  for idx, aa in enumerate(protein_sequence):
    res_id = (' ', idx + 1, ' ')
    res = Residue.Residue(res_id, aa, idx + 1)
    coords = (idx * 1.5, 0.0, 0.0)
    atom_ca = Atom.Atom(name="CA",
                        coord=coords,
                        bfactor=0,
                        occupancy=1,
                        altloc=' ',
                        fullname=' CA ',
                        serial_number=atom_serial,
                        element='C')
    atom_serial += 1
    res.add(atom_ca)
    chain.add(res)

  dna_chain_a = Chain.Chain("B")
  dna_chain_b = Chain.Chain("C")
  coords_a, coords_b = generate_dna_coords(len(dna_sequence))
  for idx, nucleotide in enumerate(dna_sequence):
    res_id_a = (' ', idx + 1, ' ')
    res_id_b = (' ', idx + 1, ' ')

    res_a = Residue.Residue(res_id_a, nucleotide, idx + 1)
    res_b = Residue.Residue(res_id_b, RNA_TO_DNA[DNA_TO_RNA[nucleotide]],
                            idx + 1)

    atom_p_a = Atom.Atom(name="P",
                         coord=coords_a[idx],
                         bfactor=0,
                         occupancy=1,
                         altloc=' ',
                         fullname=' P  ',
                         serial_number=atom_serial,
                         element='P')
    atom_serial += 1

    atom_p_b = Atom.Atom(name="P",
                         coord=coords_b[idx],
                         bfactor=0,
                         occupancy=1,
                         altloc=' ',
                         fullname=' P  ',
                         serial_number=atom_serial,
                         element='P')
    atom_serial += 1

    res_a.add(atom_p_a)
    res_b.add(atom_p_b)

    dna_chain_a.add(res_a)
    dna_chain_b.add(res_b)

  model.add(chain)
  model.add(dna_chain_a)
  model.add(dna_chain_b)
  structure.add(model)

  if not os.path.exists('results'):
    os.makedirs('results')

  io = PDBIO()
  io.set_structure(structure)
  file_name = f"results/protein_dna_structure_{iteration}.pdb"
  io.save(file_name)
  log(f"Saved protein structure with DNA to {file_name}")


def visualize_protein(pdb_files):
  view = py3Dmol.view(width=800, height=600, linked=True)
  for i, pdb_file in enumerate(pdb_files):
    with open(pdb_file) as f:
      pdb_code = f.read()
    view.addModel(pdb_code, 'pdb')
    view.setStyle({'model': i + 1}, {
        'cartoon': {
            'color': 'spectrum'
        },
        'stick': {}
    })
  view.zoomTo()
  return view.show()
