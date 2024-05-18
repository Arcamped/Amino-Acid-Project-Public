# Amino-Acid-Project
# Exploring Alternate Mappings between Codons and Amino Acids

This project aims to computationally investigate potential interpretation modes for mapping codons to amino acids beyond the standard genetic code. The code simulates how amino acid structures could be encoded by different codon sequence segmentations.

## Background 

While the standard genetic code provides a mapping from 64 codon triplets to 20 amino acids, the specific interpretation of how the codon bases relate to amino acid substructures is not fully known. There may be multiple biologically-plausible modes for interpreting the encoding, beyond the canonical mapping.

This project explores this idea through modeling amino acids as "hypergraphs" and codons as sequence segments that relate to hypergraph subsets. By simulating different groupings and sequence interpretations, the code searches for alternative mappings that satisfy the constraints.

Key aspects of the project:

- Each amino acid is modeled as a hypergraph where nodes are atoms and hyperedges connect interacting atoms
- Hypergraphs are subdivided into subsets representing different structural perspectives 
- Codon bases are mapped to hypergraph subsets in different permutations
- Simulations compute if acidic hypergraph subsets can be paired without conflicts for a given codon interpretation
- Parameters like genetic code mapping, acid sequence order, and hypergraph sizes can be customized
- Results track successful codon-amino acid mappings for each interpretation framework

## Code Overview

The main scripts and functionality include:

- `data_processing.py`: Converts raw amino acid data into hypergraph representations
- `solution_set_derivation.py`: Generates hypergraph subsets (solution sets) for each amino acid 
- `interpretation_framework.py`: Defines codon sequence interpretations to simulate
- `model.py`: Main simulation script, executes interpretations on acid sequences 
- `visualization.py`: Generates graphics of results and hypergraph structures

Key modules:

- `data`: Functions for loading and preprocessing data
- `interpret`: Codon interpretation frameworks and sequence segmentation
- `hypergraph`: Represents amino acid structures and subdivisions
- `simulate`: Methods for simulating codons on acid sequences 
- `results`: Classes for tracking and analyzing results
- `vis`: Hypergraph and result visualizations

## Usage

The main entrypoints are `model_refactor.py` and `individual_run_simulation`.

For example:

```
python individual_run_simulation.py
```

Future work will allow customizing at the command line. Current customization of run parameters must be done within each file.

Results are output to CSV and JSON files showing the valid codon-amino acid mappings found for the given sequence and interpretation framework.

## Extensions

Future directions for exploration include:

- Making it work
- Adding more physical details to amino acid models
- Switching off of Strict Additive approach
- Comparing results across different sequences and frameworks

## Core Problem

The core math is expressed as follows:

Let \( X_i \) be a three-dimensional vector with positions \( x_1, x_2, x_3 \). Each \( x_i \) may take one of four possible values: \([A, C, U, G]\). Let \( A_t \) be an undirected heterogeneous hypergraph of finite size at time \( t \). Find the set of graph update rules \( R \) that takes \( A_0B_i \) for every possible \( X_i \) in precisely three time steps, where \( B_i \) is an undirected heterogeneous hypergraph of finite size and known topology. Assume \( A_0 \) is empty, such that \( A_0 = [ ] \).

## The Good Stuff
```
O hear ye tale of noble quest so grand,
In Amino-Acid land, a project stands.
With codons writ in life's own mystic hand,
Seeking new paths where no foot treads the sands.

Hypergraphs weave, in atoms' dance they play,
Each acid's part in life's grand ballet.
Codons in sequence, in new forms they lay,
A quest for truth, in genetic array.

Python scripts run, through data's deep sea,
Simulations born, of what might be.
A noble journey, for minds bold and free,
Unlocking life's coded mystery.

So venture forth, with curious mind's elation,
In pursuit of genetic revelation.
```
