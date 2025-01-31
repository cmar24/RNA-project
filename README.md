# Creation of an objective function for the RNA folding problem
## Project description
This project implements an objective function to estimate the Gibbs free energy for RNA folding, a critical factor in predicting the native structure of a ribonucleotide chain. The native fold corresponds to the conformation with the lowest Gibbs free energy, and this project models it using interatomic distance distributions.

The repository contains three main Python scripts:

- `training.py`: Trains the objective function by deriving interatomic distance distributions from experimentally determined RNA 3D structures.
- `scoring.py`: Visualizes scoring profiles, displaying scores as a function of interatomic distances.
- `evaluation.py`: Evaluates RNA structures by estimating their Gibbs free energy using the trained objective function.

This project aims to assist RNA structure prediction by providing a systematic approach to calculate the folding energy, facilitating the identification of the native structure.

## Installation
Before running the scripts, ensure you have the required dependencies installed. You can install them using pip:
``` python
pip install -r requirements.txt
```
## Data
The training data consists of RNA 3D structures downloaded from the Protein Data Bank (PDB). To generate the dataset:

1. Visit the PDB advanced search: RCSB PDB Advanced Search.
https://www.rcsb.org/search/advanced 
2. Apply the following search criteria:
- Polymer Entity Type: RNA
- Experimental Method: X-Ray Diffraction
3. Download the first 100 structures from the search results.

Save the downloaded files into a directory (e.g., data/pdb_files), which will be used for training the objective function.

## Usage
Run the `main.py` script to process the PDB files, compute distance distributions, and train the objective function:

```bash
python main.py
```
Process the RNA 3D structures in the data/pdb_files directory.
Compute interatomic distances for base pairs.
Calculate observed and reference frequencies of distances.
Compute scores (log-ratios) and save them in data/output.

## Sources of Predicted RNA Structures
To evaluate the objective function on predicted RNA structures, you can use datasets from:

- RNA Puzzles: RNA Puzzles Website
An open competition that provides RNA structure predictions and their experimentally determined counterparts for benchmarking.
- CASP (Critical Assessment of Structure Prediction): CASP Website.
A similar competition focused on protein structure predictions but with some RNA challenges.

I went to the second website and downloaded a set of RNA predictions from CAP15 to evaluate them (included in the repository).