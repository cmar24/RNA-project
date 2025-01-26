# Creation of an objective function for the RNA folding problem
## Project description
This project implements an objective function to estimate the Gibbs free energy for RNA folding, a key aspect of predicting the native structure of a ribonucleotide chain. The native fold corresponds to the conformation with the lowest Gibbs free energy, and this project models it using interatomic distance distributions.

The repository contains three Python scripts to:

1. **Train the objective** function using interatomic distance distributions derived from experimentally determined RNA 3D structures.
2. **Visualize scoring profiles**, displaying scores as a function of interatomic distances.
3. **Evaluate RNA structures** by estimating their Gibbs free energy using the objective function

## Install dependencies

## Data
The project includes a set of data to allow the training process. The data was downloaded at https://www.rcsb.org/search/advanced where I performed an advanced search with the following criteria:
- Polymer Entity Type is RNA
- Experimental Methos is X-Ray Diffraction

Then I downloaded the first 200 files.

## Files and directories


## Usage

Were to find predicted RNAs
- RNA ouzzles
- CASP
