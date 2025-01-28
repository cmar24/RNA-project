# -*- coding: utf-8 -*-
"""
Created on Fri Sep 2024
@project: Bioinformatics of RNA. GENIOMHE M2
@author: Cristina Marianini
@year: 2024-2025
"""
import math
import numpy as np
from collections import defaultdict
import os

# Import the required functions from other scripts
from training import compute_distances, update_counts, compute_frequencies, compute_scores, save_scores
from utils.pdb_processing import process_directory
from scoring import  read_scores, plot_interaction_profiles

if __name__ == "__main__":
    # Directory containing PDB files
    directory = input("Enter the path to the directory containing PDB files: ")
    
    # Process all files in the directory
    print(f"Processing PDB files in the directory: {directory}")
    c3_atoms = process_directory(directory)
    print(f"Finished processing files. Total C3' atoms processed: {len(c3_atoms)}")

    # Continue with distance calculations and scoring
    print("Starting distance calculations...")
    total_distances = 0
    for base1, base2, distance in compute_distances(c3_atoms):
        update_counts(base1, base2, distance)
        total_distances += 1

    print(f"Distance calculations complete. Total distances processed: {total_distances}")

    # Compute frequencies and scores
    print("Computing frequencies and scores...")
    observed_freqs, reference_freqs = compute_frequencies()
    scores = compute_scores(observed_freqs, reference_freqs)

    # Save results
    print("Saving scoring results...")
    save_scores(scores)
    print("Training complete. Scoring files saved.")

    # Directory containing the score files
    score_directory = "data/output"
    plot_output_dir = "output/plots"  # Directory to save plots
    
    # Create the output directory if it doesn't exist
    os.makedirs(plot_output_dir, exist_ok=True)
    
    # Read scores and create a DataFrame
    print(f"Reading scores from {score_directory}...")
    scores_df = read_scores(score_directory)
    
    # Plot interaction profiles
    print("Generating and saving plots...")
    plot_interaction_profiles(scores_df, plot_output_dir)
    
    print(f"All plots saved in: {plot_output_dir}")