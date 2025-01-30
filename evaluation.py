import os
import math
import numpy as np
from utils.pdb_processing import read_c3_pdb_file
from training import compute_distances

# Define constants
SCORE_DIRECTORY = "data/output"
DISTANCE_BINS = np.linspace(0, 20, 21)  # 20 bins (0-20 Å)

def load_scores():
    """Load scoring profiles from training output files."""
    scores = {}
    for filename in os.listdir(SCORE_DIRECTORY):
        if filename.endswith("_scores.txt"):
            base_pair = filename.replace("_scores.txt", "")
            with open(os.path.join(SCORE_DIRECTORY, filename), 'r') as f:
                score_values = [float(line.strip()) for line in f]
            scores[base_pair] = np.array(score_values)
    return scores

def interpolate_score(base1, base2, distance, scores):
    """Find the interpolated score for a given base pair and distance."""
    pair_type = ''.join(sorted([base1, base2]))  # Sort bases alphabetically
    if pair_type not in scores:
        return 0  # Default score if no data for this pair type
    
    # Get the corresponding score array
    score_values = scores[pair_type]
    
    # Find the two closest bins for interpolation
    bin_idx = np.digitize(distance, DISTANCE_BINS) - 1
    if bin_idx < 0:
        return score_values[0]
    if bin_idx >= len(score_values) - 1:
        return score_values[-1]
    
    # Linear interpolation
    d1, d2 = DISTANCE_BINS[bin_idx], DISTANCE_BINS[bin_idx + 1]
    s1, s2 = score_values[bin_idx], score_values[bin_idx + 1]
    return s1 + (s2 - s1) * (distance - d1) / (d2 - d1)

def evaluate_structure(pdb_file):
    """Evaluate an RNA structure by computing the estimated Gibbs free energy."""
    print(f"Evaluating structure: {pdb_file}")
    
    # Load C3' atomic coordinates
    c3_atoms = read_c3_pdb_file(pdb_file)
    
    # Load trained scores
    scores = load_scores()
    
    # Compute total estimated energy
    total_score = 0
    interaction_count = 0

    for base1, base2, distance in compute_distances(c3_atoms):
        energy = interpolate_score(base1, base2, distance, scores)
        total_score += energy
        interaction_count += 1

    print(f"Evaluation complete. Total interactions: {interaction_count}")
    print(f"Estimated Gibbs Free Energy: {total_score:.4f}")
    return total_score
