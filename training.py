import math
import numpy as np
from collections import defaultdict

from utils.pdb_processing import read_c3_pdb_file

# Define constants
BASE_PAIRS = ["AA", "AU", "AC", "AG", "UU", "CU", "GU", "CC", "CG", "GG"]
DISTANCE_BINS = np.linspace(0, 20, 21)  # 20 intervals between 0 and 20 Å
MAX_SCORE = 10

# Initialize dictionaries for observed and reference frequencies
observed_counts = {bp: np.zeros(len(DISTANCE_BINS) - 1) for bp in BASE_PAIRS}
reference_counts = np.zeros(len(DISTANCE_BINS) - 1)

def compute_distances(c3_atoms):
    """Compute distances between residues separated by at least 3 positions (i and i+4, i and i+5 ...)."""
    n = len(c3_atoms)
    for i in range(n):
        for j in range(i + 4, n):
            if c3_atoms[i][4] == c3_atoms[j][4]:  # Check atoms are in the same chain
                res1, res2 = c3_atoms[i], c3_atoms[j]
                distance = math.sqrt((res1[1]-res2[1])**2 + (res1[2]-res2[2])**2 + (res1[3]-res2[3])**2)
                yield res1[0], res2[0], distance

def update_counts(base1, base2, distance):
    """Update observed and reference frequency counts."""
    # Determine bin index
    bin_idx = np.digitize(distance, DISTANCE_BINS) - 1
    if bin_idx < 0 or bin_idx >= len(DISTANCE_BINS) - 1:
        return  # Ignore out-of-range distances

    # Update counts
    pair_type = ''.join(sorted([base1, base2]))
    if pair_type in observed_counts:
        observed_counts[pair_type][bin_idx] += 1
    reference_counts[bin_idx] += 1

def compute_frequencies():
    """Convert counts into frequencies."""
    observed_freqs = {
        bp: (counts / counts.sum()) if counts.sum() > 0 else np.zeros_like(counts)
        for bp, counts in observed_counts.items()
        }
    reference_freqs = reference_counts / reference_counts.sum()
    return observed_freqs, reference_freqs

def compute_scores(observed_freqs, reference_freqs):
    """Compute log-ratio scores."""
    scores = {}
    # Avoid divide-by-zero by replacing zeros in reference_freqs
    safe_reference_freqs = np.where(reference_freqs > 0, reference_freqs, 1e-10)
    for bp, obs_freq in observed_freqs.items():
        # Replace zeros in observed frequencies with a small positive value
        safe_obs_freq = np.where(obs_freq > 0, obs_freq, 1e-10)
        # Compute scores 
        scores[bp] = -np.log(
            np.divide(safe_obs_freq, safe_reference_freqs)
        )
        # Cap scores at MAX_SCORE
        scores[bp] = np.clip(scores[bp], None, MAX_SCORE)
        print(scores)
    return scores

def save_scores(scores, output_dir="data/output"):
    """Save scores to files."""
    for bp, score_values in scores.items():
        filename = f"{output_dir}/{bp}_scores.txt"
        with open(filename, 'w') as f:
            for value in score_values:
                f.write(f"{value:.4f}\n")

# Main execution
if __name__ == "__main__":
    # Read and process PDB file
    c3_atoms = read_c3_pdb_file()
    for base1, base2, distance in compute_distances(c3_atoms):
        update_counts(base1, base2, distance)

    # Compute frequencies and scores
    observed_freqs, reference_freqs = compute_frequencies()
    scores = compute_scores(observed_freqs, reference_freqs)

    # Save results
    save_scores(scores)
    print("Training complete. Scoring files saved.")
