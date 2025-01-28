import os
import pandas as pd
import matplotlib.pyplot as plt

def read_scores(directory):
    """
    Reads all score files in the specified directory and returns a combined DataFrame.
    Assumes each file is named <BasePair>_scores.txt and contains scores per distance bin.
    """
    all_scores = []
    
    # Loop through all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith("_scores.txt"):
            base_pair = filename.replace("_scores.txt", "")  # Extract base pair name
            filepath = os.path.join(directory, filename)
            
            # Read scores and create a DataFrame for this base pair
            with open(filepath, 'r') as f:
                scores = [float(line.strip()) for line in f]
            
            df = pd.DataFrame({
                "Distance_Bin": range(len(scores)),  # Bin index (e.g., 0-1 Å, 1-2 Å, etc.)
                "Score": scores,
                "Base_Pair": base_pair
            })
            all_scores.append(df)
    
    # Combine all DataFrames into one
    return pd.concat(all_scores, ignore_index=True)

def plot_interaction_profiles(scores_df, output_dir):
    """
    Plots interaction profiles (score vs. distance) for all base pairs.
    Saves a general plot and individual plots for each base pair.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # General Plot: All base pairs
    print("Generating general interaction profile plot...")
    plt.figure(figsize=(10, 6))
    for base_pair, group in scores_df.groupby("Base_Pair"):
        plt.plot(group["Distance_Bin"], group["Score"], label=base_pair, linewidth=2)
    
    plt.title("Interaction Profiles: Score vs Distance", fontsize=16)
    plt.xlabel("Distance (Å)", fontsize=14)
    plt.ylabel("Score", fontsize=14)
    plt.legend(title="Base Pair", fontsize=12)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    general_plot_path = os.path.join(output_dir, "general_interaction_profile.png")
    plt.savefig(general_plot_path, dpi=300)
    plt.close()
    print(f"General plot saved to {general_plot_path}")
    
    # Individual Plots: One per base pair
    print("Generating individual plots for each base pair...")
    for base_pair, group in scores_df.groupby("Base_Pair"):
        plt.figure(figsize=(8, 5))
        plt.plot(group["Distance_Bin"], group["Score"], label=base_pair, color="b", linewidth=2)
        
        plt.title(f"Interaction Profile: {base_pair}", fontsize=16)
        plt.xlabel("Distance (Å)", fontsize=14)
        plt.ylabel("Score", fontsize=14)
        plt.legend(fontsize=12)
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.tight_layout()
        
        # Save the individual plot
        individual_plot_path = os.path.join(output_dir, f"{base_pair}_interaction_profile.png")
        plt.savefig(individual_plot_path, dpi=300)
        plt.close()
        print(f"Plot for {base_pair} saved to {individual_plot_path}")
