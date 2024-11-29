# -*- coding: utf-8 -*-
"""
Created on Fri Sep 2024
@project: Bioinformatics of RNA. GENIOMHE M2
@author: Cristina Marianini
@year: 2024-2025
"""

## 1. Training

import math
import numpy as np

# Download data and move it to your working directory

def read_c3_pdb_file():
    # Ask the user for the pdb filename
    filename = input("Please enter the PDB filename: ")
    
    # List to store the coordinates and residue type of C3 atoms
    c3_atoms = []
    
    try:
        # Open and read the file
        with open(filename, 'r') as file:
            for line in file:
                # Check if the line starts with "ATOM" and contains " C3' "
                if line.startswith("ATOM") and " C3' " in line:
                    # Split the line by spaces to extract all elements
                    elements = line.split()
                    
                    # Extract residue name and coordinates 
                    residue_name = elements[3]  # Residue name
                    x, y, z = map(float, elements[6:9])  # Coordinates xyz
                    chain = elements[4]

                    # Append the tuple to the list
                    c3_atoms.append((residue_name, x, y, z, chain))
        
        # Return the list of C3 atom coordinates and residue types
        return c3_atoms

    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return []
    
    except Exception as e:
        print(f"Error: An error occurred while reading the file. {e}")
        return []

# Run the function
c3_atoms = read_c3_pdb_file()
print(c3_atoms)

# Define base pair types
base_pairs = ["AA", "AU", "AC", "AG", "UU", "UC", "UG", "CC", "CG", "GG"]
distance_counts = {bp: [] for bp in base_pairs} # Dictionary to store distances

def order_base_pair(res1, res2):
    """ Order the pairs to store in the dictionnary. """
    return ''.join(sorted([res1, res2]))

def compute_intrachain_distances(c3_atoms):
    n = len(c3_atoms)
    print(n)
    for i in range(n):
        for j in range( i + 4 , n): # Consider only i and j where j >= i + 4
            print(c3_atoms[i][4],c3_atoms[j][4])
            if c3_atoms[i][4] == c3_atoms[j][4]: # Consider only intrachain distances
                res1, res2 = c3_atoms[i], c3_atoms[j]
                distance = math.sqrt((res1[1]-res2[1])**2 + (res1[2]-res2[2])**2 + (res1[3]-res2[3])**2)
                base_pair_type = order_base_pair(res1[0], res2[0])
                
                if base_pair_type in distance_counts: # confirm that the pair exits in the dict
                    distance_counts[base_pair_type].append(distance)

compute_intrachain_distances(c3_atoms)
print(distance_counts)







