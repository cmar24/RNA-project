import os
import gzip

def process_directory(directory):
    """Processes all PDB or compressed PDB files in a directory."""
    c3_atoms_list = []
    
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            # Check if the file is a compressed PDB (.pdb.gz) or plain PDB (.pdb)
            if file.endswith(".pdb"):
                c3_atoms = read_c3_pdb_file(file_path)
            elif file.endswith(".pdb.gz"):
                c3_atoms = read_compressed_pdb_file(file_path)
            else:
                print(f"Skipping unsupported file: {file}")
                continue
            
            if c3_atoms:
                c3_atoms_list.extend(c3_atoms)  # Aggregate results
    return c3_atoms_list

def read_compressed_pdb_file(filepath, atom_name="C3'"):
    """Reads a compressed PDB file (.pdb.gz) and extracts atom information."""
    c3_atoms = []
    try:
        with gzip.open(filepath, 'rt') as file:
            for line in file:
                if line.startswith("ATOM") and line[12:16].strip() == atom_name:
                    try:
                        residue_name = line[17:20].strip()
                        chain = line[21].strip()
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        c3_atoms.append((residue_name, x, y, z, chain))
                    except ValueError:
                        print(f"Error parsing line in {filepath}: {line.strip()}")
                        continue
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
    return c3_atoms

def read_c3_pdb_file(atom_name="C3'"):
    """
    Reads a PDB file and extracts information for the specified atom type (default: C3').

    Returns:
        List of tuples: Each tuple contains (residue_name, x, y, z, chain_id)
        where:
        - residue_name (str): Residue name (e.g., A, U, C, G)
        - x, y, z (float): Coordinates of the atom
        - chain_id (str): Chain identifier
    """
    # Ask the user for the PDB filename
    filename = input("Please enter the PDB filename: ")

    # List to store extracted atom information
    c3_atoms = []

    try:
        with open(filename, 'r') as file:
            for line in file:
                # Check for ATOM lines and match the atom name
                if line.startswith("ATOM") and line[12:16].strip() == atom_name:
                    try:
                        # Extract residue name, chain, and coordinates using fixed-width slicing
                        residue_name = line[17:20].strip()  # Residue name (A, U, C, G)
                        chain = line[21].strip()            # Chain ID
                        x = float(line[30:38].strip())      # X coordinate
                        y = float(line[38:46].strip())      # Y coordinate
                        z = float(line[46:54].strip())      # Z coordinate

                        # Append as a tuple
                        c3_atoms.append((residue_name, x, y, z, chain))
                    except ValueError:
                        print(f"Error parsing coordinates in line: {line.strip()}")
                        continue

        if not c3_atoms:
            print(f"No {atom_name} atoms found in the file.")
        return c3_atoms

    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return []
    except Exception as e:
        print(f"Error: An unexpected error occurred: {e}")
        return []
        