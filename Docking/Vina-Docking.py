from vina import Vina
import os

# Initialize Vina
v = Vina()

# Set the receptor file
v.set_receptor(rigid_pdbqt_filename="/content/drive/MyDrive/Colab Notebooks/fwddockingfiles/IL15RA.pdbqt")

# Path to the directory containing ligand files
ligand_dir = '/content/drive/MyDrive/Colab Notebooks/fwddockingfiles/ligands'
output_dir = '/content/drive/MyDrive/Colab Notebooks/fwddockingfiles/docking_results'

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# List all ligand files in the directory
ligand_files = [f for f in os.listdir(ligand_dir) if f.endswith('.pdbqt')]

# Docking box parameters
center = [10.628, 3.674, -3.025]
box_size = [126, 76, 100]

# Loop through each ligand file and perform docking
for ligand_file in ligand_files:
    ligand_path = os.path.join(ligand_dir, ligand_file)
    print(f"Processing ligand: {ligand_file}")

    # Set the ligand file
    v.set_ligand_from_file(ligand_path)

    # Compute Vina maps
    v.compute_vina_maps(center=center, box_size=box_size)

    # Score and optimize
    score = v.score()
    print(f"Score for {ligand_file}: {score}")
    optimize_result = v.optimize()
    print(f"Optimize result for {ligand_file}: {optimize_result}")

    # Perform docking with exhaustiveness
    v.dock(exhaustiveness=32)  # Removed `num_modes` and `energy_range` arguments

    # Write docking results
    output_file = os.path.join(output_dir, f"docking_results_{ligand_file}")
    v.write_poses(pdbqt_filename=output_file)
    print(f"Docking results written to: {output_file}")