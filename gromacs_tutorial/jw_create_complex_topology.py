import argparse

def create_complex_topology(ligand_file, protein_file,  temp_directory, output_directory):
    """
    Writes a new .gro file which is a combination of the protein and ligand structures.
    Also converts ligand.top file to ligand.itp file.
    :param ligand_file:
    :param protein_file:
    :param temp_directory:
    :param output_directory:
    :return:
    """

    # Open protein and ligand structure files
    protein_structure = open(f'{output_directory}/{protein_file}.gro', 'r')
    ligand_structure = open(f'{output_directory}/{ligand_file}.gro', 'r')

    # Create new file for the complex
    complex_structure = open(f'{output_directory}/complex.gro', 'w')

    # Open ligand topology files for reading and writing
    ligand_topology_file = open(f'{output_directory}/{ligand_file}.top', 'r')
    new_ligand_topology_file = open(f'{output_directory}/{ligand_file}.itp', 'w')

    # Read lines from protein and ligand structure files
    protein_lines = protein_structure.readlines()
    ligand_lines = ligand_structure.readlines()

    # Determine the total number of atoms in the new structure
    protein_num_atoms = int(protein_lines[1].strip())
    ligand_num_atoms = int(ligand_lines[1].strip())
    total_num_atoms = protein_num_atoms + ligand_num_atoms

    # Get the coordinates of the protein structure
    protein_boundary_coordinates = protein_lines[-1]

    # Write header and structure data to the new complex structure file
    header = "Protein Complex for ABFE\n"
    new_protein_structure = protein_lines[2:-1]
    new_ligand_structure = ligand_lines[2:-1]
    complex_structure.write(header)
    complex_structure.write(str(total_num_atoms) + "\n")
    complex_structure.writelines(new_protein_structure)
    complex_structure.writelines(new_ligand_structure)
    complex_structure.write(protein_boundary_coordinates)
    complex_structure.close()

    # Read and write the ligand topology file
    badlines = [
        "[ defaults ]",
        "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ",
        "1               2               yes             0.5          0.83333333"
    ]
    ligand_topology = []

    # Copy all lines from the input ligand topology file except for the [ defaults ], nonbonded_function,
    # and nonbonded_parameters sections
    for line in ligand_topology_file.readlines():
        if any([badline in line for badline in badlines]):
            continue
        else:
            ligand_topology.append(line)

    # Get the name of the ligand from the last line of the topology file and remove the last 7 lines
    # (which contain the [ system ], [ molecules ], and [ position_restraints ] sections)
    ligand_name = ligand_topology[-1]
    new_ligand_topology = ligand_topology[0:-7]
    new_ligand_topology_file.writelines(new_ligand_topology)
    new_ligand_topology_file.close()

    # Read and write the protein topology file
    force_field_include = '#include "./amber99sb-ildn_zn.ff/forcefield.itp"'
    ligand_topology_include = f'\n; Include Ligand Topology\n#include "{ligand_file}.itp"\n'
    new_topology = []
    with open(f'{output_directory}/topol.top', 'r') as f:
    # Copy all lines from the input topology file except for the line with the force field include,
    # and add a line to include the ligand topology file
        for i in f.readlines():
            if force_field_include in i:
                new_topology.append(i + ligand_topology_include)
            else:
                new_topology.append(i)
    new_topology.append(ligand_name)
    with open(f'{output_directory}/topol.top', 'w') as f:
        f.writelines(new_topology)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--ligand', type=str, help='name of ligand')
    parser.add_argument('-p', '--protein', type=str, help='name of protein')
    parser.add_argument('-t', '--temp',  default='.', type=str, help='name of temp directory')
    parser.add_argument('-o', '--output', default='.', type=str, help='name of output directory')
    args = parser.parse_args()
    create_complex_topology(args.ligand, args.protein, args.temp, args.output)
