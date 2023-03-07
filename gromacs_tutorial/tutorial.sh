# download the charmm forcefield and
# cgenff_charmm2gmx_py3_nx2.py from here:
# http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs

# write the topology for the protein using pdb2gmx
export PATH=/Applications/gromacs-2022/build/bin:$PATH
pwd
# we use the charm36 FF which must be in the local dir, tip3p is default water model
gmx pdb2gmx -f 3HTB_clean.pdb -o 3HTB_processed.gro -ff charmm36_mar2019 -water tip3p

# In Daniel's script we can either use mol2 or sdf - pdbqt must be converted
obabel jz4.pdb -O jz4.sdf -p 7.4

# prepare the ligand with forcefields
antechamber -i path/to/mol.sdf -fi sdf -o path/to/output/mol.mol2 -fo mol2 -pf y -s 0 -j 5 -c bcc
# If command above fails consider adding following argument `-nc $MOLCHARGE`

# create topology for ligand


