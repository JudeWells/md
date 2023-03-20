# download the charmm forcefield and
# cgenff_charmm2gmx_py3_nx2.py from here:
# http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs

# write the topology for the protein using pdb2gmx
export PATH=/Applications/gromacs-2022/build/bin:$PATH
LIGNAME=jz4
PROTNAME=3HTB_processed
pwd
# we use the charm36 FF which must be in the local dir, tip3p is default water model
gmx pdb2gmx -f 3HTB_clean.pdb -o ${PROTNAME}.gro -ff charmm36_mar2019 -water tip3p

# In Daniel's script we can either use mol2 or sdf - pdbqt must be converted
obabel ${LIGNAME}.pdb -O ${LIGNAME}.sdf -p 7.4

# prepare the ligand with forcefields

antechamber -i ${LIGNAME}.sdf -fi sdf -o ${LIGNAME}.mol2 -fo mol2 -pf y -s 0 -j 5 -c bcc
# If command above fails consider adding following argument `-nc $MOLCHARGE`

# Start create topology #
echo "source oldff/leaprc.ff99SB
source leaprc.gaff
loadamberparams ${LIGNAME}.frcmod
lig = loadmol2 ${LIGNAME}.mol2
check lig
saveoff lig ${LIGNAME}.lib
saveamberparm lig ${LIGNAME}.prmtop ${LIGNAME}.rst7
quit" >> ligprep_${LIGNAME}.in

parmchk2 -i ${LIGNAME}.mol2 -f mol2 -o ${LIGNAME}.frcmod
tleap -s -f ligprep_${LIGNAME}.in > ${LIGNAME}_ligprep.out
# End create topology #

# Start convert to Gromacs #
# this should create .gro file and .top file
python jw_convert_gmx.py -i ${LIGNAME}.prmtop
# End convert to Gromacs #

# Create complex file
python jw_create_complex_topology.py --ligand ${LIGNAME} --protein ${PROTNAME}
# Move the following files to the working directory
# complex.gro
# topol.top
# ${LIGNAME}.itp






