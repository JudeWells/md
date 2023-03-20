# download the charmm forcefield and
# cgenff_charmm2gmx_py3_nx2.py from here:
# http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs

# This script runs the full pipeline for calculating MMGBSA
# for a protein-ligands complex
# the two files required for input are
# 1. the protein in pdb format
# 2. the ligand in pdb format
# example usage:
# ./prepare.sh 3HTB jz4
# on the cluster conda environments can be found here: /SAN/orengolab/nsp13/.conda/envs
#

## These paths are for macbook
#export PATH=/Applications/gromacs-2022/build/bin:$PATH
#source /usr/local/gromacs/bin/GMXRC
## end macbook paths

### these paths are for the cluster:
export PATH=/SAN/orengolab/nsp13/
export PROJECT_USER=shared
export PROJECT_DIR=/SAN/orengolab/nsp13
export PROJECT_HOME=${PROJECT_DIR}/${PROJECT_USER}
source $PROJECT_HOME/source_files/conda.source
conda activate gmxMMPBSA/
### end cluster paths

PROTNAME=$1
LIGNAME=$2
OUTDIR=output
pwd
# write the topology for the protein using pdb2gmx
# we use the charm36 FF which must be in the local dir, tip3p is default water model
gmx pdb2gmx -f ${PROTNAME}.pdb -o ${OUTDIR}/${PROTNAME}.gro -p ${OUTDIR}/topol.top -i ${OUTDIR}/posre.itp -ff charmm36_mar2019 -water tip3p

# In Daniel's script we can either use mol2 or sdf - pdbqt must be converted
obabel ${LIGNAME}.pdb -O ${OUTDIR}/${LIGNAME}.sdf -p 7.4

# prepare the ligand with forcefields
antechamber -i ${OUTDIR}/${LIGNAME}.sdf -fi sdf -o ${OUTDIR}/${LIGNAME}.mol2 -fo mol2 -pf y -s 0 -j 5 -c bcc
# If command above fails consider adding following argument `-nc $MOLCHARGE`

# Start create topology #
echo "source oldff/leaprc.ff99SB
source leaprc.gaff
loadamberparams ${OUTDIR}/${LIGNAME}.frcmod
lig = loadmol2 ${OUTDIR}/${LIGNAME}.mol2
check lig
saveoff lig ${OUTDIR}/${LIGNAME}.lib
saveamberparm lig ${OUTDIR}/${LIGNAME}.prmtop ${OUTDIR}/${LIGNAME}.rst7
quit" >> ${OUTDIR}/ligprep_${LIGNAME}.in

parmchk2 -i ${OUTDIR}/${LIGNAME}.mol2 -f mol2 -o ${OUTDIR}/${LIGNAME}.frcmod
tleap -s -f ${OUTDIR}/ligprep_${LIGNAME}.in > ${OUTDIR}/${LIGNAME}_ligprep.out
# End create topology #

# Start convert to Gromacs #
# this should create .gro file and .top file
python jw_convert_gmx.py -i ${OUTDIR}/${LIGNAME}.prmtop
# End convert to Gromacs #

# Create complex file
python jw_create_complex_topology.py --ligand ${LIGNAME} --protein ${PROTNAME} --output ${OUTDIR}
# Move the following files to the working directory
# complex.gro
# topol.top
# ${LIGNAME}.itp

# DF's prepare.sh script starts here
$BASE=.
cd $TMPDIR
# Create new box
gmx editconf -f ${OUTDIR}/complex.gro -o ${OUTDIR}/box.gro -bt dodecahedron -c

# Solvate
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p $topol.top

# Add Ions
gmx grompp -f $$BASE/source/ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# Load Fast Gromacs Modules

export GMX_FORCE_UPDATE_DEFAULT_GPU=true

# Energy Minimization
gmx grompp -f $BASE/source/em_steep.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn -1
gmx mdrun -nb gpu -deffnm em -ntmpi 1 -ntomp 10 -pin on

# NVT Equilibration
gmx grompp -f $BASE/source/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn -1
gmx mdrun -nb gpu -deffnm nvt -ntmpi 1 -ntomp 10 -pin on

# NPT Equilibration
gmx grompp -f $BASE/source/npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn -1
gmx mdrun -deffnm npt -pin on -nb gpu -bonded gpu -pme gpu -nstlist 400 -ntmpi 1 -ntomp 10

# Production MD (1 ns)
gmx grompp -f $BASE/source/md.mdp -c npt.gro -r npt.gro -p topol.top -o md.tpr -maxwarn -1
gmx mdrun -deffnm md -pin on -nb gpu -bonded gpu -pme gpu -nstlist 400 -ntmpi 1 -ntomp 10

# Make Index and Center TRJ for MMGBSA
echo -e "1 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 \n q" | gmx make_ndx -f md.gro -o $TMPDIR/index.ndx
echo 1 0 | gmx trjconv -f md.xtc -o md_c.xtc -s md.tpr -pbc mol -center

# Save files
RESULTSDIR=$BASE/results/$LIGNAME
mkdir -p $RESULTSDIR
cp $TMPDIR/md.gro $RESULTSDIR
cp $TMPDIR/md.tpr $RESULTSDIR
cp $TMPDIR/md.xtc $RESULTSDIR
cp $TMPDIR/md_c.xtc $RESULTSDIR
cp $TMPDIR/topol.top $RESULTSDIR
cp $TMPDIR/index.ndx $RESULTSDIR
cp $TMPDIR/*.itp $RESULTSDIR
cp -r $TMPDIR/amber99sb-ildn_zn.ff $RESULTSDIR/
# DF's prepare.sh script ends here

# DF's prepare_mm.sh script starts here
#Convert TRJ
cd $TMPDIR
echo -e "1 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 \n q" | gmx make_ndx -f $res/md.gro -o $TMPDIR/index.ndx

#MMGBSA (mpirun is message passing interface for parallelization)
#mpirun gmx_MMPBSA MPI -O -i $base/mmgbsa.in -cs $res/md.tpr -ci $TMPDIR/index.ndx -cg 42 22 -ct $res/md_c.xtc -cp $res/topol.top -o $TMPDIR/lig_$d.dat -nogui
gmx_MMPBSA MPI -O -i $base/mmgbsa.in -cs $res/md.tpr -ci $TMPDIR/index.ndx -cg 42 22 -ct $res/md_c.xtc -cp $res/topol.top -o $TMPDIR/lig_$d.dat -nogui

cp $TMPDIR/lig_$d.dat $RESULTSDIR
# DF's prepare_mm.sh script ends here



