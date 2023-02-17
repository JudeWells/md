#!/bin/bash -l

# Request 60 minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=0:50:0

# CPU cores
#$ -pe smp 10

# GPU
#$ -l gpu=1

# Memory
#$ -l mem=1G

# Directory
#$ -l tmpfs=15G

# Job name
#$ -N MMGBSA_new_

#$ -wd /home/ucapfdf/Scratch/last_jude_set/logs

# JOB STARTS

# Load Modules
module purge
module load beta-modules                    
module load gcc-libs/10.2.0                 
module load compilers/gnu/10.2.0            
module load python/3.9.10                   
module load openblas/0.3.7-serial/gnu-4.9.2 
module load python3/3.9                     
module load python3/recommended
module load cuda/11.3.1/gnu-10.2.0
module load numactl/2.0.12
module load binutils/2.36.1/gnu-10.2.0
module load ucx/1.9.0/gnu-10.2.0
module load mpi/openmpi/4.0.5/gnu-10.2.0
module load gromacs/2021.5/cuda-11.3

# Molecule Name:
d=

# Source Directory
base=/home/ucapfdf/Scratch/last_jude_set
res=/home/ucapfdf/Scratch/last_jude_set/ready/$d
cd $TMPDIR

# Setting up working directory
cp $res/* .
cp $base/source/posre.itp .
cp -r $base/amber99sb-ildn.ff amber99sb-ildn_zn.ff

# Create new box
gmx_cuda editconf -f complex.gro -o box.gro -bt dodecahedron -c

# Solvate
gmx_cuda solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

# Add Ions
gmx_cuda grompp -f $base/source/ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo SOL | gmx_cuda genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# Load Fast Gromacs Modules

export GMX_FORCE_UPDATE_DEFAULT_GPU=true

# Energy Minimization
gmx_cuda grompp -f $base/source/em_steep.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn -1
gmx_cuda mdrun -nb gpu -deffnm em -ntmpi 1 -ntomp 10 -pin on

# NVT Equilibration
gmx_cuda grompp -f $base/source/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn -1
gmx_cuda mdrun -nb gpu -deffnm nvt -ntmpi 1 -ntomp 10 -pin on

# NPT Equilibration
gmx_cuda grompp -f $base/source/npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn -1
gmx_cuda mdrun -deffnm npt -pin on -nb gpu -bonded gpu -pme gpu -nstlist 400 -ntmpi 1 -ntomp 10

# Production MD (1 ns)
gmx_cuda grompp -f $base/source/md.mdp -c npt.gro -r npt.gro -p topol.top -o md.tpr -maxwarn -1
gmx_cuda mdrun -deffnm md -pin on -nb gpu -bonded gpu -pme gpu -nstlist 400 -ntmpi 1 -ntomp 10

# Make Index and Center TRJ for MMGBSA
echo -e "1 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 \n q" | gmx_cuda make_ndx -f md.gro -o $TMPDIR/index.ndx
echo 1 0 | gmx_cuda trjconv -f md.xtc -o md_c.xtc -s md.tpr -pbc mol -center

# Save files
resultsdir=$base/results/$d
mkdir -p $resultsdir
cp $TMPDIR/md.gro $resultsdir
cp $TMPDIR/md.tpr $resultsdir
cp $TMPDIR/md.xtc $resultsdir
cp $TMPDIR/md_c.xtc $resultsdir
cp $TMPDIR/topol.top $resultsdir
cp $TMPDIR/index.ndx $resultsdir
cp $TMPDIR/*.itp $resultsdir
cp -r $TMPDIR/amber99sb-ildn_zn.ff $resultsdir/
