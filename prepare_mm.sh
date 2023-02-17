#!/bin/bash -l

# Request 60 minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=0:30:0

# CPU cores
#$ -pe smp 20

# Memory
#$ -l mem=1G

# Directory
#$ -l tmpfs=15G

# Job name
#$ -N MMGBSA_new_

#$ -wd /home/ucapfdf/Scratch/last_jude_set/logs

# JOB STARTS

# Load Modules

d=
base=/home/ucapfdf/Scratch/last_jude_set/results
res=$base/$d

cd $TMPDIR

#MMGBSA Modules
module purge
module load beta-modules
module load gcc-libs/7.3.0
module load python/3.9.10 
module load openblas/0.3.7-serial/gnu-4.9.2
module load python3/3.9        
module load python3/recommended
module load compilers/gnu/7.3.0
module load mpi/openmpi/3.1.4/gnu-7.3.0
module load gromacs/2021.2/gnu-7.3.0

#Convert TRJ
echo -e "1 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 \n q" | gmx make_ndx -f $res/md.gro -o $TMPDIR/index.ndx
#echo 1 0 | gmx trjconv -f $res/md.xtc -o $res/md_c.xtc -s $res/md.tpr -pbc mol -center

#MMGBSA 
mpirun gmx_MMPBSA MPI -O -i $base/mmgbsa.in -cs $res/md.tpr -ci $TMPDIR/index.ndx -cg 42 22 -ct $res/md_c.xtc -cp $res/topol.top -o $TMPDIR/lig_$d.dat -nogui

cp $TMPDIR/lig_$d.dat $res 
