#!/bin/bash
  
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000
#SBATCH --time=14-00:00:00
#SBATCH --qos=standard
#SBATCH --partition=main
#SBATCH --job-name=opt_25

# enable autoswap modules
LMOD_DISABLE_SAME_NAME_AUTOSWAP=no


#Load modules and libs
#module load ORCA
#module load Serenity
module load GROMACS # GROMACS
module unload python
module load Python/3.6.6-foss-2018b
module add HDF5/1.12.1-gompi-2021b
module add imkl/2021.4.0
#source /scratch/spetry/Gromacs_bin/bin/GMXRC
source /home/nicoa96/serenity-master/serenity.sh
module load CUDA/9.2.88-GCC-7.3.0-2.30
export GMXLIB=/home/nicoa96/students_start/GROLIB/top

module add SciPy-bundle
# make temporary directory on /scratch/

export TMPDIR=/scratch/$USER/qmmm/tmp.$SLURM_JOBID
if [ -d $TMPDIR ]; then
  echo "$TMPDIR exists; double job start; exit"
  exit 1
fi

mkdir -p $TMPDIR

## Your project goes here
export PROJECT=`pwd`

set -e

# Copy Inputs to TMPDIR and run calculation
cp -r -p $PROJECT/* $TMPDIR
rm $TMPDIR/slurm-*.out
cd $TMPDIR

python /home/nicoa96/qmmm_revision/gmx2qmmm.py -c gro_150.gro -p topol.top -act active.ndx -n index.ndx -path path.dat

# remove leftovers and get back to the project directory

cp -r * $PROJECT
#cd ../
#mv -r $TMPDIR completed_temp/.


