#!/bin/bash
  
#SBATCH --mail-user=jansea92@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=1000
#SBATCH --time=0-00:10:00
#SBATCH --qos=hiprio
#SBATCH --partition=main
#SBATCH --job-name=cap_scan

#Load modules and libs
module load gaussian/g16_A03
module load GROMACS/2019-foss-2018b # GROMACS
module load Python/3.6.6-foss-2018b
source /scratch/spetry/Gromacs_bin/bin/GMXRC
module load CUDA/9.2.88-GCC-7.3.0-2.30
export GMXLIB=/home/jansea92/GROLIB/top

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

#cd /home/jansea92/qmmm_revision/example/opt
#echo 'Changed DIR'

#python /home/jansea92/qmmm_revision/gmx2qmmm.py

python ~/qmmm_revision/gmx2qmmm.py

# remove leftovers and get back to the project directory

cp -r * $PROJECT
cd ../
mv $TMPDIR completed_temp/.

