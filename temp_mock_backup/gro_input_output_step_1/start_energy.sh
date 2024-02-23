#!/bin/bash

#SBATCH --mail-user=jansea92@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=1000
#SBATCH --time=0-00:10:00
#SBATCH --qos=hiprio
#SBATCH --partition=main
#SBATCH --job-name=testopt

#Load modules and libs
module load GROMACS/2019-foss-2018b
module load CUDA/9.2.88-GCC-7.3.0-2.30
source /scratch/spetry/Gromacs_bin/bin/GMXRC
export GMXLIB=/home/jansea92/GROLIB/top

gmx energy -f testopt.edr -o testopt.edr.xvg -backup no
