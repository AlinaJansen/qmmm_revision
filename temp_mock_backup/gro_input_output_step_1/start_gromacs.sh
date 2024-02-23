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
module load GROMACS/2019-foss-2018b
source /scratch/spetry/Gromacs_bin/bin/GMXRC
export GMXLIB=/home/jansea92/GROLIB/top


gmx grompp -p testopt.qmmm.top -c testopt.boxlarge.g96 -n testopt.qmmm.top.ndx -f testopt.mdp -o testopt.tpr -backup no
