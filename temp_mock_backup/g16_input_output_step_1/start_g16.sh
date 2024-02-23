#!/bin/bash
  
#SBATCH --mail-user=jansea92@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=1000
#SBATCH --time=0-00:10:00
#SBATCH --qos=hiprio
#SBATCH --partition=main
#SBATCH --job-name=bug_1

#Load modules and libs
module load gaussian/g16_A03

g16 testopt.1.gjf
