#!/bin/bash
#SBATCH -p hernquist
#SBATCH -J fdm
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -o OUTPUT.lsf
#SBATCH -e ERROR.lsf
#SBATCH -t 8640 # 4 day in min
#SBATCH --mail-user=philip.mocz@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=4000

module purge
module load matlab

srun -n 1 -c 8 matlab -nojvm -nodisplay -nosplash -nodesktop -r "axionDarkMatter_SI"

#EOF
