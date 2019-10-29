#!/bin/bash
#SBATCH --time=00:05:00

module load python/3.7.3
virtualenv --no-download $SLURM-TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip

pip install --no-index -r requirements.txt

python Torus.py
