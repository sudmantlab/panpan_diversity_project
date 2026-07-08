#!/bin/bash
#SBATCH --account=XXX
#SBATCH --partition=XXX
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source activate syri_env
set -euo pipefail # safe mode
set -x # logging

echo "Reference: $ref"
echo "Query:     $qry"
echo "PAF File:  $paf"
echo "Label:     $label"

mkdir -p $label
cd $label

syri -F P                       \
    -c $paf                     \
    -r $ref                     \
    -q $qry                     \
    --prefix $label             \
    --nc $SLURM_CPUS_PER_TASK   \
    --nosnp
