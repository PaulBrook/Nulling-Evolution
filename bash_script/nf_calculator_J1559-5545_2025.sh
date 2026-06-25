#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=J1559_nf_calc
#SBATCH --account=vecchioa-gw-pta
#SBATCH --qos=bbdefault
#SBATCH --time=03-00:00:00
#SBATCH --open-mode=truncate
#SBATCH --output=%x_%j.o

module purge
module load bluebear
module load Miniconda3/4.9.2

eval "$(${EBROOTMINICONDA3}/bin/conda shell.bash hook)"
conda activate IPTA_Env

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"

python -u "${BASE_DIR}/python_code/nf_calculator.py" \
    -pulsar J1559-5545 \
    -fb 114 \
    -lb 139 \
    -outdir "${BASE_DIR}/pulsar_results/" \
    -datadir "${BASE_DIR}/pulsar_data/" \
    -bins 512
