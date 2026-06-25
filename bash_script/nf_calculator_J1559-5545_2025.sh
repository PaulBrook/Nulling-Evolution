#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=J1559_nf_calc
#SBATCH --account=vecchioa-gw-pta
#SBATCH --qos=bbdefault
#SBATCH --time=03-00:00:00
##SBATCH --array=1
#SBATCH --open-mode=truncate
#SBATCH --output=/rds/projects/v/vecchioa-gw-pta/brookp/clean/nulling_evolution/slurm_output/test250626.o

module purge
module load bluebear
module load Miniconda3/4.9.2

eval "$(${EBROOTMINICONDA3}/bin/conda shell.bash hook)"

#conda activate /rds/homes/b/brookp/.conda/envs/enterprise
conda activate IPTA_Env

cd /rds/projects/v/vecchioa-gw-pta/brookp/clean/nulling_evolution/clean_nf_evolution/python_code

python -u /rds/projects/v/vecchioa-gw-pta/brookp/clean/nulling_evolution/clean_nf_evolution/python_code/nf_calculator.py -pulsar J1559-5545 -fb 114 -lb 139 -outdir /rds/projects/v/vecchioa-gw-pta/brookp/clean/nulling_evolution/clean_nf_evolution/pulsar_results/ -datadir /rds/projects/v/vecchioa-gw-pta/brookp/clean/nulling_evolution/clean_nf_evolution/pulsar_data/ -bins 512
