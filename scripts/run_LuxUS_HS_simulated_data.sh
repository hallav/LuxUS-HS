#!/bin/sh -l 
#SBATCH --time=05:00:00 
#SBATCH --mem-per-cpu=10G 
#SBATCH --array=0-99

cd /../LuxUS_HS/scripts
module load gcc/6.3.0

source /path/to/virtual/env/bin/activate

ALGORITHM=1
PLOTS=0
OUTPUTFOLDER=/../LuxUS_HS/results/N_deviating_0/muB_minus1_4_2_3
INPUTFOLDER=/../LuxUS_HS/generated_data/N_deviating_0/muB_minus1_4_2_3

for Q in 6 12 24
do 
    for R in 6 12 24
    do
        for DIFF in 0 1
        do
            TIMEFILE=/../LuxUS_HS/results/N_deviating_0/muB_minus1_4_2_3/computation_times_ADVI_Q"$Q"_R"$R"_DIFF"$DIFF".txt
            INPUTFILE=LuxUS_HS_simulated_15102019_C10_Q"$Q"_R"$R"_W1000_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".pickle
            OUTPUTFILE=BF_LuxUS_HS_ADVI_Q"$Q"_R"$R"_DIFF"$DIFF".txt
            IDENTIFIER=Q"$Q"_R"$R"_DIFF"$DIFF"
            python run_LuxUS_HS.py -a $ALGORITHM -p $PLOTS -t $TIMEFILE -d $INPUTFILE -i $INPUTFOLDER -o $OUTPUTFOLDER -j $OUTPUTFILE -x 1 -w $SLURM_ARRAY_TASK_ID -k $IDENTIFIER
        done
    done

done

