#!/bin/sh -l 
#SBATCH --time=05:00:00 
#SBATCH --mem-per-cpu=10G 
#SBATCH --array=0-99

module load gcc/6.3.0

cd /../LuxUS_HS/scripts

source /path/tp/virtual/env/bin/activate

ALGORITHM=0
PLOTS=0
OUTPUTFOLDER=/../LuxUS_HS/results/N_deviating_0/muB_minus1_4_2_3
INPUTFOLDER=/../LuxUS_HS/generated_data/N_deviating_0/muB_minus1_4_2_3

for Q in 6 12 24 
do
    for R in 6 12 24
    do
        for DIFF in 0 1
        do
            for C in {0..9}
            do
                TIMEFILE=/../LuxUS_HS/results/N_deviating_0/muB_minus1_4_2_3/computation_times_sep_HMC_Q"$Q"_R"$R"_DIFF"$DIFF".txt
                INPUTFILE=LuxUS_simulated_sep_15102019_C10_Q"$Q"_R"$R"_W1000_cyt"$C"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".pickle
                OUTPUTFILE=BF_LuxUS_sep_HMC_Q"$Q"_R"$R"_DIFF"$DIFF".txt
                python run_LuxUS_HS_sep.py -a $ALGORITHM -p $PLOTS -t $TIMEFILE -d $INPUTFILE -i $INPUTFOLDER -o $OUTPUTFOLDER -j $OUTPUTFILE -x 1 -w $SLURM_ARRAY_TASK_ID -s $C
            done
        done
    done
done
