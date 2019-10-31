#!/bin/bash -l
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=200M
#SBATCH --array=0-99


module load GSL/2.4-goolf-triton-2017a

ND=0
MUB=muB_minus1_4_2_3
ID=15102019

OUTFOLDER=/../LuxUS_HS/results/N_deviating_"$ND"/"$MUB"
CHUNKFOLDER=/../LuxUS_HS/generated_data/N_deviating_"$ND"/"$MUB"
DESIGNFOLDER=/../LuxUS_HS/generated_data/N_deviating_"$ND"/"$MUB"
TIMERESULTFOLDER=/../LuxUS_HS/results/N_deviating_"$ND"/"$MUB"
TIMEF="runtimes_radmeth"
WIDTH=1000

for REPS in 6 12 24
do
  for READS in 6 12 24
  do
    for DIFF in 0 1
    do


      START=$(date +%s.%N)
      /../software/methpipe-3.4.3/bin/radmeth regression -factor IsCase "$DESIGNFOLDER"/design_matrix_"$ID"_C10_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt "$CHUNKFOLDER"/proportion_table_"$ID"_C10_Q"$READS"_R"$REPS"_W"$WIDTH"_set"$SLURM_ARRAY_TASK_ID"_diff"$DIFF".txt > "$OUTFOLDER"/radmeth_cpgs_Q"$READS"_R"$REPS"_diff"$DIFF"_"$SLURM_ARRAY_TASK_ID".bed
      /../software/methpipe-3.4.3/bin/radmeth adjust -bins 1:200:1 "$OUTFOLDER"/radmeth_cpgs_Q"$READS"_R"$REPS"_diff"$DIFF"_"$SLURM_ARRAY_TASK_ID".bed > "$OUTFOLDER"/radmeth_cpgs_Q"$READS"_R"$REPS"_diff"$DIFF"_"$SLURM_ARRAY_TASK_ID".adjusted.bed
      RUNTIME=$(echo "$(date +%s.%N) - $START" | bc)
      printf "%s\t%s\n" "$RUNTIME" "$SLURM_ARRAY_TASK_ID" >> "$TIMERESULTFOLDER"/"$TIMEF"_C10_Q"$READS"_R"$REPS"_D"$DIFF".txt

      done
    done
done


 
