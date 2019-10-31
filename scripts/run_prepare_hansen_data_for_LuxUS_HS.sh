#!/bin/sh -l

#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=1GB

source /path/to/virtual/env/bin/activate

INPUT_FILE=wgbs_methylationlevels_chr22.txt
INPUT_FOLDER=/../LuxUS/ROC/data
OUTPUT_FOLDER=/../LuxUS_HS/results/Hansen_data/chr22/preanalysis
DESIGN_MATRIX=/../LuxUS/ROC/scripts/design_matrix_test_fullrank.txt
MEANCOVFILE="$OUTPUT_FOLDER"/mean_coverage_chr22.txt
CYTNFILE="$OUTPUT_FOLDER"/number_of_cytosines_in_window_chr22.txt

#The BSEFF values have been estimated as a mean of the two flowcell rates from Hansen et al. 2011, supplement table 14. Rounded to 4 decimals.
BSEFF='[0.9978,0.9977,0.9979,0.9975,0.9979,0.9979]'
BSBEFF='[0,0,0,0,0,0]'
SEQERR='[0,0,0,0,0,0]'
TEST_COV_IND=1
WINDOW_WIDTH=2000
REQUIRED_DIFF=0.1


python /../LuxUS_HS/scripts/prepare_data_for_LuxUS_HS.py -i "$INPUT_FOLDER"/"$INPUT_FILE" -d $DESIGN_MATRIX -x 1 -o $OUTPUT_FOLDER -r 6  -t $TEST_COV_IND -w $WINDOW_WIDTH -u $REQUIRED_DIFF -y $MEANCOVFILE -z $CYTNFILE
