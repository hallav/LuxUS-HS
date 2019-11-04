# LuxUS-HS

LuxUS-HS is a tool for differential methylation analysis, which combines binomial observation model with a GLMM with spatial correlation structure. The spatial correlation structure 

## Outline
* Requirements
* Simulating data from the LuxUS-HS model
* Running LuxUS-HS analysis
* LuxUS-HS analysis workflow
* References

## Simulating data from the LuxUS-HS model
The script *generate_data_from_LuxUS_HS.py* (called by *run_generate_data_from_LuxUS_HS.sh*) can be used to generate data from LuxUS-HS model. The argument for the script are:

```
usage: generate_data_from_LuxUS_HS.py [-h] -q READS -r REPLICATES -c
                                      N_CYTOSINES -l L -w WIDTH -n N_DATASETS
                                      -f FOLDER -g ID -e SIGMAE2 -v SIGMAB2 -d
                                      SIGMAR2 -m MEAN_B -b MEAN_DEVIATING_B -j
                                      SIGMAB2_STAN -x SAVE_LUXUS -y
                                      SAVE_PROPORTION_TABLE -z SAVE_LUXUS_SEP
                                      -k N_DEVIATING

Simulates bisulphite sequencing data from the LuxUS model.

optional arguments:
  -h, --help            show this help message and exit
  -q READS, --reads READS
                        Number of BS calls i.e. (total) reads.
  -r REPLICATES, --replicates REPLICATES
                        The total number of replicates for cases and controls.
                        1/2 of the replicates will be cases and 1/2 will be
                        controls.
  -c N_CYTOSINES, --cytosines N_CYTOSINES
                        Number of cytosines.
  -l L, --lengthscale L
                        Lengthscale parameter for covariance matrix for the
                        simulations.
  -w WIDTH, --width WIDTH
                        Width of the genomic region where the N_cytosines lie
                        (in basepairs).
  -n N_DATASETS, --datasets N_DATASETS
                        Number of datasets to be created for both cases (no
                        differential methylation, differential methylation).
                        Total number of resulting datasets is 2xN_datasets.
  -f FOLDER, --folder FOLDER
                        Folder where to store the generated data files. Format
                        /level1/level2
  -g ID, --ID ID        Identifier used in the resulting file names.
  -e SIGMAE2, --sigmae2 SIGMAE2
                        Variance for residual term e.
  -v SIGMAB2, --sigmaB2 SIGMAB2
                        Variance for B, which is inserted into covariance
                        matrix diagonal (SigmaB). Used for simulations.
  -d SIGMAR2, --sigmar2 SIGMAR2
                        sigma_r2 used for simulations.
  -m MEAN_B, --mean_B MEAN_B
                        Mean for B0 (coefficient for the intercept term) and
                        value for B1 (coefficient for case/control covariate)
                        for the cases with differential methylation. Should be
                        given as string in format [B1, B2]
  -b MEAN_DEVIATING_B, --mean_deviating_B MEAN_DEVIATING_B
                        Mean for B0 (coefficient for the intercept term) and
                        value for B1 (coefficient for case/control covariate)
                        for the cases with differential methylation FOR THE
                        DEVIATING CYTOSINES. Should be given as string in
                        format [B1, B2]
  -j SIGMAB2_STAN, --sigmaB2_stan SIGMAB2_STAN
                        Variance for B, the value which will be used in LuxUS
                        analysis. (Can be different from what is used for
                        simulations).
  -x SAVE_LUXUS, --save_LuxUS SAVE_LUXUS
                        0 or 1 indicating whether the seimulated data should
                        be saved in a format supported by LuxUS.
  -y SAVE_PROPORTION_TABLE, --save_proportion_table SAVE_PROPORTION_TABLE
                        0 or 1 indicating whether the seimulated data should
                        be saved in proportion table format, which can be used
                        with eg. Radmeth. Also a design matrix is saved into a
                        text file.
  -z SAVE_LUXUS_SEP, --save_LuxUS_sep SAVE_LUXUS_SEP
                        0 or 1 indicating whether the seimulated data should
                        be saved in a format supported by LuxUS analysis where
                        each cytosine is analysed separately.
  -k N_DEVIATING, --N_deviating N_DEVIATING
                        Number of cytosines with deviating methylation state.

```

## Running LuxUS-HS analysis

## LuxUS-HS analysis workflow
The LuxUS-HS analysis workflow for a bisulfite sequencing data set (starting from raw data) is following:
* Align reads to a reference genome.
* Transform the aligned reads to count data.
* Estimate the experimental parameters such as bisulfite conversion efficiency if spike-in samples are available
* Run preanalysis step to get genomic windows with F-test p-value below the desired threshold (script *prepare_data_for_LuxUS_HS.py*).
* Run LuxUS-HS analysis for each chosen genomic window, this step is parallelizable (script ).
* Set BF threshold to determine DMR, possibly combine DMRs and filter them further.
* Investigation of estimated indicator variable d.

The preanalysis step and its parameters are explained in detail in https://github.com/hallav/LuxUS

## Requirements
- Python 3.6
- NumPy
- SciPy
- Matplotlib
- PyStan [3]
- CmdStan [4] (for running ADVI)

The versions used were: NumPy 1.14.5, SciPy 1.1.0, Matplotlib 2.2.2, PyStan 2.17.1.0, CmdStan 2.12.0. The tool has been tested in Linux environment.

## References

[1] Äijö, T., Yue, X., Rao, A., Lähdesmäki, H (2016) LuxGLM: a probabilistic covariate model for quantification of DNA methylation modifications with complex experimental designs. *Bioinformatics*, 32(17), i511-i519.

[2] Carpenter, B., Gelman, A., Hoffman,M. D., Lee, D., Goodrich, B., Betancourt, M., Brubaker, M., Guo, J., Li, P., Riddell, A. (2017) Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1).

[3] Stan Development Team (2018) CmdStan: the command-line interface to Stan, Version 2.18.0.   http://mc-stan.org

[4] Stan Development Team (2017) PyStan: the Python interface to Stan, Version 2.16.0.0.   http://mc-stan.org

[5] Halla-aho, V. and Lähdesmäki, H. (2019) LuxUS: Detecting differential DNA methylation using generalized linear mixed model with spatial correlation structure. doi: https://doi.org/10.1101/536722
