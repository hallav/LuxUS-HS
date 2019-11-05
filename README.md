# LuxUS-HS

LuxUS-HS is a tool for differential methylation analysis, which combines binomial observation model with a GLMM with spatial correlation structure. The spatial correlation structure includes indicator variables, which tell whether a cytosine (in the genomic window of interest) follows the same spatial correlation pattern as its neighboring cytosines. The statistical testing is done for each cytosine separately. The observation model and the statistical testing method are the same as in LuxGLM [1]. The preprocessing step used for LuxUS-HS is the same as in LuxUS [5].  

## Outline
* Requirements
* Simulating data from the LuxUS-HS model
* Running LuxUS-HS analysis
* LuxUS-HS analysis workflow
* Producing results for *title*
* References

## Requirements
- Python 3.6
- NumPy
- SciPy
- Matplotlib
- PyStan [3]
- CmdStan [4] (for running ADVI)

The versions used were: NumPy 1.17.0, SciPy 1.3.0, Matplotlib 3.1.1, PyStan 2.19.0.0, CmdStan 2.20.0. The tool has been tested in Linux environment.

## Simulating data from the LuxUS-HS model
The script *generate_data_from_LuxUS_HS.py* (called by *run_generate_data_from_LuxUS_HS.sh*) can be used to generate data from LuxUS-HS model. The script supports fixed effect with intercept and case/control indicator variable. The arguments for the script are:

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
LuxUS-HS analysis can be run with script *run_LuxUS_HS.py*, which is called by script *run_LuxUS_HS.sh* (which shows example of running LuxUS-HS for a simulated data set). The script depends on the Stan model files *luxus_HS.stan* and *luxus_1cytosine.stan* (they have to be stored in the same folder as this script, and when using ADVI the model has to be compiled with CmdStan beforehand) and *savagedickey.py*. The script writes the calculated Bayes factors for each cytosine in the defined output file. The means of the samples for the methylation proportions and indicator variable **d**s will be written in files with specified file name identifier. The arguments for the script are
```
usage: run_LuxUS_HS.py [-h] [-a ALGORITHM] [-p DIAGNOSTIC_PLOTS]
                       [-g N_GRADSAMPLES] [-e N_ELBOSAMPLES]
                       [-v N_OUTPUTSAMPLES_VI] [-m N_OUTPUTSAMPLES_HMC]
                       [-c N_HMC_CHAINS] [-t TIMEFILE] -d INPUT_DATA -i
                       INPUTFOLDER -o OUTPUTFOLDER -j OUTPUTFILE -x
                       TEST_COVARIATE [-y TEST_COVARIATE2] [-w WINDOW_INDEX]
                       [-b SIGMAB2] -k FILEIDENTIFIER

Runs LuxUS HS model for the given input data and returns a BF for the whole
window.

optional arguments:
  -h, --help            show this help message and exit
  -a ALGORITHM, --algorithm ALGORITHM
                        Give value 0 (use HMC, default) or 1 (use VI).
  -p DIAGNOSTIC_PLOTS, --diagnostic_plots DIAGNOSTIC_PLOTS
                        Give value 0 (do not plot sample diagnostics for HMC)
                        or 1 (plot sample diagnostics for HMC). Default value
                        is 0.
  -g N_GRADSAMPLES, --N_gradsamples N_GRADSAMPLES
                        Number of gradient samples used in VI. Default value
                        10 is used if not specified.
  -e N_ELBOSAMPLES, --N_elbosamples N_ELBOSAMPLES
                        Number of ELBO samples used in VI. Default value 200
                        is used if not specified.
  -v N_OUTPUTSAMPLES_VI, --N_outputsamples_VI N_OUTPUTSAMPLES_VI
                        Number of posterior samples used in VI. Default value
                        2000 is used if not specified.
  -m N_OUTPUTSAMPLES_HMC, --N_outputsamples_HMC N_OUTPUTSAMPLES_HMC
                        Number of posterior samples per chain used in HMC (the
                        burn-in will be removed from this sample number).
                        Default value 1000 is used if not specified.
  -c N_HMC_CHAINS, --N_HMC_chains N_HMC_CHAINS
                        Number of chains in HMC sampling. Default value 4 is
                        used if not specified.
  -t TIMEFILE, --timeFile TIMEFILE
                        File name (and path) for storing computation time. If
                        no file name is given the computation times will not
                        be stored into a file.
  -d INPUT_DATA, --input_data INPUT_DATA
                        Name of the data file with path.
  -i INPUTFOLDER, --inputFolder INPUTFOLDER
                        Folder where the input file is stored.
  -o OUTPUTFOLDER, --outputFolder OUTPUTFOLDER
                        Folder where to store the results.
  -j OUTPUTFILE, --outputFile OUTPUTFILE
                        File into which the BFs are written. Will be located
                        in folder specified in -o.
  -x TEST_COVARIATE, --test_covariate TEST_COVARIATE
                        Covariate to be tested. Give index (in design matrix)
                        starting from 0.
  -y TEST_COVARIATE2, --test_covariate2 TEST_COVARIATE2
                        Type 2 test: the covariate to be compared to the
                        covariate defined by argument -x. If not provided,
                        type 1 test will be performed. Give index (in design
                        matrix) starting from 0.
  -w WINDOW_INDEX, --window_index WINDOW_INDEX
                        The index of the window being analysed. If value is
                        not given the BF is saved without window index into
                        the defined output file.
  -b SIGMAB2, --sigmaB2 SIGMAB2
                        Prior for sigmaB2, used in BF calculation. Default
                        value 15 is used if not specified.
  -k FILEIDENTIFIER, --fileIdentifier FILEIDENTIFIER
                        File identifier for the d and theta mean files.
```

## LuxUS-HS analysis workflow
The LuxUS-HS analysis workflow for a bisulfite sequencing data set (starting from raw sequencing data) is following:
* Align reads to a reference genome.
* Transform the aligned reads to count data.
* Estimate the experimental parameters such as bisulfite conversion efficiency if spike-in samples are available
* Run preanalysis step to get genomic windows with F-test p-value below the desired threshold (script *prepare_data_for_LuxUS_HS.py*).
* Run LuxUS-HS analysis for each chosen genomic window, this step is parallelizable (script *run_LuxUS_HS.sh*).
* Set BF threshold to determine DMR, possibly combine DMRs and filter them further.
* Investigation of estimated indicator variable **d**.

The preanalysis step and its parameters are explained in detail in https://github.com/hallav/LuxUS

## Producing results for *LuxUS-HS: spatial correlation structure with horseshoe prior enables analysis of cytosines with deviating methylation states*

The scripts for producing the results presented in *LuxUS-HS: spatial correlation structure with horseshoe prior enables analysis of cytosines with deviating methylation states* by Halla-aho and Lähdesmäki [9] are available in this repository. The scripts for running BiSeq [6] and RADMeth [7] tools for simulation comparisons and scripts for preparing the colon cancer data set [8] for LuxUS-HS analysis are also included.

## References

[1] Äijö, T., Yue, X., Rao, A., Lähdesmäki, H (2016) LuxGLM: a probabilistic covariate model for quantification of DNA methylation modifications with complex experimental designs. *Bioinformatics*, 32(17), i511-i519.

[2] Carpenter, B., Gelman, A., Hoffman,M. D., Lee, D., Goodrich, B., Betancourt, M., Brubaker, M., Guo, J., Li, P., Riddell, A. (2017) Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1).

[3] Stan Development Team (2018) CmdStan: the command-line interface to Stan, Version 2.18.0.   http://mc-stan.org

[4] Stan Development Team (2017) PyStan: the Python interface to Stan, Version 2.16.0.0.   http://mc-stan.org

[5] Halla-aho, V. and Lähdesmäki, H. (2019) LuxUS: Detecting differential DNA methylation using generalized linear mixed model with spatial correlation structure. doi: https://doi.org/10.1101/536722

[6] Hebestreit, K., Dugas, M. and Klein, H-U. (2013) Detection of significantly differentially methylated regions in targeted bisulfite sequencing data. Bioinformatics, 29(13):1647–1653

[7] Dolzhenko, E. and Smith, A. D. (2014) Using beta-binomial regression for high-precision differential methylation analysis in multifactor whole-genome bisulfite sequencing experiments. BMC bioinformatics, 15(1):215, 2014.

[8] Hansen, K. D. (2016) bsseqData: Example whole genome bisulfite data for the bsseq package. 

[9] Halla-aho, V. and Lähdesmäki, H. (2019) LuxUS-HS: spatial correlation structure with horseshoe prior enables analysis of cytosines with deviating methylation states. 
