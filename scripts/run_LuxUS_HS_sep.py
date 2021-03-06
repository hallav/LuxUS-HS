import argparse
import pystan
import numpy
import matplotlib.pyplot as pyplot

import pickle
from hashlib import md5
import logging

import time
import subprocess
import csv


from savagedickey import calculate_savagedickey_kde_1d

def stan_cache(model_name, **kwargs):
    f=open(model_name, 'rb')
    model_code=f.read()
    f.close()
    code_hash = md5(model_code.decode('utf-8').encode('ascii')).hexdigest()
    cache_fn = 'cached-{}-{}.pkl'.format(model_name, code_hash)
    try:
        sm = pickle.load(open(cache_fn, 'rb'))
    except:
        sm = pystan.StanModel(file=model_name)
        with open(cache_fn, 'wb') as f:
            pickle.dump(sm, f)
    else:
        logging.info("Using cached StanModel")
    return sm.sampling(**kwargs)



def csvIntoExtractDict(csvFile,isDiagnostic):
    #Inputs:
    #csvfile: the name of the file from which data is extracted
    #isDiagnostic: logical, use when retrieving ELBO values from diagnostic files
    outputdict={}
    with open(csvFile) as cmdfile:
        reader=csv.reader(cmdfile)
        isFirstrow=0
        for row in reader:
            if row[0].startswith("#"):
                continue
            N_rowitems=len(row)
            if isFirstrow==0:
                variableNames=[]
                if isDiagnostic==1: #if the input file is a diagnostic file, the dictionary keys can't be found
                    outputdict['iter']=[] #Makes a dummy row of zeros as the first row
                    outputdict['time_in_seconds']=[]
                    outputdict['ELBO']=[]
                    variableNames=['iter','time_in_seconds','ELBO']
                else:
                    #If the file type is output file, the dictionary keys can be extracted from the file it$
                    for i in range(0,N_rowitems):
                        outputdict[row[i]]=[] #Makes a dummy row of zeros as the first row
                        variableNames.append(row[i])
                isFirstrow=1 #all dictionary keys have been inserted to the dict
            else:
                for i in range(0,N_rowitems):
                    outputdict[variableNames[i]].append(float(row[i]))

    return outputdict;


#Parsing samples of theta from extracted dictionary (from variational analysis)
def extractVariable2dim(outputdict,varName,N,M,N_samples):
    #Input arguments:
    #outputdict= the dictionary object from which we should parse the output
    #varName= the name of the variable in the dictionary, a string object
    #N= size in dimension 1, integer
    #M= size in dimension 2, integer, variable is originally a NxM matrix
    #N_samples= number of samples per variable in outputdict
    #This function is for extracting matrices (e.g. 2-dimensional arrays)
    extractVar=numpy.zeros((N_samples+1,M,N)) #Make a 3D-matrix where all samples of element (i,j) are stored into the thir$
    for i in range(0,N): #Go through all variable indices, for which the samples have been stored under separate dictionary$
        for j in range(0,M):
            t_i_j=varName+'.'+str(j+1)+'.'+str(i+1) #Forming the key
            extractVar[:,j,i]=outputdict[t_i_j] #Store the samples under the key into the matrix

    return extractVar;

#Parsing samples of theta from extracted dictionary (from variational analysis)
def extractVariable1dim(outputdict,varName,N,N_samples):
    #Input arguments:
    #outputdict= the dictionary object from which we should parse the output
    #varName= the name of the variable in the dictionary, a string object
    #N= size in dimension 1, integer, variable is originally a Nx1 vector
    #N_samples= number of samples per variable in outputdict
    #This function is for extracting vectors (e.g. 1-dimensional arrays)!!
    extractVar=numpy.zeros((N_samples+1,N)) #Make a 2D-matrix where all samples are stored
    for i in range(0,N): #Go through all variable indices, for which the samples have been stored under separate dictionary$
        if N==1:
            t_i=varName #If the original variable was a scalar, this is used as a key
        else:
            t_i=varName+'.'+str(i+1) #Forming the key
           
        extractVar[:,i]=outputdict[t_i] #Store the samples under the key into the matrix

    return extractVar;

def run_luxus_HMC(sigmaB2, lux_data, stan_model_file, plot_file_name, covariate_to_be_tested, N_outputsamples, N_chains, diagnostic_plots,testtype,test2cov,N_cytosines):

    time_start=time.time()
    #sm = StanModel_cache(stan_model_file)
    #fit = sm(data=lux_data, iter=N_outputsamples, chains=N_chains)
    fit=stan_cache(stan_model_file, data=lux_data, iter=N_outputsamples, chains=N_chains)
    samples=fit.extract(permuted=True)
    
    bf=numpy.zeros(N_cytosines)

    if testtype==1:
        for c_ind in range(0,N_cytosines):
            bf[c_ind]=calculate_savagedickey_kde_1d(numpy.zeros((1,)),sigmaB2*numpy.eye(1),samples['B'][:,2*c_ind+covariate_to_be_tested])
    else:
        for c_ind in range(0,N_cytosines):
            bf[c_ind]=calculate_savagedickey_kde_1d(numpy.zeros((1,)),sigmaB2*numpy.eye(1),samples['B'][:,2*c_ind+covariate_to_be_tested]-samples['B'][:,2*c_ind+test2cov])

    time_end_full=time.time()
    print(fit)

    if diagnostic_plots==1:
        if N_cytosines ==1:
            vars_to_be_printed=['sigmaE2','sigmaR2','B']
        else:
            vars_to_be_printed=['l','sigmaE2','sigmaR2','B','d']
        
        fit.plot(vars_to_be_printed)
        pyplot.tight_layout()
        pyplot.savefig(plot_file_name)

    sigmaR2_mean=numpy.mean(samples['sigmaR2'])
    sigmaE2_mean=numpy.mean(samples['sigmaE2'])
    theta_median=numpy.median(samples['theta'],axis=0)

    if N_cytosines>1:
        d_median=numpy.median(samples['d'],axis=0)
    else:
        d_median=numpy.array([0])

    runtime=time_end_full-time_start
	
    return bf, runtime, sigmaR2_mean, sigmaE2_mean, d_median, theta_median;	

def run_luxus_VI(lux_data, model_name, N_gradsamples, N_elbosamples, N_outputsamples,temp_input_data_file_name, temp_output_file_name, sigmaB2,diagnostic_plots,plot_file_name, N_predictors, N_reps, covariate_to_be_tested,testtype, test2cov,N_cytosines):

    time_start=time.time()
    pystan.misc.stan_rdump(lux_data,temp_input_data_file_name)
    subprocess.call("./%s variational grad_samples=%s elbo_samples=%s output_samples=%s data file=%s output file=%s"%(model_name,N_gradsamples,N_elbosamples,N_outputsamples,temp_input_data_file_name,temp_output_file_name),shell=True)
    samples=csvIntoExtractDict(temp_output_file_name,0)

    samples_B=extractVariable1dim(samples,'B',N_predictors*N_cytosines,N_outputsamples)

    bf=numpy.zeros(N_cytosines)

    if testtype==1:
        for c_ind in range(0,N_cytosines):
            bf[c_ind]=calculate_savagedickey_kde_1d(numpy.zeros((1,)),sigmaB2*numpy.eye(1),samples_B[:,2*c_ind+covariate_to_be_tested])
    else:
        for c_ind in range(0,N_cytosines):
            bf[c_ind]=calculate_savagedickey_kde_1d(numpy.zeros((1,)),sigmaB2*numpy.eye(1),samples_B[:,2*c_ind+covariate_to_be_tested]-samples_B[:,2*c_ind+test2cov])

    time_end_full=time.time()

    samples_sigmaR2=extractVariable1dim(samples,'sigmaR2',1,N_outputsamples)
    samples_sigmaE2=extractVariable1dim(samples,'sigmaE2',1,N_outputsamples)
    samples_theta=extractVariable1dim(samples,'theta',N_cytosines*N_reps,N_outputsamples)

    sigmaR2_mean=numpy.mean(samples_sigmaR2)
    sigmaE2_mean=numpy.mean(samples_sigmaE2)
    theta_median=numpy.median(samples_theta,axis=0)
    
    if N_cytosines>1:
        samples_d=extractVariable1dim(samples,'d',N_predictors*N_cytosines,N_outputsamples)
        d_median=numpy.median(samples_d,axis=0)
    else:
        d_median=numpy.array([0])

    if diagnostic_plots==1:

        pyplot.subplot(5,1,1)
        pyplot.hist(samples_sigmaR2,bins=20)
        pyplot.title(r'$\sigma_R^2$')
        pyplot.show()

        pyplot.subplot(5,1,2)
        pyplot.hist(samples_sigmaE2,bins=20)
        pyplot.title(r'$\sigma_E^2$')
        pyplot.show()

        pyplot.subplot(5,1,3)
        for b_ind in range(0,N_predictors*N_cytosines):
            pyplot.hist(samples_B[b_ind],bins=20,label="b_%s"%(b_ind))
        pyplot.legend(loc='upper right')

        if N_cytosines>1:
            samples_d=extractVariable1dim(samples,'d',N_cytosines,N_outputsamples)
            samples_l=extractVariable1dim(samples,'l',1,N_outputsamples)
            
            pyplot.subplot(5,1,4)
            pyplot.hist(samples_l,bins=20)
            pyplot.title(r'$\ell$')
            pyplot.show()

            pyplot.subplot(5,1,5)
            for c_ind in range(0,N_cytosines):
                pyplot.hist(samples_d[c_ind],bins=20,label="d_%s"%(c_ind))
                
            pyplot.legend(loc='upper right')

        pyplot.title(r'$\mathbf{b}$')
        pyplot.show()
        pyplot.savefig(plot_file_name)


    subprocess.call("rm %s"%(temp_input_data_file_name),shell=True)
    subprocess.call("rm %s"%(temp_output_file_name),shell=True)

    runtime=time_end_full-time_start

    return bf, runtime, sigmaR2_mean, sigmaE2_mean, d_median, theta_median; 



if __name__ == '__main__':


    default_sigmaB2 = 15
    default_diagnostic_plots = 0

    default_N_gradsamples = 10
    default_N_elbosamples = 200
    default_N_outputsamples_VI = 2000
    default_N_outputsamples_HMC = 1000
    default_algorithm=0
    default_N_HMC_chains=4
    
    parser = argparse.ArgumentParser(description='Runs LuxUS HS model for the given input data and returns a BF for the whole window.')

    parser.add_argument('-a','--algorithm',action='store',dest='algorithm',type=int,required=False,help='Give value 0 (use HMC, default) or 1 (use VI).')
    parser.add_argument('-p','--diagnostic_plots',action='store',dest='diagnostic_plots',type=int,required=False,help="Give value 0 (do not plot sample diagnostics for HMC) or 1 (plot sample diagnostics for HMC). Default value is %s."%(default_diagnostic_plots))
    parser.add_argument('-g','--N_gradsamples',action='store',dest='N_gradsamples',type=int,required=False,help="Number of gradient samples used in VI. Default value %s is used if not specified."%(default_N_gradsamples))
    parser.add_argument('-e','--N_elbosamples',action='store',dest='N_elbosamples',type=int,required=False,help="Number of ELBO samples used in VI. Default value %s is used if not specified."%(default_N_elbosamples))
    parser.add_argument('-v','--N_outputsamples_VI',action='store',dest='N_outputsamples_VI',type=int,required=False,help="Number of posterior samples used in VI. Default value %s is used if not specified."%(default_N_outputsamples_VI))
    parser.add_argument('-m','--N_outputsamples_HMC',action='store',dest='N_outputsamples_HMC',type=int,required=False,help="Number of posterior samples per chain used in HMC (the burn-in will be removed from this sample number). Default value %s is used if not specified."%(default_N_outputsamples_HMC))
    parser.add_argument('-c','--N_HMC_chains',action='store',dest='N_HMC_chains',type=int,required=False,help="Number of chains in HMC sampling. Default value %s is used if not specified."%(default_N_HMC_chains))
    parser.add_argument('-t','--timeFile',action='store',dest='timeFile',type=str,required=False,help='File name (and path) for storing computation time. If no file name is given the computation times will not be stored into a file.')
    parser.add_argument('-d','--input_data',action='store',dest='input_data',type=str,required=True,help='Name of the data file with path.')
    parser.add_argument('-i','--inputFolder',action='store',dest='inputFolder',type=str,required=True,help='Folder where the input file is stored.')	
    parser.add_argument('-o','--outputFolder',action='store',dest='outputFolder',type=str,required=True,help='Folder where to store the results.')	
    parser.add_argument('-j','--outputFile',action='store',dest='outputFile',type=str,required=True,help='File into which the BFs are written. Will be located in folder specified in -o.')
    parser.add_argument('-x','--test_covariate',action='store',dest='test_covariate',type=int,required=True,help='Covariate to be tested. Give index (in design matrix) starting from 0.')
    parser.add_argument('-y','--test_covariate2',action='store',dest='test_covariate2',type=int,required=False,help='Type 2 test: the covariate to be compared to the covariate defined by argument -x. If not provided, type 1 test will be performed. Give index (in design matrix) starting from 0.')
    parser.add_argument('-w','--window_index',action='store',dest='window_index',type=int,required=False,help='The index of the window being analysed. If value is not given the BF is saved without window index into the defined output file.')
    parser.add_argument('-b','--sigmaB2',action='store',dest='sigmaB2',type=float,required=False,help="Prior for sigmaB2, used in BF calculation. Default value %s is used if not specified."%(default_sigmaB2))
    parser.add_argument('-s','--sep_cytosine_index',action='store',dest='cytosine',type=int,required=False,help='If LuxUS analysis is performed separately for each cytosine (e.g. in case of analysis of simulated data), the index of the current cytosine should be given for storing the results properly. Give index starting from 0.')

    options = parser.parse_args()
    
    if  options.sigmaB2 is None:
        print("sigmaB2 was not specified. Default value %s is used."%(default_sigmaB2))
        sigmaB2=default_sigmaB2
    else:
        sigmaB2=options.sigmaB2

    if options.algorithm is None:
        print("Algorithm was not specified. Default value (HMC) is used.")
        algorithm=0
    else: 
        algorithm=options.algorithm

    if  options.test_covariate2 is None:
        test_type=1
        test_type2_cov=0
    else:
        test_type=2
        test_type2_cov=options.test_covariate2
    


    inputf=open("%s/%s"%(options.inputFolder,options.input_data),'rb')
    luxus_data=pickle.load(inputf)
    inputf.close()
    N_covariates=luxus_data['n_predictors']
    N_replicates=luxus_data['n_replicates']

    currenttime=time.localtime()
    
    input_data_id=options.input_data.split('.')[0]    
    
    if options.cytosine is None:
        d_median_file="d_median_algorithm%s_%s"%(algorithm,input_data_id)
        theta_median_file="theta_median_algorithm%s_%s"%(algorithm,input_data_id)
    else:
        d_median_file="d_median_sep_algorithm%s_%s"%(algorithm,input_data_id)#Although this is never saved..
        theta_median_file="theta_median_sep_algorithm%s_%s"%(algorithm,input_data_id)    
    
    if algorithm==0:

        if  options.N_outputsamples_HMC is None:
            print("N_outputsamples_HMC was not specified. Default value %s is used."%(default_N_outputsamples_HMC))
            N_outputsamples_HMC=default_N_outputsamples_HMC
        else:
            N_outputsamples_HMC=options.N_outputsamples_HMC

        if  options.N_HMC_chains is None:
            print("N_HMC_chains was not specified. Default value %s is used."%(default_N_HMC_chains))
            N_HMC_chains=default_N_HMC_chains
        else:
            N_HMC_chains=options.N_HMC_chains


        if  options.diagnostic_plots is None:
            print("diagnostic_plots was not specified. Default value %s is used."%(default_diagnostic_plots))
            diagnostic_plots=default_diagnostic_plots
        else:
            diagnostic_plots=options.diagnostic_plots

        plot_file_name="%s/%s_diagnostic_plots_HMC.png"%(options.outputFolder,input_data_id)

        if luxus_data['n_cytosines']>1:
            stan_file='luxus_HS.stan'
        else:
            stan_file='luxus_1cytosine.stan'

        print("Estimating the model parameters with HMC.")
        BF, runtime, sigmaR2_mean, sigmaE2_mean, d_median, theta_median = run_luxus_HMC(sigmaB2, luxus_data, stan_file, plot_file_name, options.test_covariate, N_outputsamples_HMC,N_HMC_chains,diagnostic_plots,test_type,test_type2_cov,luxus_data['n_cytosines'])

    if options.algorithm==1:

        if  options.N_gradsamples is None:
            print("N_gradsamples was not specified. Default value %s is used."%(default_N_gradsamples))
            N_gradsamples=default_N_gradsamples
        else:
            N_gradsamples=options.N_gradsamples

        if  options.N_elbosamples is None:
            print("N_elbosamples was not specified. Default value %s is used."%(default_N_elbosamples))
            N_elbosamples=default_N_elbosamples
        else:
            N_elbosamples=options.N_elbosamples

        if  options.N_outputsamples_VI is None:
            print("N_outputsamples_VI was not specified. Default value %s is used."%(default_N_outputsamples_VI))
            N_outputsamples_VI=default_N_outputsamples_VI
        else:
            N_outputsamples_VI=options.N_outputsamples_VI

        if  options.diagnostic_plots is None:
            print("diagnostic_plots was not specified. Default value %s is used."%(default_diagnostic_plots))
            diagnostic_plots=default_diagnostic_plots
        else:
            diagnostic_plots=options.diagnostic_plots

        if luxus_data['n_cytosines']>1:
            stan_file="luxus_HS"
        else:
       	    stan_file="luxus_1cytosine"

        plot_file_name="%s/%s_diagnostic_plots_ADVI.png"%(options.outputFolder,input_data_id)


        print("Estimating the model parameters with VI.")

        variational_temp_data_file="%s/TEMP_store_variational_input_%s_%s_%s_%s_%s_%s_%s.txt"%(options.outputFolder,input_data_id,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])
        variational_temp_output_file="%s/TEMP_store_variational_results_%s_%s_%s_%s_%s_%s_%s.csv"%(options.outputFolder,input_data_id,currenttime[0],currenttime[1],currenttime[2],currenttime[3],currenttime[4],currenttime[5])

        BF, runtime, sigmaR2_mean, sigmaE2_mean, d_median, theta_median = run_luxus_VI(luxus_data,stan_file, N_gradsamples, N_elbosamples, N_outputsamples_VI, variational_temp_data_file, variational_temp_output_file, sigmaB2,diagnostic_plots,plot_file_name,N_covariates, N_replicates, options.test_covariate,test_type,test_type2_cov,luxus_data['n_cytosines'])

    print("Calculated Bayes factor is %s"%(BF))
    print("Bayes factor calculation took %s seconds."%(runtime))
    
    if options.window_index is None:
        with open("%s/%s"%(options.outputFolder,options.outputFile),'a+') as fw:
            if options.cytosine is None:
                for c_ind in range(0,luxus_data['n_cytosines']):
                    fw.write("%f\t%s\n"%(BF[c_ind],c_ind))
            else:
                fw.write("%f\t%s\n"%(BF,options.cytosine))

    else:
        with open("%s/%s"%(options.outputFolder,options.outputFile),'a+') as fw:
            if options.cytosine is None:
                for c_ind in range(0,luxus_data['n_cytosines']):
                    fw.write("%f\t%s\t%s\n"%(BF[c_ind],options.window_index,c_ind))
            else:
                fw.write("%f\t%s\t%s\n"%(BF,options.window_index,options.cytosine))
                
    if  options.timeFile is None:
        print("File for storing computation time was not specified and computation time will not be stored.")
    else:
        with open(options.timeFile,'a+') as tf:
            if options.cytosine is None:
                tf.write("%f\t%s\n"%(runtime,options.window_index))
            else:
                tf.write("%f\t%s\t%s\n"%(runtime,options.window_index,options.cytosine))
    
    #SAVING THE THETA SAMPLE MEDIANS    
                
    with open("%s/%s"%(options.outputFolder,theta_median_file),'a+') as fw:
        for t_ind in range(0,N_replicates*luxus_data['n_cytosines']):
            if options.cytosine is None:
                fw.write("%f\t%s\t%s\n"%(theta_median[t_ind],options.window_index,t_ind))
            else:
                fw.write("%f\t%s\t%s\t%s\n"%(theta_median[t_ind],options.window_index,t_ind,options.cytosine))
    
