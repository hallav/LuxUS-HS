import numpy
import pickle
import argparse
#from matplotlib import pyplot
import json



#Sigmoid function
def sigmoid(x):
    if x.shape==():
        out=1.0/(1.0+numpy.exp(-x))
    else:
        out = numpy.zeros([x.shape[0]])
        for i in range(0,x.shape[0]):
            out[i]=1.0/(1.0+numpy.exp(-x[i]))
    return out


def generate_coordinates(Wrange,cytosines,predictors):
    #Generate random coordinates
    coords=numpy.zeros([1,cytosines])    
    while numpy.unique(coords).shape!=(cytosines,):
        coords = numpy.random.randint(1,Wrange,cytosines)
    
    coords = numpy.sort(coords)
    #Replicate the coordinates for later use
    rep_coords = numpy.kron(coords,numpy.ones(predictors))
    
    return rep_coords


def generate_D(n_predictors,n_replicates):
    #n_replicates should be an even number, n_predictors is set to 2 in this script
    
    n_replicates_half = n_replicates/2

    D=numpy.zeros((n_replicates,n_predictors))
    ind=0
    for i in range(0,numpy.int(n_replicates_half)):
        for condition in [0,1]:
            if condition==1:
                D_row=[1,1]
            else:
                D_row=[1,0]
            D[ind,:]=D_row
            ind+=1

    return D


def generate_u_R(N_replicates,sigmar2):

    u_R=numpy.random.normal(loc=0,scale=numpy.sqrt(sigmar2),size=N_replicates)

    return u_R



def generate_indicator_vector(n_cytosines,n_replicates):

    ind_rep=numpy.zeros(n_replicates*n_cytosines)

    for i in range(0,n_cytosines):
        ind_rep[range(i*n_replicates,(i+1)*n_replicates)]=range(0,n_replicates)

    return ind_rep


def generate_Y(X,B,Z_R,u_R,sigmae2,n_replicates,n_cytosines):

    #In this implementation, Z_R is indicator vector, not actual design matrice

    #Generate Y
    Y=numpy.zeros(n_cytosines*n_replicates)

    for i in range(0,n_replicates*n_cytosines):
        Y[i]=numpy.random.normal(numpy.dot(X[i,:],B)+u_R[Z_R[i]],sigmae2)

    return Y

def generate_counts(n_reads,Y,n_replicates,n_cytosines,SE,BSE,BSBE):

    theta=numpy.zeros(n_replicates*n_cytosines)
    counts=numpy.zeros(n_cytosines*n_replicates)
    for i in range(0,n_cytosines*n_replicates):
        theta[i]=sigmoid(Y[i])
        p = theta[i]*((1-SE[i])*(1-BSE[i])+SE[i]*BSE[i])+(1-theta[i])*((1-SE[i])*(1-BSBE[i])+SE[i]*BSBE[i])
        counts[i] = numpy.random.binomial(n_reads,p)

    return counts,theta




def generate_B(cytosines,predictors,replicates,sigmaB2,coords,differential,notCorrelated,l,mean_B,deviating_mean_B):
    #Currently this function only supports two covariates (assumed to be constant and case/control variable)

    #If notCorrelated=/=0 some of the d's elements will become 0 and the correlation between the cytosine(s) with that index will be zero
    d = numpy.ones(cytosines)
    #toBeRemoved = numpy.random.randint(0,cytosines-1,notCorrelated) THIS VERSION HAD DRAWS WITH REPLACEMENT!!
    toBeRemoved = numpy.random.choice(a=cytosines,size=notCorrelated,replace=False)
    d[toBeRemoved] = 0
    #Replicate the vector so that it can be easily accessed below
    d_kron = numpy.kron(d,numpy.ones(predictors)) ##IS THIS CORRECT??

    #Covariance matrix for B
    SigmaB = numpy.zeros([predictors*cytosines,predictors*cytosines])

    for i in range(0,predictors*cytosines):
        for j in range(0,predictors*cytosines):

            if i==j:
                SigmaB[i,j] = sigmaB2

            else:
                if (i%predictors == j%predictors): 
                    SigmaB[i,j] = sigmaB2*numpy.exp(-abs(coords[i]-coords[j])/pow(l,2))*d_kron[i]*d_kron[j] 
                    ##IS THE d_kron VECTOR CORRECTLY INDEXED? 

    if differential==1:
        muB = mean_B.copy()
    else:
        muB = mean_B.copy()
        muB[1] = 0 

    muB_kron = numpy.kron(numpy.ones((cytosines)),muB)

    if notCorrelated>0:
        for k in range(0,cytosines):
            if d[k]==0:
                if differential==1:
                    muB_kron[(k*predictors):((k+1)*predictors)] = [deviating_mean_B[0],0]
                else:
                    muB_kron[(k*predictors):((k+1)*predictors)] = deviating_mean_B

    B = numpy.random.multivariate_normal(muB_kron,SigmaB)

    #make an indicator vector that indicates whether the cytosines are differentially methylated or not
    if differential==1:
        differential_cytosines = numpy.ones(cytosines)
        differential_cytosines[toBeRemoved] = 0
    else:
        differential_cytosines = numpy.zeros(cytosines)
        differential_cytosines[toBeRemoved] = 1
    
    return B,d,differential_cytosines


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Simulates bisulphite sequencing data from the LuxUS-HS model.')

    parser.add_argument('-q','--reads',action='store',dest='reads',type=int,required=True,help='Number of BS calls i.e. (total) reads.')
    parser.add_argument('-r','--replicates',action='store',dest='replicates',type=int,required=True,help='The total number of replicates for cases and controls. 1/2 of the replicates will be cases and 1/2 will be controls.')
    parser.add_argument('-c','--cytosines',action='store',dest='N_cytosines',type=int,required=True,help='Number of cytosines.')
    parser.add_argument('-l','--lengthscale',action='store',dest='l',type=int,required=True,help='Lengthscale parameter for covariance matrix for the simulations.')
    parser.add_argument('-w','--width',action='store',dest='width',type=int,required=True,help='Width of the genomic region where the N_cytosines lie (in basepairs).')
    parser.add_argument('-n','--datasets',action='store',dest='N_datasets',type=int,required=True,help='Number of datasets to be created for both cases (no differential methylation, differential methylation). Total number of resulting datasets is 2xN_datasets.')
    parser.add_argument('-f','--folder',action='store',dest='folder',type=str,required=True,help='Folder where to store the generated data files. Format /level1/level2')
    parser.add_argument('-g','--ID',action='store',dest='ID',type=str,required=True,help='Identifier used in the resulting file names.')
    parser.add_argument('-e','--sigmae2',action='store',dest='sigmae2',type=float,required=True,help='Variance for residual term e.')
    parser.add_argument('-v','--sigmaB2',action='store',dest='sigmaB2',type=float,required=True,help='Variance for B, which is inserted into covariance matrix diagonal (SigmaB). Used for simulations.')
    parser.add_argument('-d','--sigmar2',action='store',dest='sigmar2',type=float,required=True,help='sigma_r2 used for simulations.')
    parser.add_argument('-m','--mean_B',action='store',dest='mean_B',type=str,required=True,help='Mean for B0 (coefficient for the intercept term) and value for B1 (coefficient for case/control covariate) for the cases with differential methylation. Should be given as string in format [B1, B2]')
    parser.add_argument('-b','--mean_deviating_B',action='store',dest='mean_deviating_B',type=str,required=True,help='Mean for B0 (coefficient for the intercept term) and value for B1 (coefficient for case/control covariate) for the cases with differential methylation FOR THE DEVIATING CYTOSINES. Should be given as string in format [B1, B2]')
    parser.add_argument('-j','--sigmaB2_stan',action='store',dest='sigmaB2_stan',type=float,required=True,help='Variance for B, the value which will be used in LuxUS analysis. (Can be different from what is used for simulations).')
    parser.add_argument('-x','--save_LuxUS',action='store',dest='save_LuxUS',type=int,required=True,help='0 or 1 indicating whether the seimulated data should be saved in a format supported by LuxUS.')
    parser.add_argument('-y','--save_proportion_table',action='store',dest='save_proportion_table',type=int,required=True,help='0 or 1 indicating whether the seimulated data should be saved in proportion table format, which can be used with eg. Radmeth. Also a design matrix is saved into a text file.')
    parser.add_argument('-z','--save_LuxUS_sep',action='store',dest='save_LuxUS_sep',type=int,required=True,help='0 or 1 indicating whether the seimulated data should be saved in a format supported by LuxUS analysis where each cytosine is analysed separately.')
    parser.add_argument('-k','--N_deviating',action='store',dest='N_deviating',type=int,required=True,help='Number of cytosines with deviating methylation state. Should be an integer with >=0 value.')

    options = parser.parse_args()


    mean_B = numpy.array(json.loads(options.mean_B))
    mean_deviating_B = numpy.array(json.loads(options.mean_deviating_B))

    N_reads=options.reads
    N_replicates=options.replicates
    N_cytosines=options.N_cytosines
    N_predictors=2 #Number of predictors is fixed in this script
    N_deviating=options.N_deviating
    l=options.l
    width=options.width

    N_datasets=options.N_datasets
    outputFolder=options.folder
    
    #variances and standard deviations
    sigmaB2 = options.sigmaB2
    sigmar2 = options.sigmar2
    sigmae2 = options.sigmae2


    fileID="%s_C%s_Q%s_R%s_W%s"%(options.ID,N_cytosines,N_reads,N_replicates,width)

    #Prior parameters for variance terms (to be used in the analysis with LuxUS-HS), not used in simulations!
    #Chosen based on real data analysis, see LuxUS tool paper by Halla-aho and Laehdesmaeki (2019).

    a_C=6
    b_C=3

    a_E=5
    b_E=5

    a_R=98
    b_R=143



    for i in range(0,N_datasets):
        print("Dataset %s"%(i)) 

        #Simulating experimental parameters
        bsEff_true=numpy.random.beta(99,1,size=N_replicates)
        bsbEff_true=numpy.random.beta(1,999,size=N_replicates)
        seqErr_true=numpy.random.beta(1,999,size=N_replicates)

        bsEff_true_rep = numpy.kron(numpy.ones((N_cytosines)),bsEff_true)
        bsbEff_true_rep = numpy.kron(numpy.ones((N_cytosines)),bsbEff_true)
        seqErr_true_rep = numpy.kron(numpy.ones((N_cytosines)),seqErr_true)

        bsbEff_mean=1.0/(1.0+999.0)
        seqErr_mean=1.0/(1.0+999.0)

        simulate_conversion_efficiency=numpy.random.binomial(10,numpy.tile(bsEff_true,(3066,1)))
        bsEff_estimate=numpy.sum(simulate_conversion_efficiency,axis=0)/(10.0*3066)

        print("The simulated values for (bsEff,bsbEff,seqErr)=(%s,%s,%s)"%(bsEff_true,bsbEff_true,seqErr_true))

        bsEff_rep = numpy.kron(numpy.ones((N_cytosines)),bsEff_estimate)
        bsbEff_rep = numpy.kron(numpy.ones((N_cytosines*N_replicates)),bsbEff_mean)
        seqErr_rep = numpy.kron(numpy.ones((N_cytosines*N_replicates)),seqErr_mean)
        reads_rep = numpy.kron(numpy.ones((N_cytosines*N_replicates)),N_reads)

        bsEff_rep_sep = bsEff_estimate.copy()
        bsbEff_rep_sep = numpy.kron(numpy.ones((N_replicates)),bsbEff_mean)
        seqErr_rep_sep = numpy.kron(numpy.ones((N_replicates)),seqErr_mean)
        reads_rep_sep = numpy.kron(numpy.ones((N_replicates)),N_reads)

        #Generate coordinates    
        coordinates = generate_coordinates(width,N_cytosines,N_predictors)
        
        for diff in [0,1]:
            
            print("Differential=%s"%(diff))
            
            B,d,differential_cytosines = generate_B(N_cytosines,N_predictors,N_replicates,sigmaB2,coordinates,diff,N_deviating,l,mean_B,mean_deviating_B)
            D = generate_D(N_predictors,N_replicates)
            Z_R = generate_indicator_vector(N_cytosines,N_replicates)
            u_R = generate_u_R(N_replicates,sigmar2)
            Y = generate_Y(numpy.kron(numpy.eye(N_cytosines),D),B,Z_R.astype(int),u_R,sigmae2,N_replicates,N_cytosines)
            counts,thetas = generate_counts(N_reads,Y,N_replicates,N_cytosines,seqErr_true_rep,bsEff_true_rep,bsbEff_true_rep)

            print("Differential cytosines (deviating cytosines have opposite state than the majority of cytosines):")
            print(differential_cytosines)
            print("Generated B")
            print(B)
            print("Generated replicate effect")
            print(u_R)

            #Save input for Stan when all cytosines are run together
            #Notice that the gamma distribution has different parameterization in Python and in Stan

            ind_R_p1=Z_R+1

            luxus_HS_data={'n_cytosines': N_cytosines,'n_replicates': N_replicates,'n_predictors':N_predictors,
                'bsEff': bsEff_rep, 'bsBEff': bsbEff_rep, 'seqErr': seqErr_rep,
                'bsC': counts.astype(int), 'bsTot': reads_rep.astype(int),
                'X': numpy.kron(numpy.eye(N_cytosines),D).astype(int), 'Z_R': ind_R_p1.astype(int), 'alpha': a_E, 'beta': b_E,
                'alphaR': a_R, 'betaR': b_R, 'alphal': 38, 'betal': 1, 'sigmaB2': options.sigmaB2_stan,'coordinates': coordinates}

     
            if options.save_LuxUS==1:
                output10=open("%s/LuxUS_HS_simulated_%s_set%s_diff%s.pickle"%(outputFolder,fileID,i,diff),'ab+')
                pickle.dump(luxus_HS_data,output10)
                output10.close()
                
            with open("%s/differential_cytosines_%s_set%s_diff%s.txt"%(outputFolder,fileID,i,diff),'a+') as dc:
                for k in range(0,N_cytosines):
                    dc.write("%s\n"%(differential_cytosines[k].astype(int)))



            #Form proportion table and design matrix
            

            #The proportion table header preparation
            label1='control'
            label2='case'
            labels=list()
            
            indl1=1
            indl2=1


            for k1 in range(0,N_replicates):
                if D[k1,1]==1:
                    labels.append("%s%d"%(label2,indl1))
                    indl1+=1
                else:
                    labels.append("%s%d"%(label1,indl2))
                    indl2+=1

            #Add a newline in the end of the list     
            labels.append('\n')
            

            coordinates_pt=coordinates[::N_predictors]

            if options.save_proportion_table==1:
                with open("%s/proportion_table_%s_set%s_diff%s.txt"%(outputFolder,fileID,i,diff),'a+') as pt:

                    #Write header line
                    pt.write('\t'.join(map(str,labels)))
                  
                    #Write data for each cytosine line by line
                    for k in range(0,N_cytosines):

                        pt_row=list()
                        pt_row.append("chr1:%d:+:CpG"%(coordinates_pt[k]))
  
                        for k1 in range(0,N_replicates):
                            pt_row.append(N_reads)
                            pt_row.append(counts[k*N_replicates+k1].astype(int))

                        pt_row.append('\n')
                        pt.write('\t'.join(str(v) for v in pt_row))            
                        

                with open("%s/design_matrix_%s_set%s_diff%s.txt"%(outputFolder,fileID,i,diff),'a+') as dm:
  
                    #Write header
                    dm.write("%s\t%s\n"%("Intercept","IsCase"))

                    #Write experimental design for each replicate separately
                    for k in range(0,N_replicates):
                        dm.write("%s\t%d\t%d\n"%(str(labels[k]),1,D[k,1].astype(int)))

            if options.save_LuxUS_sep==1:
                for k in range(0,N_cytosines):
    
                    counts_sep=counts[range(k*N_replicates,(k+1)*N_replicates)]

                    luxsc_data_sep = {'n_cytosines': 1,'n_replicates': N_replicates,'n_predictors':N_predictors,
                        'bsEff': bsEff_rep_sep, 'bsBEff': bsbEff_rep_sep, 'seqErr': seqErr_rep_sep,
                        'bsC': counts_sep.astype(int), 'bsTot': reads_rep_sep.astype(int),
                        'X': D, 'alpha': a_E, 'beta': b_E,
                        'sigmaB2': options.sigmaB2_stan, 'alphaR': a_R, 'betaR': b_R}

                    #Saving the generated datasets
                    output2=open("%s/LuxUS_simulated_sep_%s_cyt%s_set%s_diff%s.pickle"%(outputFolder,fileID,k,i,diff),'ab+')
                    pickle.dump(luxsc_data_sep,output2)
                    output2.close()
