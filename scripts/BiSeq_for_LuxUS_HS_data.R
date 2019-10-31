#This script is for running BiSeq tool (with default settings) for simulated data sets 
#(simulated from LuxUS HS model) at the same time. 

#Setting up the environment 

.libPaths("/path/to/R")
Sys.setenv(R_LIBS="/path/to/R")

#Arguments to the script

input_args <- commandArgs(TRUE)

#The input arguments should be 
#1: The input data file name identifier (with path!)
#2: The number of simulated data files to be used as input.
#3: Output filename identifier (with path!)
#4: Number of cytosines
#5: Number of replicates (total number, cases+controls)


#Test that the number of arguments is correct
if (length(input_args)!=5){
  stop("The number of input arguments is incorrect.",call.=FALSE)
}


N_data_files=as.numeric(input_args[2])
N_cytosines=as.numeric(input_args[4])
N_replicates=as.numeric(input_args[5])/2

library(BiSeq)


#Go through the given input files and save the data into a single BSRaw object


#Objects for storing the window ranges
range_starts <- rep(0,N_data_files*2)
range_ends <- rep(0,N_data_files*2)
range_chromosomes <- rep(0,N_data_files*2) 
  
#Objects for storing the windows into BSRaw objects
window_starts <- rep(0,N_data_files*N_cytosines*2)
window_ends <- rep(0,N_data_files*N_cytosines*2)
window_chromosomes <- rep(0,N_data_files*N_cytosines*2)

totalReads_all <- matrix(NA,nrow=N_data_files*N_cytosines*2,ncol=N_replicates*2)
methReads_all <-matrix(NA,nrow=N_data_files*N_cytosines*2,ncol=N_replicates*2)

ind_i=1

for (d in 0:1){

  for (i in 1:N_data_files){

    print(paste("Loading data set",i,"diff",d,sep=" "))
  
    luxus_simulated <- read.table(paste(input_args[1],i-1,"_diff",d,".txt", sep = ""),skip=1,header = FALSE, sep = "",row.names=1)
    sample_names <- read.table(paste(input_args[1],i-1,"_diff",d,".txt", sep = ""),nrows=1,header = FALSE, sep = "",stringsAsFactors =FALSE)
  
    coordinate_info <- unlist(strsplit(rownames(luxus_simulated),":"))
    coordinates <- as.numeric(coordinate_info[c(FALSE,TRUE,FALSE,FALSE)])
  
    totalReads_window <- as.matrix(luxus_simulated[,c(TRUE,FALSE)])
    methReads_window <- as.matrix(luxus_simulated[,c(FALSE,TRUE)])
  

    #The replicates had to be sorted (cases and controls in certain order) for the M3D to work properly
    totalReads_sorted <- totalReads_window[,c(seq(1,ncol(totalReads_window),by=2),seq(2,ncol(totalReads_window),by=2))]
    methReads_sorted <- methReads_window[,c(seq(1,ncol(methReads_window),by=2),seq(2,ncol(methReads_window),by=2))]
  
    totalReads_all[((ind_i-1)*N_cytosines+1):(ind_i*N_cytosines),] <- totalReads_sorted
    methReads_all[((ind_i-1)*N_cytosines+1):(ind_i*N_cytosines),] <- methReads_sorted
  
    range_starts[ind_i] <- coordinates[1]
    range_ends[ind_i] <- coordinates[N_cytosines]
    range_chromosomes[ind_i] <- paste("chr",ind_i,sep="")
  
    window_starts[((ind_i-1)*N_cytosines+1):(ind_i*N_cytosines)] <-coordinates
    window_ends[((ind_i-1)*N_cytosines+1):(ind_i*N_cytosines)] <-coordinates
    window_chromosomes[((ind_i-1)*N_cytosines+1):(ind_i*N_cytosines)] <- rep(paste("chr",ind_i,sep=""),N_cytosines)

    ind_i=ind_i+1
  }
}

#Construct the BSRaw object

print("Constructing BSRaw object.")

metadata <- list(Sequencer = "Instrument", Year = "2019")

colData_all <- DataFrame(group = factor(c(rep(0,N_replicates),rep(1,N_replicates))),row.names = sample_names[c(seq(1,ncol(totalReads_window),by=2),seq(2,ncol(totalReads_window),by=2))])
colnames(totalReads_all)<-sample_names[c(seq(1,ncol(totalReads_window),by=2),seq(2,ncol(totalReads_window),by=2))]
colnames(methReads_all)<-sample_names[c(seq(1,ncol(methReads_window),by=2),seq(2,ncol(methReads_window),by=2))]
rowRanges_all <- GRanges(seqnames = window_chromosomes,ranges = IRanges(start = window_starts, end = window_ends))

luxus_simulated_BSraw<-BSraw(metadata = metadata,
                             rowRanges = rowRanges_all,
                             colData = colData_all,
                             totalReads = totalReads_all,
                             methReads = methReads_all)

#Run BiSeq
#Filtering is not needed for the simulated data set

data.clust.unlim <- clusterSites(object = luxus_simulated_BSraw,groups = colData(luxus_simulated_BSraw)$group,perc.samples = 1,min.sites = 10,max.dist = 1000)
quant <- quantile(totalReads(data.clust.unlim), 0.9)
predictedMeth <- predictMeth(object = data.clust.unlim)


betaResults <- betaRegression(formula=~group, link="probit",object=predictedMeth,type="BR")


#Null distribution

predictedMethNull<-predictedMeth
colData(predictedMethNull)$group.null<-rep(c(0,1),N_replicates)
betaResultsNull <- betaRegression(formula=~group.null, link="probit",object=predictedMethNull,type="BR")

#Do clustering and cluster trimming

vario<-makeVariogram(betaResultsNull)
vario.sm<-smoothVariogram(vario,sill=0.9)
vario.aux<-makeVariogram(betaResults,make.variogram=FALSE)
vario.sm$pValsList<-vario.aux$pValsList


locCor<-estLocCor(vario.sm)

clusters.rej<- testClusters(locCor,FDR.cluster=0.05)

clusters.trimmed<-trimClusters(clusters.rej,FDR.loc=0.2)


print("Saving p-values into files.")
write.table(betaResults, file=paste(input_args[3],"diff","both","_betaResults_default.txt",sep=""), row.names=TRUE, col.names=TRUE)
write.table(clusters.trimmed, file=paste(input_args[3],"diff","both","_trimmedClusters_default.txt",sep=""), row.names=TRUE, col.names=TRUE)
write.table(clusters.rej$clusters.not.reject, file=paste(input_args[3],"diff","both","_notRejectedClusters_default.txt",sep=""), row.names=TRUE, col.names=TRUE)
write.table(clusters.rej$clusters.reject, file=paste(input_args[3],"diff","both","_rejectedClusters_default.txt",sep=""), row.names=TRUE, col.names=TRUE)

