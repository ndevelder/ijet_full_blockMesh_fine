library(ggplot2)

homedir = "/home/nate/"
#homedir = "C:/Users/ndeve/"

########### Set case info here ################################

casename = "ijet_full_blockMesh2"
casedir = paste(homedir,"foam/nate-3.2/run/",casename,"/",sep="")
savedir = paste(casedir,"analysis/",sep="")
ercodir = paste(savedir,"jetdata/Case_data/",sep="")

########### Note sample must have been run already ############
###############################################################

resultsdir = paste(savedir,"results",sep="")
dir.create(resultsdir, showWarnings = FALSE)
setwd(resultsdir)

dataIDs = c("00","10","15","20","25","30")
dlen = length(dataIDs)

ercodata <- list()
plotercodata <- data.frame()
xoffset = 0
dx = 3

for(did in dataIDs){

  #Format for meanU: y/D  U/UBulk
  ercofile = paste(ercodir,'ij2hr-',did,'-sw-mu.dat',sep="");
  ercodata[[did]] <- read.table(file=ercofile,sep = "",head = FALSE,skip = 4)

  ploteDataTemp <- data.frame(y=ercodata[[did]][,1],U=ercodata[[did]][,2],colname = rep(did,length(ercodata[[did]][,1])),xoff = rep(xoffset,length(ercodata[[did]][,1])))
  
  plotercodata <- rbind.data.frame(plotercodata,ploteDataTemp)
  
  xoffset = xoffset + dx
}


# Simulation Data Section

setsdir = paste(casedir,"sets",sep="")
testSample = dir.exists(setsdir) 

if(!testSample){
  cat("Execution halted...sets directory does not exist")
  stop
}

sampleTimes = list.dirs(setsdir, full.names = FALSE, recursive = FALSE)
lastSample = as.character(max(as.double(sampleTimes)))

sampleDir = paste(setsdir,"/",lastSample,sep="")
sampleFiles = list.files(sampleDir,full.names = TRUE)

simK <- list()
simY <- list()
simUx <- list()
simUy <- list()
simUz <- list()
simU <- list()

plotData <- data.frame()
xoffset = 0
uBulk = 10.185
pipeD = 0.1

for(did in dataIDs){
  
  scalarIndex = grep(paste("line",did,"_p",sep=""),sampleFiles)
  vectorIndex = grep(paste("line",did,"_U",sep=""),sampleFiles)

  scalarData = read.csv(file=sampleFiles[scalarIndex],sep = "",head = FALSE)
  vectorData = read.csv(file=sampleFiles[vectorIndex],sep = "",head = FALSE)

  sampleFileNames = list.files(sampleDir,full.names = FALSE)

  scalarNames = gsub(".xy","",sampleFileNames[scalarIndex])
  scalarNameArray <- strsplit(scalarNames,"_")
  vectorNames = gsub(".xy","",sampleFileNames[vectorIndex])
  vectorNameArray <- strsplit(vectorNames,"_")

  # Get simulation data
  simK[[did]] <- scalarData[,which(scalarNameArray[[1]] == "k")]
  simY[[did]] <- scalarData[,which(scalarNameArray[[1]] == paste("line",did,sep=""))]
  
  simUx[[did]] <- vectorData[,which(vectorNameArray[[1]] == "U")]
  simUy[[did]] <- vectorData[,which(vectorNameArray[[1]] == "U")+1]
  simUz[[did]] <- vectorData[,which(vectorNameArray[[1]] == "U")+2]
  simU[[did]] <- sqrt(simUx[[did]]^2 + simUy[[did]]^2 + simUz[[did]]^2)
  
  
  plotDataTemp <- data.frame(y=simY[[did]]/pipeD,U=simU[[did]]/uBulk,colname = rep(did,length(simY[[did]])),xoff = rep(xoffset,length(simY[[did]])))
  
  plotData <- rbind.data.frame(plotData,plotDataTemp)
  
  xoffset = xoffset + dx

}


  
# Construct Log-y Plot 
termplot <- ggplot()+
     geom_point(data=plotData[seq(1, nrow(plotData), 4),],aes(x=U,y=y,color=colname),shape=21,size=2)+
     geom_path(data=plotercodata,aes(x=U,y=y,color=colname))+
     xlab("U/(U bulk)")+
     ylab("y/D")+
     xlim(-1,2)+
    labs(color="r/D")+
       scale_y_continuous(minor_breaks=seq(0,0.34,by=0.05),limits = c(0,0.34))+
  facet_grid(~ colname)

   ggsave(filename="Umean_ijet.pdf", plot=termplot, width=9, height=9)

  
