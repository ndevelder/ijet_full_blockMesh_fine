library(ggplot2)
library(plyr)
library(dplyr)
library(numDeriv)
library(reshape2)
library(zoo)
library(grid)

homedir = "/home/nate/"
#homedir = "C:/Users/ndeve/"
savedir = "Dropbox/Research/Data/Channel/"

datadir = paste(homedir,savedir,"Files/",sep="")
workingdir = paste(homedir,savedir,sep="")

########### Set case info here ################################
caseTime = 40
casename = "channel_3d_rough_phiOverK_new8_85_0.34759627_MGHPCC"
casedir = paste("/home/nate/foam/nate-3.2/run/",casename,"/",sep="")


########### Note sample must have been run already ############
###############################################################

resultsdir = paste(casedir,"results",sep="")
dir.create(resultsdir, showWarnings = FALSE)
setwd(resultsdir)

chanMeansLoc = paste(datadir,'chan590_clean.csv',sep="");
chanMeans <- read.csv(file=chanMeansLoc,sep = ",",head = TRUE)

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

scalarIndex = grep("epsilon",sampleFiles)
vectorIndex = grep("vorticity",sampleFiles)

scalarData = read.csv(file=sampleFiles[scalarIndex],sep = "",head = FALSE)
vectorData = read.csv(file=sampleFiles[vectorIndex],sep = "",head = FALSE)

sampleFileNames = list.files(sampleDir,full.names = FALSE)

scalarNames = gsub(".xy","",sampleFileNames[scalarIndex])
scalarNameArray <- strsplit(scalarNames,"_")
vectorNames = gsub(".xy","",sampleFileNames[vectorIndex])
vectorNameArray <- strsplit(vectorNames,"_")

# Get simulation data
simK = scalarData[,which(scalarNameArray[[1]] == "k")]
simY = scalarData[,which(scalarNameArray[[1]] == "lineY")]
simP = scalarData[,which(scalarNameArray[[1]] == "p")]
simEps = scalarData[,which(scalarNameArray[[1]] == "epsilon")]
simPhiOverK = scalarData[,which(scalarNameArray[[1]] == "tpphi")]
simPhi = scalarData[,which(scalarNameArray[[1]] == "tpphi")]*simK
simNut = scalarData[,which(scalarNameArray[[1]] == "nut")]
simProd = scalarData[,which(scalarNameArray[[1]] == "tpProd")]*simK
simEH = scalarData[,which(scalarNameArray[[1]] == "epsHat")]*simK

simUx = vectorData[,which(vectorNameArray[[1]] == "U")]
simPsi = vectorData[,which(vectorNameArray[[1]] == "tppsi")+4]*simK
simVort = vectorData[,which(vectorNameArray[[1]] == "vorticity")+6]
simPsiProd = scalarData[,which(scalarNameArray[[1]] == "psiProd")]
simSProd = scalarData[,which(scalarNameArray[[1]] == "sProd")]

simAlpha = 1/(1+1.5*(simPhi/simK))

simNutFrac = simNut/(simNut + 5*0.0017)

# Get channel dns data and constants
kappa = 0.41
R11_dns = chanMeans$R_uu
R22_dns = chanMeans$R_vv
R12_dns = chanMeans$R_uv
R33_dns = chanMeans$R_ww
k_dns = 0.5*(R11_dns + R22_dns + R33_dns)
k_dns[0] = 1e-15
tpphi_dns = R22_dns/k_dns
phi_dns = R22_dns
psi_dns = R12_dns
gradU_dns = chanMeans$dUmean.dy
y_dns = chanMeans$y
U_dns = chanMeans$Umean
nu_dns = 1/max(gradU_dns)
yPlus_dns = y_dns/nu_dns
tauW_dns = nu_dns*max(gradU_dns)
uTau_dns = sqrt(tauW_dns)
yPlus_dns = y_dns*uTau_dns/nu_dns
eps_dns = -1*chanMeans$dissip/nu_dns
prod_dns = abs(gradU_dns*psi_dns)
alpha_dns = 1/(1+1.5*(phi_dns/k_dns))
phiOverK_dns = phi_dns/k_dns
nut_dns = 0.21*phi_dns*k_dns/eps_dns

nu = nu_dns
ihalf = which(simY<=0)
cy = (simY[ihalf] + 1.0)
  
cases = cbind.data.frame(u=simUx[ihalf]/max(U_dns),k=simK[ihalf]/max(k_dns),eps=simEps[ihalf]/max(eps_dns),phi=simPhi[ihalf],psi=simPsi[ihalf],prod=simProd[ihalf]/max(prod_dns),y=cy,c=rep("simulation",times=length(cy)))
casem = melt(cases,id=c("y","c"))
  
dns = cbind.data.frame(u=U_dns/max(U_dns),k=k_dns/max(k_dns),eps=eps_dns/max(eps_dns),phi=phi_dns,psi=psi_dns,prod=prod_dns/max(prod_dns),y=y_dns,c=rep("kmm dns",times=length(y_dns)))
dnsm = melt(dns,id=c("y","c"))
  
combData = rbind(casem,dnsm)
  
# Construct Log-y Plot 
termplot <- ggplot()+
    geom_point(data=casem,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsm,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_y_continuous(minor_breaks=seq(-1.3,1.3,by=0.1))+
    scale_x_log10()+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_allvar_comparison_log.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)

  
  # Construct Plot 
  termplot <- ggplot()+
    geom_point(data=casem,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsm,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(-0.001,1.01)+
    scale_y_continuous(minor_breaks=seq(0.0,1.3,by=0.1))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_allvar_comparison.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)

  # Construct Plot 
  termplot <- ggplot()+
    geom_point(data=casem,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsm,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(-0.001,0.3)+
    scale_y_continuous(minor_breaks=seq(0.0,1.3,by=0.1))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_allvar_nearwall_comparison.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
# Velocity Only Plot
  
  casesU = cbind.data.frame(u=simUx[ihalf],y=cy,c=rep("simulation",times=length(cy)))
  caseU = melt(casesU,id=c("y","c"))
  
  dnssU = cbind.data.frame(u=U_dns,y=y_dns,c=rep("kmm dns",times=length(y_dns)))
  dnsU = melt(dnssU,id=c("y","c"))
  
  # Construct U Plot 
  termplot <- ggplot()+
    geom_point(data=caseU,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsU,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_y_continuous(minor_breaks=seq(0.0,22,by=0.5))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_U_comparison.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
  # Construct U Log-y Plot 
  termplot <- ggplot()+
    geom_point(data=caseU,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsU,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_x_log10()+
    scale_y_continuous(minor_breaks=seq(0.0,22,by=0.5))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_U_comparison_log.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
# Phi Only Plot
  casesPhi = cbind.data.frame(phi=simPhi[ihalf],y=cy,c=rep("simulation",times=length(cy)))
  casePhi = melt(casesPhi,id=c("y","c"))
  
  dnssPhi = cbind.data.frame(phi=phi_dns,y=y_dns,c=rep("kmm dns",times=length(y_dns)))
  dnsPhi = melt(dnssPhi,id=c("y","c"))
  
  # Construct U Plot 
  termplot <- ggplot()+
    geom_point(data=casePhi,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsPhi,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_y_continuous(minor_breaks=seq(0.0,1.5,by=0.1))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_phi_comparison.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
  # Construct U Log-y Plot 
  termplot <- ggplot()+
    geom_point(data=casePhi,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsPhi,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_x_log10()+
    scale_y_continuous(minor_breaks=seq(0.0,22,by=0.5))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_phi_comparison_log.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
  
  # PhiOverK Only Plot
  casesPhiOK = cbind.data.frame(phiOverK=simPhiOverK[ihalf],y=cy,c=rep("simulation",times=length(cy)))
  casePhiOK = melt(casesPhiOK,id=c("y","c"))
  
  dnssPhiOK = cbind.data.frame(phiOverK=phiOverK_dns,y=y_dns,c=rep("kmm dns",times=length(y_dns)))
  dnsPhiOK = melt(dnssPhiOK,id=c("y","c"))
  
  # Construct U Plot 
  termplot <- ggplot()+
    geom_point(data=casePhiOK,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsPhiOK,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_y_continuous(minor_breaks=seq(0.0,22,by=0.5))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_phiOverK_comparison.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
  
  
  # K Only Plot
  casesK = cbind.data.frame(k=simK[ihalf],y=cy,c=rep("simulation",times=length(cy)))
  caseK= melt(casesK,id=c("y","c"))
  
  dnssK = cbind.data.frame(k=k_dns,y=y_dns,c=rep("kmm dns",times=length(y_dns)))
  dnsK = melt(dnssK,id=c("y","c"))
  
# Construct k Plot 
  termplot <- ggplot()+
    geom_point(data=caseK,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsK,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_y_continuous(minor_breaks=seq(0.0,22,by=0.5))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_k_comparison.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
  # Construct U Log-y Plot 
  termplot <- ggplot()+
    geom_point(data=caseK,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsK,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_x_log10()+
    scale_y_continuous(minor_breaks=seq(0.0,22,by=0.5))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_k_comparison_log.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
  
  # Vorticity Only Plot
  casesV = cbind.data.frame(k=-1*simVort[ihalf],y=cy,c=rep("simulation",times=length(cy)))
  caseV= melt(casesV,id=c("y","c"))
  
  dnssV = cbind.data.frame(k=gradU_dns,y=y_dns,c=rep("kmm dns",times=length(y_dns)))
  dnsV = melt(dnssV,id=c("y","c"))
  
  # Construct k Plot 
  termplot <- ggplot()+
    geom_point(data=caseV,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsV,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_y_continuous(minor_breaks=seq(0.0,600,by=50))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_vorticity_comparison.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
  # Construct U Log-y Plot 
  termplot <- ggplot()+
    geom_point(data=caseV,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsV,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_x_log10()+
    scale_y_continuous(minor_breaks=seq(0.0,600,by=50))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_vorticity_comparison_log.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
  
  
  # EpsHat Only Plot
  casesEH = cbind.data.frame(v=simEH[ihalf],y=cy,c=rep("simulation",times=length(cy)))
  caseEH= melt(casesEH,id=c("y","c"))
  
  dnssEH = cbind.data.frame(v=eps_dns,y=y_dns,c=rep("kmm dns",times=length(y_dns)))
  dnsEH = melt(dnssEH,id=c("y","c"))
  
  # Construct EH Plot 
  termplot <- ggplot()+
    geom_point(data=caseEH,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsEH,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_y_continuous(minor_breaks=seq(0.0,600,by=50))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_epsHat_comparison.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
  # Construct EH Log-y Plot 
  termplot <- ggplot()+
    geom_point(data=caseEH,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsEH,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_x_log10()+
    scale_y_continuous(minor_breaks=seq(0.0,600,by=50))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_epsHat_comparison_log.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
  
  
  
  # Epsilon Only Plot
  casesE = cbind.data.frame(k=simEps[ihalf],y=cy,c=rep("simulation",times=length(cy)))
  caseE= melt(casesE,id=c("y","c"))
  
  dnssE = cbind.data.frame(k=eps_dns,y=y_dns,c=rep("kmm dns",times=length(y_dns)))
  dnsE = melt(dnssE,id=c("y","c"))
  
  # Construct eps Plot 
  termplot <- ggplot()+
    geom_point(data=caseE,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsE,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_y_continuous(minor_breaks=seq(0.0,600,by=50))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_epsilon_comparison.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
  # Construct eps Log-y Plot 
  termplot <- ggplot()+
    geom_point(data=caseE,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsE,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_x_log10()+
    scale_y_continuous(minor_breaks=seq(0.0,600,by=50))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_epsilon_comparison_log.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
  casesPE = cbind.data.frame(k=simProd[ihalf]/simEps[ihalf],y=cy,c=rep("prod/epsilon",times=length(cy)))
  casePE= melt(casesPE,id=c("y","c"))
  
  # Construct P/eps Plot 
  termplot <- ggplot()+
    geom_point(data=casePE,aes(x=y,y=value,color=c))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_y_continuous(minor_breaks=seq(0.0,2,by=0.1))+
    coord_fixed(ratio=0.5)
  filenamerp = paste(casename,"_prodepsilonratio_comparison.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  

  
  # Alpha Plot
  casesA = cbind.data.frame(alpha=simAlpha[ihalf],y=cy,c=rep("simulation",times=length(cy)))
  caseA= melt(casesA,id=c("y","c"))
  
  dnssA = cbind.data.frame(alpha=alpha_dns,y=y_dns,c=rep("kmm dns",times=length(y_dns)))
  dnsA = melt(dnssA,id=c("y","c"))
  
  # Construct alpha Plot 
  termplot <- ggplot()+
    geom_point(data=caseA,aes(x=y,y=value,color=variable))+
    geom_path(data=dnsA,aes(x=y,y=value,color=variable))+
    xlab("y")+
    xlim(0.001,1.01)+
    scale_y_continuous(minor_breaks=seq(0.0,600,by=50))+
    theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+
    coord_fixed(ratio=0.13)
  filenamerp = paste(casename,"_alpha_comparison.pdf",sep="")
  ggsave(filename=filenamerp, plot=termplot, width=6, height=6) 

  
 # Prod Compare
 casesPr = cbind.data.frame(prod=simProd[ihalf],y=cy,c=rep("prod simulation",times=length(cy)))
 casePr= melt(casesPr,id=c("y","c"))
 
#casesPrs = cbind.data.frame(k=simSProd[ihalf],y=cy,c=rep("Sprod simulation",times=length(cy)))
 #casePrs = melt(casesPrs,id=c("y","c"))
 
 dnssPr = cbind.data.frame(prod=-1*psi_dns*gradU_dns,y=y_dns,c=rep("kmm dns",times=length(y_dns)))
 dnsPr = melt(dnssPr,id=c("y","c"))
 
 # Construct Prod Compare Plot 
 termplot <- ggplot()+
   geom_point(data=casePr,aes(x=y,y=value,color=variable))+
   #geom_point(data=casePrs,aes(x=y,y=value,color=variable))+
   geom_path(data=dnsPr,aes(x=y,y=value,color=variable))+
   xlab("y")+
   xlim(0.001,1.01)+
   scale_x_log10()+
   scale_y_continuous(minor_breaks=seq(0.0,22,by=0.5))+
   theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+    coord_fixed(ratio=0.13)
 filenamerp = paste(casename,"_prod_comparison_log.pdf",sep="")
 ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  
 # Construct Prod Compare Plot 
 termplot <- ggplot()+
   geom_point(data=casePr,aes(x=y,y=value,color=variable))+
   #geom_point(data=casePrs,aes(x=y,y=value,color=variable))+
   geom_path(data=dnsPr,aes(x=y,y=value,color=variable))+
   xlab("y")+
   xlim(0.001,1.01)+
   scale_y_continuous(minor_breaks=seq(0.0,22,by=0.5))+
   theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+    coord_fixed(ratio=0.13)
 filenamerp = paste(casename,"_prod_comparison.pdf",sep="")
 ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
 
 
 
 
 # Nut Compare
 casesN = cbind.data.frame(nut=simNut[ihalf],y=cy,c=rep("nut simulation",times=length(cy)))
 caseN= melt(casesN,id=c("y","c"))
 
 dnssN = cbind.data.frame(nut=nut_dns,y=y_dns,c=rep("kmm dns",times=length(y_dns)))
 dnsN = melt(dnssN,id=c("y","c"))
 
 # Construct Nut Compare Plot 
 termplot <- ggplot()+
   geom_point(data=caseN,aes(x=y,y=value,color=variable))+
   geom_path(data=dnsN,aes(x=y,y=value,color=variable))+
   xlab("y")+
   xlim(0.001,1.01)+
   scale_x_log10()+
   scale_y_log10(breaks=c(0.001,0.005,0.01,0.05,0.1,0.5,1),limits=c(0.0001,1))+
   theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+    coord_fixed(ratio=0.13)
 filenamerp = paste(casename,"_nut_comparison_log.pdf",sep="")
 ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
 
 # Construct Nut Compare Plot 
 termplot <- ggplot()+
   geom_point(data=caseN,aes(x=y,y=value,color=variable))+
   geom_path(data=dnsN,aes(x=y,y=value,color=variable))+
   xlab("y")+
   xlim(0.001,1.01)+
   scale_y_continuous(minor_breaks=seq(0.0,22,by=0.5))+
   theme(legend.justification=c(1,0),legend.position=c(1,0),plot.margin=unit(c(0,0,0,0),"cm"),panel.margin=unit(c(0,0,0,0),"cm"),aspect.ratio=1)+    coord_fixed(ratio=0.13)
 filenamerp = paste(casename,"_nut_comparison.pdf",sep="")
 ggsave(filename=filenamerp, plot=termplot, width=6, height=6)
  