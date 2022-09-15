#library(ucminf)
library(pracma)
library(Matrix)
library(matrixStats)
library(tictoc)
library(pheatmap)

#read parameters and data for the linlog model
A<-as.matrix(read.delim("A.txt",sep = ",",header=FALSE))#,row.names = "none")
B<-as.matrix(read.delim("Bmatrix.txt",sep = ",",header=FALSE))#,row.names = "none")
C<-as.matrix(read.delim("C.txt",sep = ",",header=FALSE))#,row.names = "none")
N<-as.matrix(read.delim("Nmet_matrix.txt",sep = ",",header=FALSE))#,row.names = "none")
#reference abundances
X0<-as.matrix(read.delim("X0.txt",sep = ",",header=FALSE))
Y0<-as.matrix(read.delim("Y0.txt",sep = ",",header=FALSE))
e0<-as.matrix(read.delim("e0.txt",sep = ",",header=FALSE))
#metabolites
metabolites<-c("pep","g6p","pyr","f6p","f16p","dhap","pg3","ac-coA/coA","pg6","ru5p","r5p","s7p","kg2","suc","fum","mal","cit","gap","pgp","2pg","x5p","e4p","icit","succCoA/coA","gly","oaa")
c_metabolites<-c("GLC","ATP/ADP","ATP/AMP","NADPH/NADP","NADH/NAD","FADH2/FAD","ACE")
enzymes<-c("PtsG","Pgi","PfkA,,PfkB","FbaA,,FbaB","TpiA","GapA","Pgk","GpmA,,GpmB","Eno","PykA,,PykF","AceEF,LpdA,(PDH)","Zwf;,Pgl","Gnd","Rpe","RpiA,,RpiB","TktA,,TktB","TalA,,TalB","TktA,,TktB2","GltA,,PrpC","AcnA,,AcnB","IcdA","SucAB","SucCD","SdhABCD","FumA,,FumB,,FumC","Mdh","Ppc","MaeB,,MaeA","AceA","AceB","Edd;,Eda","Pta;,AckA,,AckB","PckA","Acs")

#define a certain number of windows to place genes
MX<-50
#and use it to derive random positions
random_pos<-MX*rand(34,1)
Y0<-as.matrix(Y0)

#define a range for beta. The fact we approximate the genome with
#50 windows makes so that the beta values are much larger than for normal genomes
ball=logspace(-1.5,-0.5,n=10);

#defining the objective function for the next optimization

objfun <- function(p,N,BB,A,C,y,ball) {
  #declaring a matrix for storing all metabolite steady state values
  Sss<-matrix(0,nrow=nrow(N),ncol=length(ball))
    #now we calculate ss levels when enzyme abundances change
  #as an effect of a change in the number of active replication forks
  #through beta
    for(i in 1:length(ball)){
        #new enzyme abundances based on gene position and the beta
        e<-1+ball[i]*p
        E<-as.matrix(Diagonal(n=34,x = e))
        #steady state solution
        tmpSss=-inv(N %*% E %*% BB) %*% (N %*% E %*% A + N %*% E %*% C %*% log(y))
        Sss[,i]<-tmpSss
      
    }
  
  
  objfinal <- 1E+3 * mean(rowSds(10^Sss) / (rowMeans(10^Sss)))
  
}


RES<-NULL
RES2<-NULL
tic()
for(i in 1:5){
  random_pos<-round(MX*rand(34,1),digits = 0)
  
optimNLMINB<-nlminb(start=random_pos,
                    objective = objfun,
                    A=A,B=B,C=C,N=N,y=Y0,ball=ball,
                    control = list(trace = 1,
                                   abs.tol = 1e-15,
                                   xf.tol=1e-12,
                                   iter.max=1000000,
                                   eval.max=1000000,
                                   x.tol=1e-12,
                                   step.min=1,
                                   step.max=MX-1),
                    lower=1,upper=MX)


RES<-cbind(RES,optimNLMINB$par)
RES2<-cbind(RES2,random_pos)
toc()
tic()
}

rownames(RES)<-enzymes
YlOrRd<-colorRampPalette(c("white",RColorBrewer::brewer.pal(12,"YlOrRd")))(33)

pheatmap(RES,color = YlOrRd)

