library(nnls)
library(Rcpp,lib.loc="/home/slcxding/R_libs")
library(RcppArmadillo,lib.loc="/home/slcxding/R_libs")
library(truncnorm,lib.loc="/home/slcxding/R_libs")
sourceCpp("/home/slcxding/File/MAGE_GMM.cpp")

################################################################################
###########################        haplotype data generation      #######################
################################################################################

set.seed(600)

ns=10000  ## Set the people number 10000
haplotype = read.table("/home/slcxding/File/haplotype.txt",header=F) #choose the 0-1 haplotype file
n.1 = dim(haplotype)[1]
haplotype = as.matrix(haplotype)
sample1=sample(1:n.1,ns,replace=TRUE)
sample2=sample(1:n.1,ns,replace=TRUE)
genotype=haplotype[sample1,]+haplotype[sample2,] #generate the SNP data
nonraredup.common=genotype[,which(colSums(genotype)>2*ns*0.05)] # rare allel frequency >0.05
nonraredup.common=as.matrix(nonraredup.common)
n.2 = dim(nonraredup.common)[2]  ## The dimension is 10000*576, 10000 people, 576SNP

print(n.2)
nonraredup.rare=genotype[,which(colSums(genotype)>2*ns*0.005 & colSums(genotype)<2*ns*0.05)] 
nonraredup.rare=as.matrix(nonraredup.rare)
n.3 = dim(nonraredup.rare)[2]  
print(n.3)

####### Initialization #########

n = 1000 # iteration numbers
ss = 2000 # sample size
sn_common = 10 # Common SNP numbers
sn_rare = 40 # Rare SNP numbers
cn_common = 4 # Common causal SNP numbers
cn_rare = 6 # Rare causal SNP numbers
n_main = 10
n_common.main = 2
n_rare.main = 8
alpha_common_main = 0.09
alpha_rare_main = 0.17
max_permutation = 10000000
initial_permutation = 1000
threshold = 100

scaler = dbeta(0.05,1,25)/dbeta(0.05,0.5,0.5)
weight_fun = function(maf) ifelse(maf>=0.05, scaler*dbeta(maf,0.5,0.5), dbeta(maf,1,25))

###### Claculation Function#########

alpha_common = 0.19
alpha_rare = 0.59

#######################################
   
for (rep in 1:n){
  x1 = rnorm(ss, mean=62.4, sd=11.5)											
  x2 = rbinom(ss, prob=0.52, size=1)											
  e1 = rbinom(ss, prob=0.5, size=1)	  
  indicator = sample(1:10000, ss, replace=F)
      
  indi_common = sample(1:n.2,sn_common,replace = F)
  none_indi_common = sample(1:sn_common,cn_common)
  G_common = rep(10,sn_common*ss)
  G_common = matrix(G_common,ncol=sn_common)
  SNP=NULL
  for (j in 1:sn_common)
  {
    SNP = nonraredup.common[,indi_common[j]]										
    G_common[,j] = SNP[indicator];
  }
      
  indi_rare = sample(1:n.3,sn_rare,replace = F)
  none_indi_rare = sample(1:sn_rare,cn_rare)
  G_rare = rep(10,sn_rare*ss)
  G_rare = matrix(G_rare,ncol=sn_rare)
  SNP=NULL
  for (j in 1:sn_rare)
  {
    SNP = nonraredup.rare[,indi_rare[j]]										
    G_rare[,j] = SNP[indicator];
  }

  indi.common.main = sample(1:sn_common, n_common.main, replace=F)
  G.common.main = G_common[, indi.common.main]
      
  indi.rare.main = sample(1:sn_rare, n_rare.main, replace=F)
  G.rare.main = G_rare[, indi.rare.main]

  G = cbind(G_common, G_rare)
      
  colnames(G)=rep("gene",dim(G)[2])
  for(i in 1:dim(G)[2]){
    colnames(G)[i] = paste("gene",i,sep="")
  }

  epsilon = rnorm(ss, mean=0, sd=1.5)
  SNP_common = G_common[,none_indi_common]
  SNP_rare = G_rare[,none_indi_rare]
  SNP = cbind(SNP_common, SNP_rare)
      
  alpha1 = runif(n_common.main, min=alpha_common_main-0.02, max=alpha_common_main+0.02)
  alpha2 = runif(n_rare.main, min=alpha_rare_main-0.02, max=alpha_rare_main+0.02)
  beta1 = runif(cn_common, min=alpha_common-0.02, max=alpha_common+0.02)
  beta2 = runif(cn_rare, min=alpha_rare-0.02, max=alpha_rare+0.02)
      
  SNP.sign = rep(c(1,-1),(cn_common+cn_rare))
  y_common = y_rare = y_common.main = y_rare.main = rep(0,ss)
  
  for (i in 1:cn_common){
    y_common = y_common + beta1[i]*SNP[,i]*SNP.sign[i]*e1
  }
  
  for (j in (cn_common+1):(cn_common+cn_rare)){
    y_rare = y_rare + beta2[j-cn_common]*SNP[,j]*SNP.sign[j]*e1
  }
  
  if (n_common.main == 1){
    y_common.main = y_common.main + alpha1*G.common.main*SNP.sign[1]
  } else {
    for (k in 1:n_common.main){
      y_common.main = y_common.main + alpha1[k]*G.common.main[,k]*SNP.sign[k]
    }
  }
  
  if(n_rare.main == 1){
    y_rare.main = y_rare.main + alpha2*G.rare.main*SNP.sign[10]
  } else{
    for (t in 1:n_rare.main){
      y_rare.main = y_rare.main + alpha2[t]*G.rare.main[,t]*SNP.sign[t+n_common.main]
    }
  }
      
  y = 0.05*x1+0.057*x2+0.64*e1+epsilon
  y = y-mean(y)
      
  X = cbind(x1,x2)
  sub.n = dim(G)[1]
  snp.n = dim(G)[2]  
  maf = colSums(G)/(sub.n*2)
  weight.num = weight_fun(maf)
      
  if (length(weight.num)==1){
    weight = weight.num
    M = ComputeProjxx(X, e1) ## cpp
    G.weight = G*weight
  }else{
    weight = diag(weight.num)
    M = ComputeProjxx(X, e1) ## cpp
    G.weight = MatMult(G,weight) ## cpp
  }
      
  S.weight = MatMult(diag(e1),G.weight) ## cpp
  eigen = sym_eigenvectors(M) ## cpp
  n.zero = dim(X)[2]+2
  A = eigen[,-c(1:n.zero)]
  KS = GetLinearKernel(S.weight)
  KG = GetLinearKernel(G.weight)
  T1 = c(kernel_product(A,KG))
  T2 = c(kernel_product(A,KS))
  T3 = c(diag(sub.n-n.zero))
  TT = cbind(T1,T2,T3)
  V = c(calculateV(A,y))
  N = length(V)
      
  fit = nnls(cbind(T1,T2,T3),V)
  original_estimate = coef(fit)[2]
  print(original_estimate)
  count_extreme = 0
  n_permutation = initial_permutation
  total_permutation_done = 0
  
  while (n_permutation <= max_permutation){
    for (i in 1:(n_permutation-total_permutation_done)){
      G_permuted = G[sample(1:nrow(G), nrow(G)),]
      maf = colSums(G_permuted)/(sub.n*2)
      weight.num = weight_fun(maf)
      
      if (length(weight.num)==1){
        weight = weight.num
        G.weight = G_permuted*weight
      }else{
        weight = diag(weight.num)
        G.weight = MatMult(G_permuted,weight) ## cpp
      }
      
      S.weight = MatMult(diag(e1),G.weight) ## cpp
      
      KS = GetLinearKernel(S.weight)
      KG = GetLinearKernel(G.weight)
      
      T1 = c(kernel_product(A,KG))
      T2 = c(kernel_product(A,KS))
      TT = cbind(T1,T2,T3)
      
      fit_permuted = nnls(cbind(T1,T2,T3),V)
      estimated_permuted = coef(fit_permuted)[2]
      if (estimated_permuted >= original_estimate){
        count_extreme = count_extreme + 1
      }
    }
    total_permutation_done = n_permutation
    if (count_extreme >= threshold){
      break
    } else{
      n_permutation = min(n_permutation*10, max_permutation)
    }
  }      
  pvalue = count_extreme/total_permutation_done
  cat(c(pvalue), "\n", file="pvalue.txt", append=TRUE)
}
    
    
q("no")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    