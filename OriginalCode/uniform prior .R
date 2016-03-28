############################################
## ShuaiQi's Project
## Date 	03-28-2016
## Aim: 	Try Baysesian on whold dataset UNIFORM
## @ authors: 	SQ
## Data source: /dscrhome/gd44/SQProject/RStan/2016/exon_level_process_v2.txt
## Data source: current directory: exon_level_process_v2.txt
## Models: 	Bayesian Stan
## Parameters:	See Stan model part
## Outputs: 	The same repository, *.Rout document
## 
## 
## Data source and copyright?
## Read in table from local hard drive:
## setup the working directory to where exon_level_process_v2.txt locates
   
##########################################
## Part One Read in the data
## create the tables
##########################################
## set current working directory
## setwd("/Users/shuaiqizhang/Desktop/project /data ")

## read-in source data from the same directory
table<-read.table("exon_level_process_v2.txt")
colnames(table)<- c("chr", "gene", "dom", "subdom", "exon", "gene.dom", 
             "gene.dom.subdom", 
             "envarp",    # pass
             "envarpf",   # pass functional
             "envarpfr",  # pass functional rare
             "emutr")     # mutation rate 

## data manipulation
table<-within(table,envarpfc<-envarpf-envarpfr)#y
table<-within(table,gene<-factor(gene))
table<-within(table,gene.dom<-factor(gene.dom))
table<-within(table,gene.dom.subdom<-factor(gene.dom.subdom))
table<-table[which(table$envarp!=0), ] #exclude exons with 0 variant site from data 
table$x=scale(table$envarp)

## for the use of counting number of gene
sumenvarp<-aggregate(table$envarp, by=list(Category=table$gene), FUN=sum)
sumenvarpfc<-aggregate(table$envarpfc, by=list(Category=table$gene), FUN=sum)[,2]
table1<-data.frame(cbind(sumenvarp,sumenvarpfc))
colnames(table1)<-c("gene","sumenvarp","sumenvarpfc")




##############################################
## RStan Model
##############################################




hiernormalwide<-"
data{ #get the data set 
	int<lower=0> N;   # number of exon level
	int<lower=0> J;   # number of gene level
	int <lower=1,upper=J> gene[N];  #index of gene
	int <lower=1,upper=N> exon[N];  #index of exon 
	vector[N] xij;   #x at exon level
	vector[N] yij; #y at exon level
}
parameters{ #specify the parameter we want to know 
	real beta;  #common slope for the exon level
	real mu;      #common intercept for the exon level
	vector[N] aij; #random intercept for the exon level
	real <lower=0> sigma_aj[J];  #variance of intercept at exon level 
	vector[J] aj; #random intercept for the gene level 
	real <lower=0> sigma_a;  #variance of intercept at gene level
	real <lower=0> sigma; #variance of yij
}
transformed parameters{ #specify the model we will use 
}
model { #give the prior distribution
        vector[N] lambdaij; #exon level model
	for (i in 1:N)
     		lambdaij[i] <- mu+beta*xij[i]+aij[exon[i]]+aj[gene[i]];#specify the gene group
   	beta ~normal(0,1);
   	sigma ~ uniform(0, 20);
   	sigma_a ~uniform(0,10);
   	sigma_aj ~uniform(0,10);
   	aj ~ normal(0, sigma_a);
   	yij ~ normal(lambdaij,sigma);
   	for (i in 1:N)
   	    aij[i]~ normal(0,sigma_aj[gene[i]]);
}
"

##############################################
## load rstan package

library("rstan")


##############################################
## pass values to parameters
J <- dim(table1)[1] #gene number 
N <- dim(table)[1]  #exon number 
xij = c(table$envarp)
xj = table1$sumenvarp
yij = table$envarpfc
yj = table1$sumenvarpfc
gene <- as.numeric(table$gene) #list
genelevel <- length(unique(gene)) #number
indexg <- match(gene, unique(gene))  #list
exon <- c(1:length(table$envarpfc))



##############################################
## M1_table
M1_table<-list(N=N, 
               J=J, 
               xij=xij,
               yij=yij,
               gene=indexg, 
               exon=exon)
control=list(adapt_delta=0.99,max_treedepth=12)

###############################################
## fit data to model
fitwide<-stan(model_code=hiernormalwide, data=M1_table,iter=100000,warmup=90000,chains=4)

print(fitwide)   ##print fitwide

answeruniform<-extract(fitwide,permuted=TRUE)

print(answeruniform) #print samples 