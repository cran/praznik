library(praznik)

data(MadelonD)
set.seed(17)

input<-list(X=MadelonD$X,Y=MadelonD$Y,k=20,threads=2)
pmr<-sample(nrow(MadelonD$X),nrow(MadelonD$X))
pmc<-sample(ncol(MadelonD$X),ncol(MadelonD$X))
inputp<-list(X=MadelonD$X[pmr,pmc],Y=MadelonD$Y[pmr],k=20,threads=2)

for(algo in c("MIM","JMIM","NJMIM","JMI","DISR","CMIM","MRMR","JIM")){
 do.call(algo,input)->j
 do.call(algo,inputp)->p
 stopifnot(all.equal(j$score,p$score))
 stopifnot(identical(names(j$selection),names(p$selection)))
}

# We need exception for CMI since attributes become indistinguishable for it
inputp<-list(X=MadelonD$X[pmr,],Y=MadelonD$Y[pmr],k=20,threads=2)
do.call(CMI,input)->j
do.call(CMI,inputp)->p
stopifnot(all.equal(j$score,p$score))
stopifnot(identical(names(j$selection),names(p$selection)))
