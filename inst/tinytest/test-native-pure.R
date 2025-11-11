source('pure.R')

data.frame(lapply(iris[,-5],cut,10))->X
X$const<-factor(rep(1,150))
X$tri<-factor(rep(1:3,50))
Y<-iris$Species
list(X=X,Y=Y,k=4)->input
c("CMI","MIM","JMIM","NJMIM","JMI","JMI3","DISR","CMIM","MRMR")->algos

for(algo in algos){
 do.call(sprintf("pure%s",algo),input)->pure
 do.call(algo,input)->native
 expect_equal(pure,native)
}

(function(){
 input$X$const<-NULL
 for(e in 1:5)
  input$X[[sprintf("const_%s",e)]]<-factor(rep(17,150))
 input$k<-10

 for(algo in algos){
  if(.Machine$sizeof.pointer!=8) skip("Unstable on i386")
  do.call(sprintf("pure%s",algo),input)->pure
  do.call(algo,input)->native
  expect_equal(pure,native)
 }
})()

(function(){
 input$X$spoiler<-factor(1:150)
 input$k<-3

 for(algo in algos){
  if(.Machine$sizeof.pointer!=8) skip("Unstable on i386")
  do.call(sprintf("pure%s",algo),input)->pure
  do.call(algo,input)->native
  expect_equal(pure,native)
 }
})()

MRMR(iris[,-5],iris[,5],positive=TRUE)->ans
expect_true(all(ans$score>0))

mutinfo<-function(x,y)
 .Call(praznik:::C_getMi,factor(x),factor(y))
 
expect_equal(
 apply(X,2,mutinfo,Y),
 miScores(X,Y)
)

condmutinfo<-function(x,y,z){
 unique(data.frame(x,y,z))->uxyz
 data.frame(t(apply(uxyz,1,function(xyz){
  c(
   pxyz=mean(x==xyz[1] & y==xyz[2] & z==xyz[3]),
   pxz=mean(x==xyz[1] & z==xyz[3]),
   pyz=mean(y==xyz[2] & z==xyz[3]),
   pz=mean(z==xyz[3])
  )
 })))->p
 sum(with(p,pxyz*log(pxyz*pz/pxz/pyz)))
}
Z<-factor((1:150)%%7)
expect_equal(
 apply(X,2,condmutinfo,Y,Z),
 cmiScores(X,Y,Z)
)

expect_equal(
 cmiScores(X,Y,factor(1:150)),
 apply(X,2,function(x) 0)
)
expect_equal(
 cmiScores(X,Y,factor(rep("Q",150))),
 miScores(X,Y)
)

entro<-function(x){
 table(x)/length(x)->p
 sum(-ifelse(p>0,p*log(p),0))
}
expect_equal(
 hScores(X),
 apply(X,2,entro)
)

expect_equal(
 hScores(data.frame(lapply(X,joinf,Y))),
  jhScores(X,Y)
 )

 Z<-factor((1:150)%%7)
 expect_equal(
  jmiScores(X,Y,factor(1:150)),
  setNames(rep(miScores(data.frame(Y),Y),ncol(X)),names(X))
 )
 expect_equal(
  jmiScores(X,Y,Z),
  cmiScores(X,Y,Z)+miScores(data.frame(Y),Z)
 )
 for(e in 1:ncol(X))
  expect_equal(
   miScores(X,Y)[e],
   jmiScores(X,Y,X[,e])[e]
  )

sapply(X,function(z) jmiScores(X,Y,z))->m
diag(m)<-0
expect_equal(
 apply(m,2,max), # by col also checks symmetry
 maxJmiScores(X,Y)
)

sapply(X,function(z) cmiScores(X,Y,z))->m
diag(m)<-NA
expect_equal(
 apply(m,1,max,na.rm=TRUE),
 maxCmiScores(X,Y)
)
expect_equal(
 apply(m,1,min,na.rm=TRUE),
 minCmiScores(X,Y)
)
expect_equal(
 apply(m,1,range,na.rm=TRUE),
 minMaxCmiScores(X,Y)
)

if(.Machine$sizeof.pointer==8)
 for(met in sapply(algos,get))
  expect_equal(
   met(iris[,rep(1:4,10)],iris$Species,threads=2),
   met(iris[,rep(1:4,10)],iris$Species,threads=1)
  )


pureImp<-function(X,Y){
 gi<-function(X,Y){
  k<-(k<-table(X,Y))/sum(k)
  sum(k^2/rowSums(k))-sum(colSums(k)^2)
 }
 apply(X,2,gi,Y)
}
expect_equal(impScores(X,Y),pureImp(X,Y))

data(MadelonD)
JIM(MadelonD$X,MadelonD$Y,20)->ans
expect_equal(pureImp(MadelonD$X[ans$selection[1]],MadelonD$Y),ans$score[1])
expect_true(all(grepl("^Rel",names(ans$selection))))
