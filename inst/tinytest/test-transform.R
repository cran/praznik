
for(n in 5:10){
 x<-runif(n)
 tx<-kTransform(x)
 expect_equal(tx,kTransform(x*2))
 expect_equal(length(tx),n*(n-1))

 cx<-cut(x,2)
 kTransform(cx)->tcx
 expect_true(is.factor(tcx))
 expect_equal(length(tcx),n*(n-1))
 expect_error(
  kTransform(cut(x,3)),
  "Unordered factor with more"
 )
 expect_equal(levels(tcx),c("<",">","="))

 kTransform(1:n)->kn
 kTransform(n:1)->kr
 expect_identical(sort(kn),sort(kr))
 expect_equal(length(table(factor(kn))),2)
 kTransform(rep(pi,n))->kz
 expect_true(all(kz==kz[1]))
}

expect_error(
 kTransform("abc"),
 "Unsupported"
)
expect_error(
 kTransform(data.frame(a=letters,stringsAsFactors=FALSE)),
 "Unsupported"
)
expect_error(
 kTransform(NULL),
 "Unsupported"
)
eo4<-factor(c("<","<","<",">","=","<",">","=","<",">",">",">"),levels=c("<",">","="))
eo3<-factor(c("<","<",">","<",">",">"),levels=c("<",">","="))
eo2<-factor(c("<",">"),levels=c("<",">","="))
eona<-factor(c(NA,"<",NA,NA,">",NA),levels=c("<",">","="))
expect_equal(kTransform(c(FALSE,TRUE)),eo2)
expect_equal(kTransform(c(-Inf,0,Inf)),eo3)
expect_equal(kTransform(c(1.1,2.2,2.2,3.3)),eo4)
expect_equal(kTransform(c(1.1,NA,3.3)),eona)
expect_equal(kTransform(c(1.1,NaN,3.3)),eona)
expect_equal(kTransform(c(FALSE,NA,TRUE)),eona)
expect_equal(kTransform(c(0,NA,1)),eona)
expect_equal(kTransform(ordered(c("a","b","b","c"))),eo4)

expect_equal(
 sort(as.character(kTransform(c(1,1,2,0)))),
 sort(c('=','>','<','=','>','<','<','<','<','>','>','>'))
)

x<-iris[,-5]
expect_equal(
 kTransform(x),
 data.frame(lapply(x,kTransform))
)

for(e in 1:20){
 set.seed(e)
 x<-sample(runif(1,1,20),runif(1,3,20),replace=TRUE)
 kx<-kTransform(x)
 expect_equal(kInverse(kx),rank(x))
}

expect_equal(
 kInverse((kTransform(1:5)==">")*pi),
 1:5
)

expect_error(kInverse(iris$Species),"Factor does not seem")
expect_error(kInverse(1:5),"Invalid size")
expect_error(kInverse(iris),"Invalid input")
expect_error(kInverse(c(1i,2i)),"Invalid input")
expect_error(kInverse(integer(0)),"Input too short")
expect_error(kInverse(1),"Input too short")
 
Z<-data.frame(a=1:3,b=3:1,c=rep(2.,3))
Zo<-data.frame(a__b=c(F,T,T),a__c=c(F,T,T),b__c=c(T,T,F))
expect_identical(tspTransform(Z),Zo)
expect_identical(names(tspTransform(Z,sep='~')),c('a~b','a~c','b~c'))
expect_identical(names(tspTransform(Z,sep=NULL)),c('T1','T2','T3'))
expect_identical(sum(names(tspTransform(Z,sample=2))%in%names(Zo)),2L)
expect_error(tspTransform(Z,sample=20),'Invalid sample')
expect_error(tspTransform(iris),'Incomparable')
Z$a<-as.numeric(Z$a)
expect_identical(tspTransform(Z),Zo)

Z<-data.frame(a=c(1,NA,3),b=c(3,NA,1),c=c(NA,3,1))
Zo<-data.frame(a__b=c(F,NA,T),a__c=c(NA,NA,T),b__c=c(NA,NA,T))
expect_identical(tspTransform(Z),Zo)
Z$b<-as.integer(Z$b)
expect_identical(tspTransform(Z),Zo)
Z$c<-as.integer(Z$c)
expect_identical(tspTransform(Z),Zo)

