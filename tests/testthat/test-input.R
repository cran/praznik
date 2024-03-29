test_that("Crazy decision works",{
 data.frame(A=factor(rep(1:2,each=20)))->X
 rep(c(TRUE,FALSE),each=20)->Y
 expect_equal(names(MIM(X,Y,1)$selection),"A")
 expect_equal(MIM(X,rep(-.7,40),1),list(selection=integer(0),score=numeric(0)))
})

test_that("Crazy attributes work",{
 badfactor<-factor(c(rep(c("z","a"),each=7),letters))[1:14]
 data.frame(
  bool=rep(c(TRUE,FALSE),each=7),
  int=as.integer(rep(c(-37,33),each=7)),
  badfactor=badfactor
 )->X
 factor(rep(letters[1:2],each=7))->Y
 MIM(X,Y,3)->S
 expect_true(all(names(S$selection)%in%names(X)))
 expect_equal(setNames(S$score,NULL),setNames(rep(S$score[1],3),NULL))
})

test_that("X must be a data.frame is thrown",{
 expect_error(MIM(list(1:3),NULL,NULL),"X must be a data.frame")
})

test_that("No rows error is thrown",{
 expect_error(MIM(iris[0,-5],iris[0,5]),"X has no rows")
 expect_error(hScores(numeric(0)),"X has a zero length")
})

test_that("Nameless data.frames work",{
 z<-iris[,-5]
 names(z)<-NULL
 JMI(z,iris$Species,4)->ans
 expect_null(names(ans$score))
 expect_null(names(ans$selection))
 expect_true(all(ans$selection%in%(1:4)))
})

test_that("Constant real features work",{
 #Throws segfault in 1.0.0
 MIM(data.frame(a=rep(1,150)),iris$Species,1)->ans
 expect_equal(length(ans$selection),0)
})

test_that("Crazy-double valued features work",{
 Y<-factor(c("a","b"))
 X<-data.frame(
  a=c(0,.1),
  b=c(.Machine$double.eps,.Machine$double.neg.eps),
  c=c(.Machine$double.neg.eps,0),
  d=c(.Machine$double.eps,0),
  e=c(.Machine$double.xmax,.Machine$double.xmin),
  f=c(.Machine$double.eps,.Machine$double.xmin),
  g=c(.Machine$double.eps,.Machine$double.xmax),
  h=c(.Machine$double.neg.eps,.Machine$double.xmin),
  i=c(.Machine$double.neg.eps,.Machine$double.xmax),
  j=c(0,.Machine$double.xmin),
  k=c(0,.Machine$double.xmax),
  l=c(-1,.Machine$double.xmax)
 )
 expect_equal(
  miScores(X,Y),
  setNames(rep(log(2),ncol(X)),names(X))
 )
})

#Following https://gitlab.com/mbq/praznik/issues/23
c(
 "CMI","MIM","JMIM","NJMIM","JMI","DISR","CMIM","MRMR",
 "cmiMatrix","cmiScores","dnmiMatrix","hScores","impScores",
 "jmiMatrix","jmiScores","miMatrix","miScores","njmiMatrix",
 "njmiScores","nmiMatrix","maxJmiScores","maxCmiScores",
 "minCmiScores"
)->algos
for(algo in algos){
 test_that(sprintf("Named vector X work with %s",algo),{
  if(grepl("^h",algo)){
   z<-do.call(algo,list(c(a=1,b=2,c=3)))
  }else if(grepl("(^mi|^imp|^dnmi|^nmi|^maxJmi|^maxCmi|^minCmi)",algo)){
   z<-do.call(algo,list(c(a=1,b=2,c=3),1:3))
  }else if(grepl("(^cmi|^jmi|^njmi)",algo)){
   z<-do.call(algo,list(c(a=1,b=2,c=3),1:3,3:1))
  }else{  
   z<-do.call(algo,list(c(a=1,b=2,c=3),1:3,k=1))
  }
  if(grepl("Matrix$",algo)){
   expect_null(rownames(z))
   expect_null(colnames(z))
  }else if(grepl("Scores$",algo)){
   expect_null(names(z))
  }else{
   expect_null(names(z$selection))
   expect_null(names(z$score))
  }
 })
}

test_that("Magical cut works like R cut",{
 expect_equal(
  miScores(iris[,-5],iris$Species),
  miScores(data.frame(lapply(iris[,-5],cut,10)),iris$Species)
 )
 ii<-iris[c(1:14,100:112),]
 expect_equal(
  miScores(ii[,-5],ii$Species),
  miScores(data.frame(lapply(ii[,-5],cut,9)),ii$Species)
 )
 ii<-iris[c(1:3,101:103),]
 expect_equal(
  miScores(ii[,-5],ii$Species),
  miScores(data.frame(lapply(ii[,-5],cut,2)),ii$Species)
 )
})

test_that("Zero-score features work",{
 expand.grid(a=c(T,F),b=c(T,F),n1=c(T,F),n2=c(T,F),n3=c(T,F))->X
 X$aub<-X$a|X$b
 X$anb<-X$a&X$b
 Y<-X$a!=X$b
 for(e in c(JMI,DISR,MRMR))
  expect_equal(sort(names(e(X,Y,ncol(X))$selection)),sort(names(X)))
 for(e in c(MIM,CMIM,JMIM,NJMIM)){
  e(X,Y,ncol(X))->ans
  expect_true(all(ans$score>0))
  expect_true(all(names(ans$selection)%in%names(X)))
  expect_false(any(duplicated(ans$selection)))
 }
})

test_that("MIM works with zero-score features",{
 #Tests for the following error reported by smilesun:
 # https://github.com/mbq/praznik/issues/9
 x<-data.frame(a=rep(1,150),b=rep(2,150),c=rep(3,150))
 y<-iris$Species
 expect_identical(MIM(x,y,3),DISR(x,y,3))
})

test_that("X and Y must be only reals, booleans, integers or factors",{
 Y<-c(TRUE,TRUE,FALSE,FALSE,FALSE)
 li<-data.frame(A=1:5)
 li$A<-list(1,1:2,1:3,1:4,1:5)
 badX<-list(
  char=data.frame(A=letters[1:5],stringsAsFactors=FALSE),
  realna=data.frame(A=c((1:4)*5.5,NA)),
  realinf=data.frame(A=c(1:4,Inf)),
  img=data.frame(A=1:5+3i),
  li=li
 )
 for(X in badX){
  expect_error(MIM(X,Y,1))
  expect_error(hScores(X$A))
 }
 for(Y in badX)
  expect_error(MIM(data.frame(A=Y),X$A,1))
})

test_that("NAs and other quirks are caught",{
 Y<-iris$Species; Y[3]<-NA
 X<-iris[,-5]; X[12,3]<-NA
 expect_error(MIM(X,iris$Species,1),"Non-finite numeric values are not allowed")
 X[12,3]<-Inf
 expect_error(MIM(X,iris$Species,1),"Non-finite numeric values are not allowed")
 X[12,3]<-NaN
 expect_error(MIM(X,iris$Species,1),"Non-finite numeric values are not allowed")
 X<-iris[,"Species",drop=FALSE]; X[17,1]<-NA
 expect_error(MIM(X,iris$Species,1),"NA values are not allowed")
 X<-iris[,"Species",drop=FALSE]
 X[,1]<-as.integer(X[,1]); X[17,1]<-NA
 expect_error(MIM(X,iris$Species,1),"NA values are not allowed")
})

test_that("scorers throws bad Z size",{
 expect_error(cmiScores(iris[,-5],iris[,5],1:10),"Z vector size mismatch")
 expect_error(jmiScores(iris[,-5],iris[,5],1:10),"Z vector size mismatch")
 expect_error(njmiScores(iris[,-5],iris[,5],1:10),"Z vector size mismatch")
})

test_that("threads argument is processed well",{
 expect_error(miScores(iris[,-5],iris[,5],-17L),"Invalid threads argument")
 expect_error(miScores(iris[,-5],iris[,5],NA),"Invalid threads argument")
 expect_warning(miScores(iris[,-5],iris[,5],1+parallel::detectCores()),"Thread count capped")
})

test_that("triScores throws on too narrow input",{
 expect_error(triScores(iris[,1,drop=FALSE]),"Cannot process")
 expect_error(triScores(1:10),"Cannot process")
 expect_error(triScores(iris[,1:2]),"Cannot process")
 expect_error(triScores(iris[1:2,rep(1,2346)]),"Too many features")
})

test_that("JMI3 works with k=1",{
 expect_equal(
  pureJMI3(iris[,-5],iris[,5],k=1),
  JMI3(iris[,-5],iris[,5],k=1)
 )
})

test_that("JMI3 works with k=2",{
 expect_equal(
  pureJMI3(iris[,-5],iris[,5],k=2),
  JMI3(iris[,-5],iris[,5],k=2)
 )
})
