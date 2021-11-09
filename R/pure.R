# Pure R implementations of selection methods 

#' @importFrom stats setNames
#' @importFrom utils tail


pureMIM<-function(X,Y,k=3){
 miScores(X,Y)->mim
 sort(mim,decreasing=TRUE)[1:k]->ans
 ans[ans>0]->ans
 list(
  selection=setNames(match(names(ans),names(X)),names(ans)),
  score=ans
 )
}

pureCMIM<-function(X,Y,k=3){
 nX<-names(X)
 X<-data.frame(X)
 ascores<-miScores(X,Y)
 selection<-names(which.max(ascores))
 fscores<-max(ascores)
 #Conceptually wrong, but as defined by Fleuret
 scores<-ascores
 if(k>1) for(e in 1:(k-1)){
  joinf(X[,tail(selection,1)])->w
  scores[colnames(X)!=tail(selection,1)]->scores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X
  newScores<-cmiScores(X,Y,w)
  scores<-pmin(
   newScores,  
   scores
  )
  if(max(scores)==0) break

  selection<-c(selection,names(which.max(scores)))
  fscores<-c(fscores,max(scores))
 }
 list(
  selection=setNames(match(selection,nX),selection),
  score=setNames(fscores,selection)
 )
}

pureJMIM<-function(X,Y,k=3){
 nX<-names(X)
 X<-data.frame(X)
 ascores<-miScores(X,Y)
 selection<-names(which.max(ascores))
 fscores<-max(ascores)
 scores<-rep(Inf,ncol(X))
 if(k>1) for(e in 1:(k-1)){
  joinf(X[,tail(selection,1)])->x
  scores[colnames(X)!=tail(selection,1)]->scores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X
  newScores<-jmiScores(X,Y,x) 
  scores<-pmin(
   newScores,
   scores
  )
  if(max(scores)==0) break

  selection<-c(selection,names(which.max(scores)))
  fscores<-c(fscores,max(scores))
 }
 list(
  selection=setNames(match(selection,nX),selection),
  score=setNames(fscores,selection)
 )
}

pureNJMIM<-function(X,Y,k=3){
 nX<-names(X)
 X<-data.frame(X)
 ascores<-miScores(X,Y)
 selection<-names(which.max(ascores))
 fscores<-max(ascores)
 scores<-rep(Inf,ncol(X))
 if(k>1) for(e in 1:(k-1)){
  joinf(X[,tail(selection,1)])->x
  scores[colnames(X)!=tail(selection,1)]->scores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X
  newScores<-njmiScores(X,Y,x) 
  scores<-pmin(
   newScores,
   scores
  )
  if(max(scores)==0) break

  selection<-c(selection,names(which.max(scores)))
  fscores<-c(fscores,max(scores))
 }
 list(
  selection=setNames(match(selection,nX),selection),
  score=setNames(fscores,selection)
 )
}

pureCMI<-function(X,Y,k=3){
 nX<-names(X)
 X<-data.frame(X)
 ascores<-miScores(X,Y)
 selection<-names(which.max(ascores))
 Z<-X[,which.max(ascores)]
 X<-X[,-which.max(ascores),drop=FALSE]
 fscores<-max(ascores)
 if(k>1) for(e in 1:(k-1)){
  ascores<-cmiScores(X,Y,Z)
  if(max(ascores)==0) break

  selection<-c(selection,names(which.max(ascores)))
  fscores<-c(fscores,max(ascores))

  Z<-joinf(Z,X[,which.max(ascores)])
  X<-X[,-which.max(ascores),drop=FALSE]
 }
 list(
  selection=setNames(match(selection,nX),selection),
  score=setNames(fscores,selection)
 )
}

pureJMI<-function(X,Y,k=3){
 nX<-names(X)
 X<-data.frame(X)
 ascores<-miScores(X,Y)
 selection<-names(which.max(ascores))
 fscores<-max(ascores)
 scores<-rep(0,ncol(X))
 if(k>1) for(e in 1:(k-1)){
  joinf(X[,tail(selection,1)])->x
  scores[colnames(X)!=tail(selection,1)]->scores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X
  scores<-scores+jmiScores(X,Y,x) 

  selection<-c(selection,names(which.max(scores)))
  fscores<-c(fscores,max(scores))
 }
 list(
  selection=setNames(match(selection,nX),selection),
  score=setNames(fscores,selection)
 )
}

pureJMI3<-function(X,Y,k=3){
 nX<-names(X)
 X<-data.frame(X)
 W<-c()
 if(k>0){
  s1<-miScores(X,Y)
  fscores<-max(s1)
  U<-X[,which.max(s1),drop=FALSE]
  X<-X[,-which.max(s1),drop=FALSE]
 }
 if(k>1){
  s2<-jmiScores(X,Y,U[,1])
  fscores<-c(fscores,max(s2))
  W<-X[,which.max(s2),drop=FALSE]
  X<-X[,-which.max(s2),drop=FALSE]
 }
 if(k>2){
  scores<-rep(0,ncol(X))
  for(e in 3:k){
   for(u in U){
    scores<-scores+jmiScores(X,Y,joinf(u,W[,1]))
   }
   U<-cbind(U,W)
   i<-which.max(scores)
   fscores<-c(fscores,scores[i])
   scores<-scores[-i]
   W<-X[,i,drop=FALSE]
   X<-X[,-i,drop=FALSE]
  }
 }
 selection<-c(names(U),names(W))
 list(
  selection=setNames(match(selection,nX),selection),
  score=setNames(fscores,selection)
 )
}

pureMRMR<-function(X,Y,k=3){
 nX<-names(X)
 X<-data.frame(X)
 rel<-miScores(X,Y)
 red<-rep(0,ncol(X))
 selection<-names(which.max(rel))
 fscores<-max(rel)
 if(k>1) for(e in 1:(k-1)){
  joinf(X[,tail(selection,1)])->x
  rel[colnames(X)!=tail(selection,1)]->rel
  red[colnames(X)!=tail(selection,1)]->red
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X

  nred<-miScores(X,x)
  red<-red+nred;
  scores<-rel-red/e

  selection<-c(selection,names(which.max(scores)))
  fscores<-c(fscores,max(scores))
 }
 list(
  selection=setNames(match(selection,nX),selection),
  score=setNames(fscores,selection)
 )
}

pureDISR<-function(X,Y,k=3){
 nX<-names(X)
 X<-data.frame(X)
 ascores<-miScores(X,Y)
 selection<-names(which.max(ascores))
 fscores<-max(ascores)
 rep(0,ncol(X))->scores
 if(k>1) for(e in 1:(k-1)){
  joinf(X[,tail(selection,1)])->x
  scores[colnames(X)!=tail(selection,1)]->scores
  X[,colnames(X)!=tail(selection,1),drop=FALSE]->X
  
  scores<-scores+njmiScores(X,Y,x)

  selection<-c(selection,names(which.max(scores)))
  fscores<-c(fscores,max(scores))
 }
 list(
  selection=setNames(match(selection,nX),selection),
  score=setNames(fscores,selection)
 )
}

