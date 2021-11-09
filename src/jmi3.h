SEXP C_JMI3(SEXP X,SEXP Y,SEXP K,SEXP Threads){
 int n,k,m,ny,*y,*nx,**x,nt;
 struct ht **hta;
 prepareInput(X,Y,K,Threads,&hta,&n,&m,&k,&y,&ny,&x,&nx,&nt);

 double bs=0.; int *cY,*ctmp,bi=0;
 initialMiScan(hta,n,m,y,ny,x,nx,&cY,&ctmp,NULL,&bs,&bi,nt);
 if(bs==0) return(makeAns(0,NULL,NULL));

 //Hold pointers to selected features
 int **u=(int**)R_alloc(sizeof(int*),(k<=2)?1:(k-2));

 //Save selected X as W and discard from further consideration
 int *w=x[bi],nw=nx[bi]; x[bi]=NULL;

 //Yet put it as a first selected attribute
 double *score; int *idx;
 SEXP Ans; PROTECT(Ans=makeAns(k,&score,&idx));
 score[0]=bs; idx[0]=bi+1; u[0]=w;
 
 //Time for an actual algorithm
 double *as=(double*)R_alloc(sizeof(double),m); //Accumulated score
 for(int e=0;e<m;e++) as[e]=0.;
 int *wxc=(int*)R_alloc(sizeof(int),n*nt),*cWXc=ctmp;
 int *uwxc=(int*)R_alloc(sizeof(int),n*nt);
 bs=0.;

 #pragma omp parallel num_threads(nt)
 for(int e=1;e<k;e++){
  double tbs=0.;
  int tbi=-1,tn=omp_get_thread_num();
  struct ht *ht=hta[tn];
  int *wx=wxc+(tn*n),*cWX=cWXc+(tn*n),*uwx=uwxc+(tn*n);
  int *cUWX=cWX; //cUWX is used exclusively with cWX

  #pragma omp for
  for(int ee=0;ee<m;ee++){
   //Ignore attributes already selected
   if(!x[ee]) continue;

   //Mix x[ee] with x making wx
   int nwx=fillHt(ht,n,nx[ee],x[ee],nw,w,wx,NULL,NULL,1);
   if(e==1){
    //Make MI of mix and Y and increase its accumulated score
    fillHt(ht,n,ny,y,nwx,wx,NULL,NULL,cWX,0);
    as[ee]=miHt(ht,cY,cWX); 
   }else{
    for(int eu=0;eu<e-1;eu++){
     int nu=nx[idx[eu]-1];
     //Mix wx with u to make uwx
     int nuwx=fillHt(ht,n,nwx,wx,nu,u[eu],uwx,NULL,NULL,1);
     fillHt(ht,n,ny,y,nuwx,uwx,NULL,NULL,cUWX,0); 
     as[ee]=as[ee]+miHt(ht,cY,cUWX);
    }
   }

   if(as[ee]>tbs){
    tbs=as[ee]; tbi=ee;
   }
  }
  #pragma omp critical
  if((tbs>bs) || (tbs==bs && tbi<bi)){
   bs=tbs;
   bi=tbi;
  }
  
  #pragma omp barrier 
  #pragma omp single
  {
   //Reset accumulated scores after the second selection
   if(e==1) for(int ee=0;ee<m;ee++) as[ee]=0.;
   //Do not expand u in the last iteration
   if(e<k) u[e-1]=w;
   w=x[bi]; nw=nx[bi]; x[bi]=NULL; 
   score[e]=bs; idx[e]=bi+1;
   bs=0.;
  }
 }

 Ans=finishAns(k,Ans,X);

 UNPROTECT(1);
 return(Ans);
}

