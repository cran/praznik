void k2ab(int e,int *a,int *b){
 int k=e+1;
 double jf=(sqrt(8.*((double)k)-7.)-1.)/2.;
 int j=floor(jf)+2;
 int i=k-(j-1)*(j-2)/2;
 a[0]=i-1;
 b[0]=j-1;
}

void cmp_two(SEXP A,SEXP B,int *out){
 int n=length(A);
 if(isInteger(A) && isInteger(B)){
  int *a=INTEGER(A),
      *b=INTEGER(B);
  for(int e=0;e<n;e++)
   if((a[e]==NA_INTEGER)||(b[e]==NA_INTEGER)){
    out[e]=NA_LOGICAL;
   }else{
    out[e]=a[e]>=b[e];
   }
 }else if(isReal(A) && isReal(B)){
  double *a=REAL(A),
         *b=REAL(B);
  for(int e=0;e<n;e++)
   if((a[e]==NA_REAL)||(b[e]==NA_REAL)||ISNAN(a[e])||ISNAN(b[e])){
    out[e]=NA_LOGICAL;
   }else{
    out[e]=a[e]>=b[e];
   }
 }else if((isReal(A) && isInteger(B))||(isInteger(A) && isReal(B))){
  double *a;
  int *b;
  Rboolean swap;
  if(isReal(A)){
   a=REAL(A);
   b=INTEGER(B);
   swap=FALSE;
  }else{
   a=REAL(B);
   b=INTEGER(A);
   swap=TRUE;
  }
  for(int e=0;e<n;e++){
   if((a[e]==NA_REAL)||ISNAN(a[e])||(b[e]==NA_INTEGER)){
    out[e]=NA_LOGICAL;
   }else{
    out[e]=(!swap)?(a[e]>=((double)b[e])):(b[e]>=((double)a[e]));
   }
  }
 }else{
  error("Incomparable types!");
 }
}

double unif_rand_pos(){
 double ans;
 do{
  ans=unif_rand(); //unif_rand is in [0..1)
 }while(ans==0.);
 return(ans);
}

SEXP C_tsp(SEXP X,SEXP Sep,SEXP Sample){
 if(!isDataFrame(X)) error("Input must be a data.frame");
 if(length(X)<2) error("Data frame with less than two columns");
 int m=length(X);
 int n=length(VECTOR_ELT(X,0));
 int nice_names=!isNull(Sep);


 const char *sep=nice_names?CHAR(STRING_ELT(Sep,0)):NULL;
 int buflen;
 SEXP Names=PROTECT(getAttrib(X,R_NamesSymbol));
 if(nice_names){
  int ncx=0;
  for(int e=0;e<m;e++){
   int c=strlen(CHAR(STRING_ELT(Names,e)));
   if(c>ncx) ncx=c;
  }
  buflen=strlen(sep)+ncx*2+1;
 }else{
  buflen=floor(log10(m*(m-1)/2+1))+3; //Plus one for digit, one for \0 and one for T 
 }
 char *buf=(char*)R_alloc(sizeof(char),buflen);

 int k=asInteger(Sample);
 if(k==NA_INTEGER || k<0 || k>(m*(m-1)/2)) error("Invalid sample argument, must be between one and m*(m-1)/2");
 int s=(k>0)?k:m*(m-1)/2; //s -- how many actually
 int *resv=(k>0)?(int*)R_alloc(sizeof(int),s):NULL;
 for(int e=0;e<k;e++) resv[e]=e;

 if(k){
  GetRNGstate();
  for(int e=0;e<k;e++) resv[e]=e;
  int e=k-1;
  int max=m*(m-1)/2;
  double w=exp(log(unif_rand_pos())/((double)k));

  while(1){
   e+=floor(log(unif_rand_pos())/log(1.-w))+1;
   if(e>=max) break;
   resv[(int)R_unif_index(k)]=e;
   w=w*exp(log(unif_rand_pos())/k);
   if(w>=1.) break; //Almost impossible
  }
  PutRNGstate();
 }
 
 SEXP Ans=PROTECT(allocVector(VECSXP,s));
 SEXP AnsN=PROTECT(allocVector(STRSXP,s));
 setAttrib(Ans,R_NamesSymbol,AnsN);

 for(int e=0;e<s;e++){
  int a,b,ec=resv?resv[e]:e;
  k2ab(ec,&a,&b);
  if(nice_names){
   snprintf(buf,buflen,"%s%s%s",CHAR(STRING_ELT(Names,a)),sep,CHAR(STRING_ELT(Names,b)));
  }else{
   snprintf(buf,buflen,"T%d",ec+1);
  }
  SET_STRING_ELT(AnsN,e,mkChar(buf));
  SEXP Col=PROTECT(allocVector(LGLSXP,n));
  SET_VECTOR_ELT(Ans,e,Col);
  UNPROTECT(1);
  int *out=LOGICAL(Col);
  cmp_two(VECTOR_ELT(X,a),VECTOR_ELT(X,b),out);
 }
 UNPROTECT(3);
 return(Ans);
}

