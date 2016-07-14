#include <stdlib.h>
//#include <malloc.h>
#include <stdio.h>
#include <math.h>

const int   NR_END     = 1;  

// std error handler
void err_handler(char *error_text)
{
  fprintf(stderr,"Run-time error...\n%s\n",error_text);
  fprintf(stderr,"now exiting to system...\n");
  exit(1);
}


// allocate an unsigned long vector with subscript range v[nl..nh]
unsigned long *lvector(long nl, long nh)
{
  void err_handler(char*);
  unsigned long *v;
  
  v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
  if (!v) err_handler("allocation failure in lvector()");
  return v-nl+NR_END;
}

// allocate a double vector with subscript range v[nl..nh]
double *dvector(long nl, long nh)
{
  void err_handler(char*);
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) err_handler("allocation failure in dvector()");
  return v-nl+NR_END;
}

// free a double dvector allocated with dvector()
void free_dvector(double *v, long nl, long nh)
{
  free((char*) (v+nl-NR_END));
}

// allocate a float vector with subscript range v[nl..nh]
float *vector(long nl, long nh)
{
  void err_handler(char*);
  float *v;

  v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
  if (!v) err_handler("allocation failure in vector()");
  return v-nl+NR_END;
}

// free an float vector allocated with vector()
void free_vector(float *v, long nl, long nh)
{
  free((char*) (v+nl-NR_END));
}

// allocate an char vector with subscript range v[nl..nh]
char *cvector(long nl, long nh)
{
  void err_handler(char*);
  char *v;
  v=(char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(char)));
  if (!v) err_handler("allocation failure in cvector()");
  return v-nl+NR_END;
}

// free an unsigned char vector allocated with cvector() 
void free_cvector(char *v, long nl, long nh)
{
  free((char*) (v+nl-NR_END));
}

// allocate an int vector with subscript range v[nl..nh]
int *ivector(long nl, long nh)
{
  void err_handler(char*);
  int *v;
  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) err_handler("allocation failure in ivector()");
  return v-nl+NR_END;
}

// free an int vector allocated with ivector()
void free_ivector(int *v, long nl, long nh)
{
  free((char *) (v+nl-NR_END));
}

// allocate a float matrix with subscript range m[nrl..nrh][ncl..nch]
float **matrix(long nrl, long nrh, long ncl, long nch)
{
  void err_handler(char*);
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m) err_handler("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl]) err_handler("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

// free a float matrix allocated by matrix()
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
{
  free((char*) (m[nrl]+ncl-NR_END));
  free((char*) (m+nrl-NR_END));
}

// allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]
double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
  void err_handler(char*);
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;
  
  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) err_handler("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) err_handler("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

// free a double matrix allocated by dmatrix()
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
  free((char*) (m[nrl]+ncl-NR_END));
  free((char*) (m+nrl-NR_END));
}

// allocate a char matrix with subscript range m[nrl..nrh][ncl..nch]
char **cmatrix(long nrl, long nrh, long ncl, long nch)
{
  void err_handler(char*);
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  char **m;
  
  // allocate pointers to rows
  m=(char **) malloc((size_t)((nrow+NR_END)*sizeof(char*)));
  if (!m) err_handler("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;
  
  
  // allocate rows and set pointers to them
  m[nrl]=(char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char)));
  if (!m[nrl]) err_handler("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  // return pointer to array of pointers to rows
  return m;
}

// free an int matrix allocated by imatrix()
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
{
  free((char*) (m[nrl]+ncl-NR_END));
  free((char*) (m+nrl-NR_END));
}

// allocate a int matrix with subscript range m[nrl..nrh][ncl..nch]
int **imatrix(long nrl, long nrh, long ncl, long nch)
{
  void err_handler(char*);
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  // allocate pointers to rows
  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) err_handler("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;


  // allocate rows and set pointers to them
  m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) err_handler("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  // return pointer to array of pointers to rows
  return m;
}

// free an int matrix allocated by cmatrix()
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
  free((int*) (m[nrl]+ncl-NR_END));
  free((int*) (m+nrl-NR_END));
}

// allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh]
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  void err_handler(char*);
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;
  
  // allocate pointers to pointers to rows
  t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
  if (!t) err_handler("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;
  
  // allocate pointers to rows and set pointers to them
  t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
  if (!t[nrl]) err_handler("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  // allocate rows and set pointers to them
  t[nrl][ncl]=(float *) 
    malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
  if (!t[nrl][ncl]) err_handler("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  
  // return pointer to array of pointers to rows
  return t;
}

// free a float f3tensor allocated by f3tensor()
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
{
  free((char*) (t[nrl][ncl]+ndl-NR_END));
  free((char*) (t[nrl]+ncl-NR_END));
  free((char*) (t+nrl-NR_END));
}

// allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh]
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  void err_handler(char*);
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;
  
  // allocate pointers to pointers to rows
  t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
  if (!t) err_handler("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;
  
  // allocate pointers to rows and set pointers to them
  t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
  if (!t[nrl]) err_handler("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  // allocate rows and set pointers to them
  t[nrl][ncl]=(double *) 
    malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
  if (!t[nrl][ncl]) err_handler("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  
  // return pointer to array of pointers to rows
  return t;
}

// free a double f3tensor allocated by f3tensor()
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
{
  free((char*) (t[nrl][ncl]+ndl-NR_END));
  free((char*) (t[nrl]+ncl-NR_END));
  free((char*) (t+nrl-NR_END));
}

// random number generator

const int   IM1=2147483563;
const int   IM2=2147483399;
const float AM=(1.0/2147483563.);
const int   IMM1=(2147483563-1);
const int   IA1=40014;
const int   IA2=40692;
const int   IQ1=53668;
const int   IQ2=52774;
const int   IR1=12211;
const int   IR2=3791;
const int   RAN2_NTAB=32;
const int   NDIV=(1+(2147483563-1)/32);
const float RAN2_EPS=1.2e-7;
const float RAN2_RNMX=(1.0-1.2e-7);

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[32];
  float temp;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1; else *idum = -(*idum);
    idum2=(*idum);
    for (j=RAN2_NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < RAN2_NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RAN2_RNMX) return RAN2_RNMX; else return temp;
}

float gasdev(long *idum)
{
  float ran1(long *idum);
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*ran2(idum)-1.0;
      v2=2.0*ran2(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void indexx(unsigned long n, double arr[], unsigned long indx[])
{
  void err_handler(char*);
  unsigned long i,indxt,ir=n,itemp,j,k,l=1;
  int jstack=0,*istack;
  double a;
  
  istack=ivector(1,NSTACK);
  for (j=1;j<=n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	indxt=indx[j];
	a=arr[indxt];
	for (i=j-1;i>=1;i--) {
	  if (arr[indx[i]] <= a) break;
	  indx[i+1]=indx[i];
	}
	indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1]);
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	SWAP(indx[l+1],indx[ir])
	  }
      if (arr[indx[l]] > arr[indx[ir]]) {
	SWAP(indx[l],indx[ir])
	  }
      if (arr[indx[l+1]] > arr[indx[l]]) {
	SWAP(indx[l+1],indx[l])
	  }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for (;;) {
	do i++; while (arr[indx[i]] < a);
	do j--; while (arr[indx[j]] > a);
	if (j < i) break;
	SWAP(indx[i],indx[j])
	  }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack > NSTACK) err_handler("NSTACK too small in indexx.");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP


float gammln(float xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

#define PI 3.141592654

//float poidev(float xm, long *idum)
//{
//  float gammln(float xx);
//  float ran2(long *idum);
//  static float sq,alxm,g,oldm=(-1.0);
//  float em,t,y;
//  
//  if (xm < 12.0) {
//    if (xm != oldm) {
//      oldm=xm;
//      g=exp(-xm);
//    }
//    em = -1;
//    t=1.0;
//    do {
//      ++em;
//      t *= ran2(idum);
//    } while (t > g);
//  } else {
//    if (xm != oldm) {
//      oldm=xm;
//      sq=sqrt(2.0*xm);
//      alxm=log(xm);
//      g=xm*alxm-gammln(xm+1.0);
//    }
//    do {
//      do {
//	y=tan(PI*ran2(idum));
//	em=sq*y+xm;
//      } while (em < 0.0);
//      em=floor(em);
//      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
//    } while (ran2(idum) > t);
//  }
//  return em;
//}
#undef PI
