#include <stdlib.h>
//#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

const double pi    = 3.1415926536; // Pi

// house keeping functions: util.c
void   err_handler(char*);
double *dvector(long,long);
void   free_dvector(double*,long,long);
float  *vector(long,long);
void   free_vector(float*,long,long);
int    *ivector(long,long);
void   free_ivector(int*,long,long);
char   *cvector(long,long);
void   free_cvector(char*,long,long);
float  **matrix(long,long,long,long);
void   free_matrix(float**,long,long,long,long);
int    **imatrix(long,long,long,long); 
void   free_imatrix(int**,long,long,long,long);
double **dmatrix(long,long,long,long);
void   free_dmatrix(double**,long,long,long,long);
int    **imatrix(long,long,long,long);
void   free_imatrix(int**,long,long,long,long);
double  ***d3tensor(long,long,long,long,long,long);
void   free_d3tensor(double***,long,long,long,long,long,long);
float  ***f3tensor(long,long,long,long,long,long);
void   free_f3tensor(float***,long,long,long,long,long,long);
float  ****f4tensor(long,long,long,long,long,long,long,long);
void   free_f4tensor(float****,long,long,long,long,long,long,long,long);
int    **vimatrix(long,long,long,int*);
void   free_vimatrix(int**,long,long,long);
float  ran2(long*);
float  gasdev(long*);
float poidev(float,long*);

// integration functions: integration.c
double qsimp(double (*func)(double),double,double);
double qsimp2(double (*func)(double),double,double);
double qmidinf(double (*func)(double,double),double,double,double);
double qmidinf2(double (*func)(double,double,double),double,double,double,double);
double qsimpmid(double (*func)(double,double,double),double,double,double,double);
double qsimpmid2(double (*func)(double,double),double,double,double);

// PS mass functions: ps_mass_func.c
void   nm_power();
double Pk(double); 
double Sigth(double,double);
double dSigthbydR(double,double);

// chisq optimisation routines: fit_param.c
double golden(double,double,double,double(*f)(double),double,double*);
void powell(double[],double**,int,double,int*,double*,double (*func)(double []));

//  matrix inversion: mat_inv.c
void mat_inv(double**,double**,int,double*);
void dsvdcmp(double**,int,int,double[],double**);

// linear algebra programs: lin_alg.c
void   evect_qrql(double**,int,double[]);
void   jacobi(double**,int,int,double[],double**,int*);
void   eigsrt(double[],double**,int,int);
int   choldc(double**,int,double[]);
void   cholsl(double**,int,double[],double[],double[],double*);
void   choldc_invertL(double**,int,int);

// spline fitting routines: spline.c
void   spline(double[],double[],int,double,double,double[]);
void   splint(double[],double[],double[],int,double,double*);
void   splind(double[],double[],double[],int,double,double*);

// SDSS table input: table_input.c
void read_gal_file(struct use_gal*,char*,long*,FILE*,double,double,int);
