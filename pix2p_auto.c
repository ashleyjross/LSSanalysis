/*
*
 */
#include "headera.h"
#include <stdlib.h>
//#include <malloc.h>
//#include <fstream>
//#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_math.h>
//#include <gsl/gsl_complex.h>
//#include <gsl/gsl_complex_math.h>
//#include <gsl/gsl_const_num.h>
//#include <gsl/gsl_const_mksa.h>
//#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_eigen.h>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_permutation.h>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_integration.h>
//#include <gsl/gsl_monte.h>
//#include <gsl/gsl_monte_plain.h>
//#include <gsl/gsl_monte_miser.h>
//#include <gsl/gsl_monte_vegas.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_odeiv.h>
const int   NR_END     = 1;

main(int argc, char *argv[]) {
	int nran = 100000;
	int Nline_max =1000;
	double a;
	double *lcd = dvector(1,1000000);
	double *lcr = dvector(1,1000000);
	double *lsd = dvector(1,1000000);
	double *lsr = dvector(1,1000000);
	double *lzw = dvector(1,1000000);
	double *lpc = dvector(1,1000000);
	int i = 0;
	int k = 0;
	int j = 0;
	int N = 0;
	double ti = time(NULL);
	int nbin = 20;
	double angm;
	angm = atof(argv[2]);
	//double minang = .05*pi/180.*angm;
	//double maxang = 10.*pi/180.*angm;
	double minang = .1*pi/180.*angm;
    double maxang = 10.*pi/180.*angm;
    double binl[nbin];
	double binlpix[nbin];
	double angl[nbin];
	double binedges[nbin+2];
	double g = minang/maxang;
	for(int nb=0;nb<nbin;nb++) {
		binl[nb] = 0;
		binlpix[nb] = 0;
	}
	for(int bn=0;bn<nbin+2;bn++) {
		double num = (nbin);
		double pn = bn/num;
		double pnu = (bn+1)/num;
		double angle = maxang * pow(g, pn);	
		double angb = (maxang * pow(g, pn)+maxang * pow(g, pnu))/2.;
		binedges[bn] = cos(angle);
		if (bn<nbin) {
			angl[bn] =angb;
		}
		//printf("%lf %lf %lf\n",angle,pn,cos(angle));
	}
	char fxi_name[200];
	char s1[200];
	sprintf(s1,argv[1]);
	strcpy(fxi_name,s1);
	sprintf(fxi_name,strcat(fxi_name,"odenspczw.dat"));
	FILE *fxi;
	printf("%s\n",s1);
	fxi=fopen(fxi_name,"r");
	if (fxi == NULL) perror ("Error opening file");
	int Nline=0;
	const int bsz=240;
	char buf[bsz];
	
	while(fgets(buf,bsz,fxi)!= NULL) {
		double cd,cr,sr,sd,zw,pc;
		//if(sscanf(buf,"%*d %f %f %f %f",&cd,&cr,&sd,&sd)!=4) 
		//	err_handler("xi input");
		//puts(buf);
		sscanf(buf,"%lf %lf %lf %lf %lf %lf",&sr,&cr,&sd,&cd,&zw,&pc);
		//cd = 1.;
		//cr = 1.;
		//sr = 1.;
		//sd = 1.;
		lcd[Nline] = cd;
		lcr[Nline] = cr;
		lsr[Nline] = sr;
		lsd[Nline] = sd;
		lzw[Nline] = zw;
		lpc[Nline] = pc;
		Nline++;
    }
	printf("%i \n",Nline);
/*	while (i<nran) {
		
		
		l[i] = rand();
		i = i+1;
		
	}
*/
	float be;
	float r;
	while (k<Nline)  {
		j = k+1;
		while (j<Nline)  {
			a = lcd[k]*lcd[j]*(lcr[k]*lcr[j] + lsr[k]*lsr[j]) + lsd[k]*lsd[j];
			//r = lsr[k];
			be = binedges[0];
			int ba = -1;
			while (a > be) {
				ba++;
				be = binedges[ba+1];
			}
			if (ba > -1 && ba < nbin){
				binlpix[ba] += lpc[k]*lpc[j];
				binl[ba] += lzw[k]*lzw[j]*lpc[k]*lpc[j];
			}
			j++;
			N+=1;
			
		}
		//printf("%f %f\n",be,r);
		k++;
	}
	double tf = time(NULL)-ti;
	printf("took %g\n",tf);
	char foname[200];
	//char sf1[200];
	strcat(s1,"2ptPixc.dat");
	

	sprintf(foname,s1);
	FILE *fo;
	fo =fopen(foname,"w");
	for (int ip=0; ip<nbin; ip++) {
		fprintf(fo,"%g %g \n",angl[ip]*180./pi,binl[ip]/binlpix[ip]);
	}
		//fclose(fxi);		
	//exit(0);
}

void err_handler(char *error_text)
{
    fprintf(stderr,"Run-time error...\n%s\n",error_text);
    fprintf(stderr,"now exiting to system...\n");
    exit(1);
}


double *dvector(long nl, long nh)
{
    void err_handler(char*);
    double *v;
    
    v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
    if (!v) err_handler("allocation failure in dvector()");
    return v-nl+NR_END;
}

