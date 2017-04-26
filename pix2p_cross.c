/*
 *  test.c
 *  
 *
 *  Created by ashleyr on 30/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include "headera.h"
#include <stdlib.h>
//#include <malloc.h>
//#include <fstream>
//#include <iomanip>
#include <stdio.h>
#include <math.h>
//#include <string.h>
//#include <string>
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
	double *lcd2 = dvector(1,1000000);
	double *lcr2 = dvector(1,1000000);
	double *lsd2 = dvector(1,1000000);
	double *lsr2 = dvector(1,1000000);
	double *lzw2 = dvector(1,1000000);
	double *lpc2 = dvector(1,1000000);
	int i = 0;
	int k = 0;
	int j = 0;
	int N = 0;
	double ti = time(NULL);
    double angm;
    angm = atof(argv[3]);
    double bsd;
    bsd = atof(argv[4]);
    double maxd;
    maxd = atof(argv[5]);
    printf("%lf %lf %lf\n",angm,bsd,maxd);
    //double minang = .2*pi/180.*angm;
    double minang = 0*pi/180.;
    double binsize = bsd*angm*pi/180.;
    double maxang = angm*maxd*pi/180.;
    //double binsize = 0.1*pi/180.;
    //double maxang = 8.*pi/180.;
    int nbin = (maxang-minang)/binsize;
    double binl[nbin];
    double angl[nbin];
    double binlpix[nbin];
    double binedges[nbin+2];
    //double g = minang/maxang;
    for(int nb=0;nb<nbin;nb++) {
        binl[nb] = 0;
        binlpix[nb] = 0;
    }
    for(int bn=0;bn<nbin+2;bn++) {
        double num = (bn);
        double angle = maxang - binsize*num;
        if (angle < 0) angle = 0;
        double angb = maxang-binsize*num-binsize/2.;
        binedges[bn] = cos(angle);
        if (bn<nbin) {
            angl[bn] =angb;
        }
        printf("%lf %lf %lf\n",angle,binedges[bn],cos(angle));
    }
	char fxi_name[200];
	char s1[200];
	char s2[200];
	sprintf(s1,argv[1]);
	sprintf(s2,argv[2]);
	//arg2 = argv[2];
	//printf("%s %s\n",s1,s2);
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
	printf("%s\n",s2);
	char f2_name[200];
	strcpy(f2_name,s2);
	sprintf(f2_name,strcat(f2_name,"odenspczw.dat"));
	FILE *f2;
	printf("%s \n",f2_name);
	f2=fopen(f2_name,"r");
	if (f2 == NULL) perror ("Error opening file");
	int Nline2=0;
	const int bsz2=240;
	char buf2[bsz2];
	
	while(fgets(buf2,bsz2,f2)!= NULL) {
		double cd2,cr2,sr2,sd2,zw2,pc2;
		//if(sscanf(buf,"%*d %f %f %f %f",&cd,&cr,&sd,&sd)!=4) 
		//	err_handler("xi input");
		//puts(buf);
		sscanf(buf2,"%lf %lf %lf %lf %lf %lf",&sr2,&cr2,&sd2,&cd2,&zw2,&pc2);
		//cd = 1.;
		//cr = 1.;
		//sr = 1.;
		//sd = 1.;
		lcd2[Nline2] = cd2;
		lcr2[Nline2] = cr2;
		lsr2[Nline2] = sr2;
		lsd2[Nline2] = sd2;
		lzw2[Nline2] = zw2;
		lpc2[Nline2] = pc2;
		Nline2++;
    }
	printf("%i %i\n",Nline,Nline2);
/*	while (i<nran) {
		
		
		l[i] = rand();
		i = i+1;
		
	}
*/
	float be;
	float r;
	while (k<Nline)  {
		j = 0;
		while (j<Nline2)  {
			a = lcd[k]*lcd2[j]*(lcr[k]*lcr2[j] + lsr[k]*lsr2[j]) + lsd[k]*lsd2[j];
			//r = lsr[k];
			be = binedges[0];
			int ba = -1;
			while (a > be) {
				ba++;
				be = binedges[ba+1];
			}
			if (ba > -1 && ba < nbin){
				binlpix[ba] += lpc[k]*lpc2[j];
				binl[ba] += lzw[k]*lzw2[j]*lpc[k]*lpc2[j];
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
	strcat(s2,"2ptPixc.dat");
	sprintf(foname,strcat(s1,s2));
	FILE *fo;
	fo =fopen(foname,"w");
	for (int ip=0; ip<nbin; ip++) {
		printf("%g %g",binl[ip],binlpix[ip]);
		fprintf(fo,"%g %g \n",angl[ip]*180./pi,binl[ip]/binlpix[ip]);
	}
	//fclose(fxi);		
	//exit(0);
}
