/*
 *  note, doesn't allow files with more than one million lines; split up larger files
 *  inputs are file1 file2 
 * (both files with pre-computed angle info; default angmin, zdmin = 0, 200 linearly spaced bins in theta, dz with max corresponding to 500 Mpc/h, 30 zbins between 0.4 and .7)
 *  
 *	Created by ashleyr on 30/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
/* compile with
 gcc -Wall -c utila.c
 gcc -Wall -c pp2pt_Dmufb.c
 gcc -Wall -o pp2pt_Dmufb pp2pt_Dmufb.o utila.o

 */

#include "headera.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_linalg.h>
const char *ddir={"/mnt/lustre/ashleyr/pcadw/"};
const char *odir={"/mnt/lustre/ashleyr/paircounts/"};

main(int argc, char *argv[]) {
	double *lcd = dvector(1,2000000);
	double *lcr = dvector(1,2000000);
	double *lsd = dvector(1,2000000);
	double *lsr = dvector(1,2000000);
	double *lw = dvector(1,2000000);
	double *lzd = dvector(1,2000000);
	double *lcd2 = dvector(1,2000000);
	double *lcr2 = dvector(1,2000000);
	double *lsd2 = dvector(1,2000000);
	double *lsr2 = dvector(1,2000000);
	double *lw2 = dvector(1,2000000);
	double *lzd2 = dvector(1,2000000);
	int i = 0;
	int k = 0;
	int j = 0;
	int N = 0;
	double ptot = 0;
	double ti = time(NULL);
	int nrbin = 250;
	int nmubin = 100;
	int nbin = nrbin*nmubin;
	double maxr = 250.;
	double minr = 0.0;
	double binl[nbin];
	for(int in=0;in<nbin;in++) {
		binl[in] = 0;
	}
	char fxi_name[200];
	char s1[200];
	char s2[200];
	sprintf(s1,argv[1]);	
	sprintf(s2,argv[2]);
	strcpy(fxi_name,s1);
	//sprintf(,strcat(fxi_name,"pcadw.dat"));
	sprintf(fxi_name,"%s%spcadw.dat",ddir,s1);
	FILE *fxi;
	printf("%s %s\n",fxi_name,s2);
	fxi=fopen(fxi_name,"r");
	if (fxi == NULL) perror ("Error opening file");
	int Nline=0;
	const int bsz=10000;
	char buf[bsz];
	
	while(fgets(buf,bsz,fxi)!= NULL) {
		int np;
		double cd,cr,sr,sd,zde,w;
		long double zd;
		sscanf(buf,"%lf %lf %lf %lf %Lf %lf ",&sr,&cr,&sd,&cd,&zd,&w);
		lcd[Nline] = cd;
		lcr[Nline] = cr;
		lsr[Nline] = sr;
		lsd[Nline] = sd;
		lw[Nline] = w;
		lzd[Nline] = zd;
		Nline++;
    }
	char st1[200];
	strcpy(st1,s2);
	char f2_name[200];
	//sprintf(f2_name,strcat(st1,"pcadw.dat"));
	sprintf(f2_name,"%s%spcadw.dat",ddir,st1);
	FILE *f2;
	printf("%s \n",f2_name);
	f2=fopen(f2_name,"r");
	if (f2 == NULL) perror ("Error opening file");
	int Nline2=0;
	const int bsz2=10000;
	char buf2[bsz2];
	
	while(fgets(buf2,bsz2,f2)!= NULL) {
		double cd2,cr2,sr2,sd2,zde2,w2;
		long double zd2;
		sscanf(buf2,"%lf %lf %lf %lf %Lf %lf ",&sr2,&cr2,&sd2,&cd2,&zd2,&w2);
		lcd2[Nline2] = cd2;
		lcr2[Nline2] = cr2;
		lsr2[Nline2] = sr2;
		lsd2[Nline2] = sd2;
		lw2[Nline2] = w2;
		lzd2[Nline2] = zd2;
		Nline2++;
		
    }
	printf("%i %i %f %f %f %f %f %f\n",Nline,Nline2,lcd2[10],lcd2[300],lzd2[5],lzd2[500],lcd[0],lcd[Nline2-1]);
	//double be;
    int be;
    double r;
	double rsq;
	double a;
	double wt2;
	double mu;
	double r1;
	double r2;
	double rrad;
	double gtot = 0;
	int bnt;
	int bmu;
	printf("%f\n",maxr);
	double maxrsq = maxr*maxr;
	double minrsq = minr*minr;
	printf("%f %f\n",maxr,maxrsq);
	while (k<Nline)  {
	//while (k<1){
		j = 0;
		while (j<Nline2)  {
		//while (j < 100) {
			wt2 = lw[k]*lw2[j];
			gtot += wt2;			
			a = lcd[k]*lcd2[j]*(lcr[k]*lcr2[j] + lsr[k]*lsr2[j]) + lsd[k]*lsd2[j];
			rsq = lzd[k]*lzd[k]+lzd2[j]*lzd2[j]-2.*lzd[k]*lzd2[j]*a;
			/*;*/
			if (rsq < maxrsq && rsq > minrsq) {
				r = sqrt(rsq);
				r1 = lzd2[j];
				r2 = lzd[k];
				rrad = fabs(r1-r2);
				mu = rrad/r;
				be = (r-minr);
				bmu = 100*mu;
				bnt = be*100+bmu;
                //printf("%d\n",bnt);
				binl[bnt] += wt2;
				//if (r < 20 and mu < .1) {
				//printf("%f %f %f %f %f \n",r,rrad,mu,r1,r2);
				//		}
									
			}
			//gtot += wt2/dde;
			j++;
			N+=1;
			
		}
		k++;
	}
	double tf = time(NULL)-ti;
	printf("took %g\n",tf);
    //printf("got through loop\n");
	//printf("%f\n",gtot);
	char foname[200];
	if (s1==s2) {
		
	printf("%s %s",s1,s2);
	strcat(s1,"2ptdmu.dat");
	
	}
	else {
		strcat(s1,s2);
		strcat(s1,"2ptdmu.dat");
	}

	//sprintf(foname,s1);
	sprintf(foname,"%s%s",odir,s1);
	FILE *fo;
	fo =fopen(foname,"w");
	fprintf(fo,"%lf \n",gtot);
	for (int pi=0; pi<nbin; pi++) {
		fprintf(fo,"%lf \n",binl[pi]);
	}
		//fclose(fxi);		
	//exit(0);
}
