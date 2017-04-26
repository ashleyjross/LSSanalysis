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
//const char *ddir={"/global/u2/a/ajross/patchypcfiles/"};
const char *ddir={""};

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
    int nangbin = 120;
	int nbin = nangbin;
    double minang = 0*pi/180.;
    double maxang = 6.*pi/180.;
    double binsize = (maxang-minang)/(1.*nbin);
    double angl[nangbin];
    double binedges[nangbin+1];
    for(int bn=0;bn<nangbin+1;bn++) {
        double num = (bn);
        double angle = maxang - binsize*num;
        if (angle < 0) angle = 0;
        double angb = maxang-binsize*num-binsize/2.;
        binedges[bn] = cos(angle);
        printf("%0.8f %0.8f %0.8f\n",angle*180./pi,binedges[bn],cos(angle));
    }

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
		double cd,cr,sr,sd,w;
		sscanf(buf,"%lf %lf %lf %lf %lf ",&sr,&cr,&sd,&cd,&w);
		lcd[Nline] = cd;
		lcr[Nline] = cr;
		lsr[Nline] = sr;
		lsd[Nline] = sd;
		lw[Nline] = w;
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
		double cd2,cr2,sr2,sd2,w2;
		sscanf(buf2,"%lf %lf %lf %lf %lf ",&sr2,&cr2,&sd2,&cd2,&w2);
		lcd2[Nline2] = cd2;
		lcr2[Nline2] = cr2;
		lsr2[Nline2] = sr2;
		lsd2[Nline2] = sd2;
		lw2[Nline2] = w2;
		Nline2++;
		
    }
	printf("%i %i %f %f %f %f \n",Nline,Nline2,lcd2[10],lcd2[300],lcd[0],lcd[Nline2-1]);
	//double be;
    /*for(int bn=0;bn<nangbin+1;bn++) {
        printf("%d %0.8f \n",bn,binedges[bn]);
    }*/
    printf("%0.8f \n",binedges[106]);
    double a;
    double maxa = 0;
	double wt2;
	double gtot = 0;
	int bnt;
    double bth;
	while (k<Nline)  {
	//while (k<2){
        j = 0;
		while (j<Nline2)  {
        //while (j<1)  {
//printf("%i %i %f %f %f %f \n",k,j,lcd2[j],lcr2[j],lcd[k],lcr[k]);
            //while (j < 100) {
            //printf("%0.8f %d\n",binedges[106],k);
            wt2 = lw[k]*lw2[j];

            gtot += wt2;

            a = lcd[k]*lcd2[j]*(lcr[k]*lcr2[j] + lsr[k]*lsr2[j]) + lsd[k]*lsd2[j];
            /*if (lcd[k] == lcd2[j]){
                printf("%i %i %f %f %f %f %f %f %f %f %f %f %f %f %f\n",k,j,lcd2[j],lcr2[j],lcd[k],lcr[k],lw[k],lw2[j],lsr[k],lsr2[j],lsd[k],lsd2[j],wt2,gtot,a);
                //lcd2[j] *= 1.000001;
            }*/
            //bth = binedges[0];
            int ba = -1;
            //while (a > bth) {
            while (a > binedges[ba+1] && ba < nbin) {
                //bth = binedges[ba+2];
                ba++;
            }

            bnt = ba;
            /*if (a > 0.99999962){
                printf("%d %0.8f %d %0.8f %0.8f %0.8f\n",bnt,bth,ba,binedges[ba+1],binedges[106],a);
            }*/

            //printf("%d\n",bnt);
            //printf("%lf \n",binedges[106]);
            if (bnt > -1){
                binl[bnt] += wt2;
                //printf("%lf %d\n",binedges[106],bnt);
            }
				//if (r < 20 and mu < .1) {
				//printf("%f %f %f %f %f \n",r,rrad,mu,r1,r2);
				//		}
									
			j++;
			N+=1;
			
		}
		k++;
	}
	double tf = time(NULL)-ti;
	//printf("%0.8f\n",maxa);
    printf("took %g\n",tf);
    //printf("got through loop\n");
	//printf("%f\n",gtot);
	char foname[200];
	if (s1==s2) {
		
	printf("%s %s",s1,s2);
	strcat(s1,"2ptth.dat");
	
	}
	else {
		strcat(s1,s2);
		strcat(s1,"2ptth.dat");
	}

	//sprintf(foname,s1);
	sprintf(foname,"%s%s",ddir,s1);
	FILE *fo;
	fo =fopen(foname,"w");
	fprintf(fo,"%lf \n",gtot);
	for (int pi=0; pi<nbin; pi++) {
		fprintf(fo,"%lf \n",binl[pi]);
	}
		//fclose(fxi);		
	//exit(0);
}
