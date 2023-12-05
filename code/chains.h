#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TYPES_DECLARATION

typedef unsigned int  Int_Type;
typedef unsigned long Long_Type;

typedef struct {
	double x;
	double y;
	double z;
} Vector;

typedef struct {
	double psi;
	double tau;
	double theta;
} Ribbon;

typedef struct {

  	Int_Type ar;

	Vector *tt;
	Vector *bb;
	Vector *bv;
	Vector *nn;
	Vector *newn;
	Vector *newb;
	Vector *dispvec;

	Vector *ttl;
	Vector *nnl;
	Vector *refr;
	
	double *f_tt;
	double *f_ttmid;
	double *f_bb;
	double *f_nn;
	double *f_newb;
	double *f_newn;
	double *f_ttsec;
	double *f_bbsec;
	double *f_nnsec;
	double *f_newbsec;
	double *f_newnsec;
	double *f_bv;
	double zeta;
	double linkg;
 	double link;
	double linknew;
	double linknewsec;
	double linksec; 
	double twist;
	double writhe;
	double rg1;
	double rg2;
	double rg3;
	double link1; 
	double wrmiddle;
	double twist1;
	double writhe1;
	double energy;	
	double distmid;
	double torque; 
	double curvature;
	double end_to_end;
	
} Statistics;

typedef struct {

  Int_Type d;
  Int_Type h;
  Int_Type m;
  Int_Type s;

} CPU_Time;


#endif
#define TYPES_DECLARATION
#define PI 3.141592653589793238462643383279502884
#define SQRT3 1.732050807568877293527446341505872366
