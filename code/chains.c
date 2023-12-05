#include "chains.h"
#include "random.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h> 

Ribbon *ribbon;

Vector *mvec;
Vector *tt;
Vector *ttmid;
Vector *bv;
Vector *nn;
Vector *bb;
Vector *newn;
Vector *newb;
Vector *dispvec;
Vector *rr;
Vector *rrmid;
Vector *nil;
Vector *refrsec;
Vector *refr;

Statistics stat;

long seed;

double kspring;
double mratio;
double temperature;
double force;
double mvert;
double mlj;
double kpla;
double mdiag;
double smod;
double bigp;
double initd;
double mdih;
double mdihsec;
double anguone;
double angutwo;
double tempangone;
double tempangtwo;
double tempangonesec;
double tempangtwosec;
double torque;
double *energynew;
double *psinew;
double *taunew;
double *thetanew;
double *initdnew;
double *anguonenew;
double *angutwonew;
double *tempangonenew;
double *tempangtwonew;
double *tempangonesecnew;
double *tempangtwosecnew;

double *energynewcheck;
double *f_tt;
double *f_bb;
double *f_nn;
double *f_newb;
double *f_newn;
double *f_ttsec;
double *len;
double *f_bbsec;
double *f_nnsec;
double *f_newbsec;
double *f_newnsec;
double *f_bv;
double *f_ttmid;
double *cos_psi;
double *cos_tau;
double *sin_psi;
double *sin_tau;
double *cos_theta;
double *sin_theta;

double tm[3][3];

double get_zeta();
double get_energy();
double getdistmid();
double get_lj();
double get_zetasec();
double wrmid();
double get_zetasecc();
double get_planar();
double get_dihedral();
double get_spring();
double getlinkdbl();
double get_stackint();
double get_stackintdiag();
double get_curvature();
double getrgone();
double getrgtwo();
double getrgthree();
double get_end_to_end();
double dot(Vector *, Vector *);
double get_link_n();
double getlinknew();
double getlinknewsec();
double get_link_nsec();
double get_twist();
double get_writhe();
Vector cproduct(Vector *, Vector *);

Int_Type num_of_iteration;
Int_Type num_of_particles;
 
void end();
void run();
void init();
void debug();
void frenet();
void measure();
void get_cf(); 
void init_random();
void write_shape();
void write_getlinkdbl();
void allocate_memory();
void read_conf(char *);
void write_conf(FILE *);
void write_wrmid(FILE *);
void one_step(Int_Type);
void one_stepinit(Int_Type);
void write_hi(FILE *, Int_Type);
void get_correlation_functions(double *, double *);

int main() 
{
	init();
	run();
	end();
	return EXIT_SUCCESS;
}

void init()
{
	char f_name[32];
	Int_Type cflag;
 
	FILE *fptr;

	if ((fptr = fopen("Input.txt","r")) == NULL){
		printf("Error! opening file");
    	exit(1);
	}
	
	fscanf(fptr,"%u%u%lg%lg%lg%ld%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg", &num_of_particles, &num_of_iteration, &temperature, &force, &torque, &seed, &kspring, &mdiag, &mvert, &kpla, &mdih, &mdihsec, &smod, &mlj, &bigp, &mratio);
	fclose(fptr); 

	stat.ar = 0;	
	stat.zeta = 0;
	stat.link = 0;
	stat.wrmiddle = 0;
	stat.rg1 = 0;
	stat.rg2 = 0;
	stat.rg3 = 0;
	stat.linkg = 0;
	stat.distmid = 0;
	stat.linksec = 0;
	stat.linknew = 0;
	stat.linknewsec = 0;
	stat.twist = 0; 
	stat.writhe = 0;
	stat.energy = 0;
	stat.curvature = 0;
	stat.end_to_end = 0;
	
	tm[0][0] = 0;
	tm[0][1] = 0;
	tm[0][2] = 0;
	tm[1][0] = 0;
	tm[1][1] = 0;
	tm[1][2] = 0;
	tm[2][0] = 0;
	tm[2][1] = 0;
	tm[2][2] = 0;

    allocate_memory();
    init_random();
}

void allocate_memory()
{
  	Int_Type i;

	ribbon = (Ribbon *) malloc(num_of_particles*sizeof(Ribbon));

	cos_psi = calloc(num_of_particles,sizeof(double));
	cos_tau = calloc(num_of_particles,sizeof(double));
	sin_psi = calloc(num_of_particles,sizeof(double));
	sin_tau = calloc(num_of_particles,sizeof(double));
	
	energynew = (double *) calloc(64, sizeof(double));
	psinew = (double *) calloc(64, sizeof(double));
	taunew = (double *) calloc(64, sizeof(double));
	thetanew = (double *) calloc(64, sizeof(double));
	energynewcheck = (double *) calloc(64, sizeof(double));

	initdnew = (double *) calloc(64, sizeof(double));

	anguonenew = (double *) calloc(64, sizeof(double));
	angutwonew = (double *) calloc(64, sizeof(double));

	tempangonenew = (double *) calloc(64, sizeof(double));
	tempangtwonew = (double *) calloc(64, sizeof(double));

	tempangonesecnew = (double *) calloc(64, sizeof(double));
	tempangtwosecnew = (double *) calloc(64, sizeof(double));

	len = (double *) malloc(num_of_particles*sizeof(double));

	f_tt = (double *) malloc(num_of_particles*sizeof(double));
	f_bb = (double *) malloc(num_of_particles*sizeof(double));
	f_nn = (double *) malloc(num_of_particles*sizeof(double));

	f_newb = (double *) malloc(num_of_particles*sizeof(double));
	f_newn = (double *) malloc(num_of_particles*sizeof(double));

	f_ttsec = (double *) malloc(num_of_particles*sizeof(double));
	f_bbsec = (double *) malloc(num_of_particles*sizeof(double));
	f_nnsec = (double *) malloc(num_of_particles*sizeof(double));

	f_newbsec = (double *) malloc(num_of_particles*sizeof(double));
	f_newnsec = (double *) malloc(num_of_particles*sizeof(double));

 	f_bv = (double *) malloc((num_of_particles/2)*sizeof(double));
	f_ttmid = (double *) malloc((num_of_particles/2)*sizeof(double));

	mvec = (Vector *) malloc(num_of_particles*sizeof(Vector));

	tt = (Vector *) malloc(num_of_particles*sizeof(Vector));	
	nn = (Vector *) malloc(num_of_particles*sizeof(Vector));
	bb = (Vector *) malloc(num_of_particles*sizeof(Vector));	
	newn = (Vector *) malloc(num_of_particles*sizeof(Vector));
	newb = (Vector *) malloc(num_of_particles*sizeof(Vector));

	bv = (Vector *) malloc((num_of_particles/2)*sizeof(Vector));
	ttmid = (Vector *) malloc((num_of_particles/2)*sizeof(Vector));
	

	dispvec = (Vector *) malloc(num_of_particles*sizeof(Vector));

	rr = (Vector *) malloc(num_of_particles*sizeof(Vector));
	rrmid = (Vector *) malloc((num_of_particles/2)*sizeof(Vector));
	
	refr = (Vector *) malloc(num_of_particles*sizeof(Vector));
	refrsec = (Vector *) malloc(num_of_particles*sizeof(Vector));
	nil = (Vector *) malloc(num_of_particles*sizeof(Vector));

	stat.tt = (Vector *) calloc(num_of_particles,sizeof(Vector));
	stat.nn = (Vector *) calloc(num_of_particles,sizeof(Vector));
	stat.bb = (Vector *) calloc(num_of_particles,sizeof(Vector));
	
	stat.f_tt = calloc(num_of_particles,sizeof(double));
	stat.f_bb = calloc(num_of_particles,sizeof(double));
	stat.f_nn = calloc(num_of_particles,sizeof(double));
	stat.f_newb = calloc(num_of_particles,sizeof(double));
	stat.f_newn = calloc(num_of_particles,sizeof(double));

	stat.f_ttsec = calloc(num_of_particles,sizeof(double));
	stat.f_bbsec = calloc(num_of_particles,sizeof(double));
	stat.f_nnsec = calloc(num_of_particles,sizeof(double));
	stat.f_newbsec = calloc(num_of_particles,sizeof(double));
	stat.f_newnsec = calloc(num_of_particles,sizeof(double));

	stat.f_ttmid = calloc((num_of_particles/2),sizeof(double));
	stat.f_bv = calloc((num_of_particles/2),sizeof(double));
}

void init_random()
{
	double psi, tau, theta;
	double const alpha = 2.0;

  	Int_Type i;

	len[0] = 0.64;
	len[num_of_particles/2] = 0.64;

 	for (i = 0; i < num_of_particles/2 - 2; i++){
		len[i+1] = 0.64;
		psi = PI*ran2(&seed);
		tau = 2*PI*ran2(&seed);
		theta = sqrt(-2*log(ran2(&seed))*temperature*0.64/alpha);

		ribbon[i].psi = psi;
		ribbon[i].tau = tau;
		ribbon[i].theta = theta;
  	}
	
	for (i = num_of_particles/2; i < num_of_particles - 2; i++){
		len[i+1] = 0.64;
		psi = PI*ran2(&seed);
		tau = 2*PI*ran2(&seed);
		theta = sqrt(-2*log(ran2(&seed))*temperature*0.64/alpha);

		ribbon[i].psi = psi;
		ribbon[i].tau = tau;
		ribbon[i].theta = theta;
	}

	initd = 2.0;

	anguone = 2.0*PI*ran2(&seed);
	angutwo = PI*ran2(&seed);

	tempangone = 2.0*PI*ran2(&seed);
	tempangtwo = PI*ran2(&seed);

	tempangonesec = 2.0*PI*ran2(&seed);
	tempangtwosec = PI*ran2(&seed);
}

double dot(Vector *v1, Vector *v2){

	return v1->x*v2->x+v1->y*v2->y+v1->z*v2->z; 

}

Vector cproduct(Vector *v1, Vector *v2){
	
	Vector v;
	v.x =v1->y*v2->z-v2->y*v1->z;
	v.y =-v1->x*v2->z+v2->x*v1->z;
	v.z =v1->x*v2->y-v2->x*v1->y;

	return v; 
}

void frenet()
{
	double cth, sth, sp, st, cp, ct, psi, tau, theta, normbvlast, normbvinit, normttmidlast;
	Vector initvec;
	Int_Type i;
	
	nn[0].x = 1.;
	nn[0].y = 0.;
	nn[0].z = 0.;
	
	bb[0].x = 0.;
	bb[0].y = 1.;
	bb[0].z = 0.;

	newn[0].x = 1.;
	newn[0].y = 0.;
	newn[0].z = 0.;
	
	newb[0].x = 0.;
	newb[0].y = 1.;
	newb[0].z = 0.;
	
	initvec.x = cos(anguone)*sin(angutwo);
	initvec.y = sin(anguone)*sin(angutwo);	
	initvec.z = cos(angutwo);

	rr[0].x = -1.0*initd*initvec.x/2.0;
	rr[0].y = -1.0*initd*initvec.y/2.0;
	rr[0].z = -1.0*initd*initvec.z/2.0;

	tt[0].x = cos(tempangone)*sin(tempangtwo);	
	tt[0].y = sin(tempangone)*sin(tempangtwo);
	tt[0].z = cos(tempangtwo);

	for (i = 1; i < num_of_particles/2 - 1; i++){
		double normbi, normti, normni, normdisp, normnewn, normnewb;

		psi = ribbon[i-1].psi;
		tau = ribbon[i-1].tau;
		theta = ribbon[i-1].theta;

		sp = sin(psi);
		st = sin(tau);
		sth = sin(theta);
		cp = cos(psi);
		ct = cos(tau);
		cth = cos(theta);

		mvec[i-1].x = ct*sp;
		mvec[i-1].y = st*sp;
		mvec[i-1].z = cp;
	
		dispvec[i-1] = cproduct(&tt[i-1], &mvec[i-1]);
		normdisp = dot(&dispvec[i-1], &dispvec[i-1]);

		dispvec[i-1].x = dispvec[i-1].x/sqrt(normdisp);
		dispvec[i-1].y = dispvec[i-1].y/sqrt(normdisp);
		dispvec[i-1].z = dispvec[i-1].z/sqrt(normdisp);
		
		tt[i].x = cth*tt[i-1].x + sth*dispvec[i-1].x;
		tt[i].y = cth*tt[i-1].y + sth*dispvec[i-1].y;
		tt[i].z = cth*tt[i-1].z + sth*dispvec[i-1].z;
		
		newn[i].x = tt[i].x - tt[i-1].x;
		newn[i].y = tt[i].y - tt[i-1].y;
		newn[i].z = tt[i].z - tt[i-1].z;
		
		normnewn = sqrt(dot(&newn[i], &newn[i]));

		newn[i].x /= normnewn;
		newn[i].y /= normnewn;
		newn[i].z /= normnewn;

		newb[i] =  cproduct(&tt[i], &newn[i]);

		normnewb = sqrt(dot(&newb[i], &newb[i]));

		newb[i].x /= normnewb;
		newb[i].y /= normnewb;
		newb[i].z /= normnewb;
			
		rr[i].x  = rr[i-1].x+len[i-1]*tt[i-1].x;
		rr[i].y  = rr[i-1].y+len[i-1]*tt[i-1].y;
		rr[i].z  = rr[i-1].z+len[i-1]*tt[i-1].z;
		
		refr[i].x = tt[i].x + tt[i-1].x;
		refr[i].y = tt[i].y + tt[i-1].y;
		refr[i].z = tt[i].z + tt[i-1].z;
		
		nn[i].x = nn[i-1].x - 2.0*dot(&refr[i], &nn[i-1])*refr[i].x/(dot(&refr[i], &refr[i]));
		nn[i].y = nn[i-1].y - 2.0*dot(&refr[i], &nn[i-1])*refr[i].y/(dot(&refr[i], &refr[i]));
		nn[i].z = nn[i-1].z - 2.0*dot(&refr[i], &nn[i-1])*refr[i].z/(dot(&refr[i], &refr[i]));
		
		normni = dot(&nn[i], &nn[i]);
		nn[i].x /= sqrt(normni);
		nn[i].y /= sqrt(normni);
		nn[i].z /= sqrt(normni);

		bb[i] = cproduct(&tt[i], &nn[i]);
		
		normbi = dot(&bb[i], &bb[i]);
		bb[i].x /= sqrt(normbi);
		bb[i].y /= sqrt(normbi);
		bb[i].z /= sqrt(normbi);

	}	

	rr[num_of_particles/2 - 1].x  = rr[num_of_particles/2 - 2].x+len[num_of_particles/2 - 2]*tt[num_of_particles/2 - 2].x;
	rr[num_of_particles/2 - 1].y  = rr[num_of_particles/2 - 2].y+len[num_of_particles/2 - 2]*tt[num_of_particles/2 - 2].y;
	rr[num_of_particles/2 - 1].z  = rr[num_of_particles/2 - 2].z+len[num_of_particles/2 - 2]*tt[num_of_particles/2 - 2].z;	

	rr[num_of_particles/2].x  = 1.0*initd*initvec.x/2.0;
	rr[num_of_particles/2].y  = 1.0*initd*initvec.y/2.0;
	rr[num_of_particles/2].z  = 1.0*initd*initvec.z/2.0;

	tt[num_of_particles/2].x = cos(tempangonesec)*sin(tempangtwosec);	
	tt[num_of_particles/2].y = sin(tempangonesec)*sin(tempangtwosec);
	tt[num_of_particles/2].z = cos(tempangtwosec);

	bv[0].x = rr[num_of_particles/2].x - rr[0].x;
	bv[0].y = rr[num_of_particles/2].y - rr[0].y;
	bv[0].z = rr[num_of_particles/2].z - rr[0].z;

	normbvinit = dot(&bv[0], &bv[0]);
		
	bv[0].x /= sqrt(normbvinit);
	bv[0].y /= sqrt(normbvinit);
	bv[0].z /= sqrt(normbvinit);
	
	rrmid[0].x = 0.5*(rr[num_of_particles/2].x + rr[0].x);
	rrmid[0].y = 0.5*(rr[num_of_particles/2].y + rr[0].y);
	rrmid[0].z = 0.5*(rr[num_of_particles/2].z + rr[0].z);

	nn[num_of_particles/2].x = 1.;
	nn[num_of_particles/2].y = 0.;
	nn[num_of_particles/2].z = 0.;

	bb[num_of_particles/2].x = 0.;
	bb[num_of_particles/2].y = 1.;
	bb[num_of_particles/2].z = 0.;

	newn[num_of_particles/2].x = 1.;
	newn[num_of_particles/2].y = 0.;
	newn[num_of_particles/2].z = 0.;

	newb[num_of_particles/2].x = 0.;
	newb[num_of_particles/2].y = 1.;
	newb[num_of_particles/2].z = 0.;


	for (i = num_of_particles/2 + 1; i < num_of_particles - 1; i++){
		double normttmid, normbi, normti, normni, normdisp, normbvi, normnewn, normnewb;

		psi = ribbon[i-1].psi;
		tau = ribbon[i-1].tau;
		theta = ribbon[i-1].theta;

		sp = sin(psi);
		st = sin(tau);
		sth = sin(theta);
		cp = cos(psi);
		ct = cos(tau);
		cth = cos(theta);

		mvec[i-1].x = ct*sp;
		mvec[i-1].y = st*sp;
		mvec[i-1].z = cp;
	
		dispvec[i-1] = cproduct(&tt[i-1], &mvec[i-1]);
		normdisp = dot(&dispvec[i-1], &dispvec[i-1]);

		dispvec[i-1].x = dispvec[i-1].x/sqrt(normdisp);
		dispvec[i-1].y = dispvec[i-1].y/sqrt(normdisp);
		dispvec[i-1].z = dispvec[i-1].z/sqrt(normdisp);
		
		tt[i].x = cth*tt[i-1].x + sth*dispvec[i-1].x;
		tt[i].y = cth*tt[i-1].y + sth*dispvec[i-1].y;
		tt[i].z = cth*tt[i-1].z + sth*dispvec[i-1].z;
		
		newn[i].x = tt[i].x - tt[i-1].x;
		newn[i].y = tt[i].y - tt[i-1].y;
		newn[i].z = tt[i].z - tt[i-1].z;
		
		normnewn = sqrt(dot(&newn[i], &newn[i]));

		newn[i].x /= normnewn;
		newn[i].y /= normnewn;
		newn[i].z /= normnewn;

		newb[i] =  cproduct(&tt[i], &newn[i]);

		normnewb = sqrt(dot(&newb[i], &newb[i]));

		newb[i].x /= normnewb;
		newb[i].y /= normnewb;
		newb[i].z /= normnewb;
		
		rr[i].x  = rr[i-1].x+len[i-1]*tt[i-1].x;
		rr[i].y  = rr[i-1].y+len[i-1]*tt[i-1].y;
		rr[i].z  = rr[i-1].z+len[i-1]*tt[i-1].z;

		bv[i-num_of_particles/2].x = rr[i].x - rr[i - num_of_particles/2].x;
		bv[i-num_of_particles/2].y = rr[i].y - rr[i - num_of_particles/2].y;
		bv[i-num_of_particles/2].z = rr[i].z - rr[i - num_of_particles/2].z;
	
		normbvi = dot(&bv[i-num_of_particles/2], &bv[i-num_of_particles/2]);
		
		bv[i-num_of_particles/2].x /= sqrt(normbvi);
		bv[i-num_of_particles/2].y /= sqrt(normbvi);
		bv[i-num_of_particles/2].z /= sqrt(normbvi);

		rrmid[i-num_of_particles/2].x = 0.5*(rr[i-num_of_particles/2].x + rr[i].x);
		rrmid[i-num_of_particles/2].y = 0.5*(rr[i-num_of_particles/2].y + rr[i].y);
		rrmid[i-num_of_particles/2].z = 0.5*(rr[i-num_of_particles/2].z + rr[i].z);

		ttmid[i-num_of_particles/2-1].x = rrmid[i-num_of_particles/2].x - rrmid[i-num_of_particles/2-1].x;
		ttmid[i-num_of_particles/2-1].y = rrmid[i-num_of_particles/2].y - rrmid[i-num_of_particles/2-1].y;
		ttmid[i-num_of_particles/2-1].z = rrmid[i-num_of_particles/2].z - rrmid[i-num_of_particles/2-1].z;

		normttmid = dot(&ttmid[i-num_of_particles/2-1], &ttmid[i-num_of_particles/2-1]);
		
		ttmid[i-num_of_particles/2-1].x /= sqrt(normttmid);
		ttmid[i-num_of_particles/2-1].y /= sqrt(normttmid);
		ttmid[i-num_of_particles/2-1].z /= sqrt(normttmid);

		refr[i].x = tt[i].x + tt[i-1].x;
		refr[i].y = tt[i].y + tt[i-1].y;
		refr[i].z = tt[i].z + tt[i-1].z;
		
		nn[i].x = nn[i-1].x - 2.0*dot(&refr[i], &nn[i-1])*refr[i].x/(dot(&refr[i], &refr[i]));
		nn[i].y = nn[i-1].y - 2.0*dot(&refr[i], &nn[i-1])*refr[i].y/(dot(&refr[i], &refr[i]));
		nn[i].z = nn[i-1].z - 2.0*dot(&refr[i], &nn[i-1])*refr[i].z/(dot(&refr[i], &refr[i]));
		
		normni = dot(&nn[i], &nn[i]);
		nn[i].x /= sqrt(normni);
		nn[i].y /= sqrt(normni);
		nn[i].z /= sqrt(normni);

		bb[i] = cproduct(&tt[i], &nn[i]);
		
		normbi = dot(&bb[i], &bb[i]);
		bb[i].x /= sqrt(normbi);
		bb[i].y /= sqrt(normbi);
		bb[i].z /= sqrt(normbi);

	}	
	
	rr[num_of_particles - 1].x  = rr[num_of_particles - 2].x+len[num_of_particles - 2]*tt[num_of_particles - 2].x;
	rr[num_of_particles - 1].y  = rr[num_of_particles - 2].y+len[num_of_particles - 2]*tt[num_of_particles - 2].y;
	rr[num_of_particles - 1].z  = rr[num_of_particles - 2].z+len[num_of_particles - 2]*tt[num_of_particles - 2].z;

	rrmid[num_of_particles/2-1].x = 0.5*(rr[num_of_particles - 1].x + rr[num_of_particles/2 - 1].x);
	rrmid[num_of_particles/2-1].y = 0.5*(rr[num_of_particles - 1].y + rr[num_of_particles/2 - 1].y);
	rrmid[num_of_particles/2-1].z = 0.5*(rr[num_of_particles - 1].z + rr[num_of_particles/2 - 1].z);

	bv[num_of_particles/2 - 1].x = rr[num_of_particles-1].x - rr[num_of_particles/2-1].x;
	bv[num_of_particles/2 - 1].y = rr[num_of_particles-1].y - rr[num_of_particles/2-1].y;
	bv[num_of_particles/2 - 1].z = rr[num_of_particles-1].z - rr[num_of_particles/2-1].z;

	normbvlast = dot(&bv[num_of_particles/2 - 1], &bv[num_of_particles/2 - 1]);

	bv[num_of_particles/2 - 1].x /= sqrt(normbvlast);
	bv[num_of_particles/2 - 1].y /= sqrt(normbvlast);
	bv[num_of_particles/2 - 1].z /= sqrt(normbvlast);

	ttmid[num_of_particles/2-2].x = rrmid[num_of_particles/2-1].x - rrmid[num_of_particles/2-2].x;
	ttmid[num_of_particles/2-2].y = rrmid[num_of_particles/2-1].y - rrmid[num_of_particles/2-2].y;
	ttmid[num_of_particles/2-2].z = rrmid[num_of_particles/2-1].z - rrmid[num_of_particles/2-2].z;

	normttmidlast = dot(&ttmid[num_of_particles/2-2], &ttmid[num_of_particles/2-2]);
		
	ttmid[num_of_particles/2-2].x /= sqrt(normttmidlast);
	ttmid[num_of_particles/2-2].y /= sqrt(normttmidlast);
	ttmid[num_of_particles/2-2].z /= sqrt(normttmidlast);

}

double get_energy()
{	
	double cp, ct, en = 0;
	double const alpha = 1.0/5.0;
	
	Int_Type i;

	for (i = 1; i < num_of_particles/2 - 1; i++){
		
		en += (ribbon[i-1].theta)*(ribbon[i-1].theta) + (ribbon[num_of_particles/2 + i-1].theta)*(ribbon[num_of_particles/2 + i-1].theta);
	}

	for (i = num_of_particles/2+1; i < num_of_particles-1; i++){
		
		en += (ribbon[i-1].theta)*(ribbon[i-1].theta) + (ribbon[num_of_particles/2 + i-1].theta)*(ribbon[num_of_particles/2 + i-1].theta);
	}

	
	if (force){
		en -= get_zetasec();
	}
	
	if (torque){
		en -= 2*PI*torque*(get_link_n() + get_link_nsec());
		
	}
	
	return en;
}


double getdistmid()
{
	double dmid = 0;
	Int_Type i;
	
	frenet();

	for (i = 1; i < num_of_particles/2; i++){
		dmid += sqrt((rrmid[i].x-rrmid[i-1].x)*(rrmid[i].x-rrmid[i-1].x) + (rrmid[i].y-rrmid[i-1].y)*(rrmid[i].y-rrmid[i-1].y) + (rrmid[i].z-rrmid[i-1].z)*(rrmid[i].z-rrmid[i-1].z));
	}

	dmid = dmid/(double)(num_of_particles/2-1);
 
	return dmid;	
}


double getlinkdbl()
{
	double lkdbl = 0;
	Int_Type i,j;
	
	frenet();

	for (i = 1; i < num_of_particles/2; i++){
		for (j = num_of_particles/2+1; j < num_of_particles; j++){
			double dprr, normrr;
			Vector rrdiffstrand, rrdiffseg, rrsecdiffseg, cprr;
			rrdiffstrand.x = rr[i].x-rr[j].x;
			rrdiffstrand.y = rr[i].y-rr[j].y;
			rrdiffstrand.z = rr[i].z-rr[j].z; 

			normrr = sqrt(dot(&rrdiffstrand, &rrdiffstrand));
			rrdiffseg.x = rr[i].x-rr[i-1].x;
			rrdiffseg.y = rr[i].y-rr[i-1].y;
			rrdiffseg.z = rr[i].z-rr[i-1].z;

			rrsecdiffseg.x = rr[j].x-rr[j-1].x;
			rrsecdiffseg.y = rr[j].y-rr[j-1].y;
			rrsecdiffseg.z = rr[j].z-rr[j-1].z;

			cprr = cproduct(&rrdiffseg, &rrsecdiffseg);
			
			dprr = dot(&rrdiffstrand, &cprr);

			lkdbl += dprr/(normrr*normrr*normrr);
		}
	}

	lkdbl = lkdbl/(4*PI);
 
	return lkdbl;	
}

double getrgone()
{
	double rgone = 0;
	Vector rmean, ri;
	Int_Type i,j;
	
	frenet();

	for (i = 0; i < num_of_particles/2; i++){
		rmean.x += rr[i].x;
		rmean.y += rr[i].y;
		rmean.z += rr[i].z;
	}

	rmean.x = rmean.x/(double)(num_of_particles/2);
	rmean.y = rmean.y/(double)(num_of_particles/2);
	rmean.z = rmean.z/(double)(num_of_particles/2);
	
	for (j = 0; j < num_of_particles/2; j++){
		ri.x = rr[i].x - rmean.x;
		ri.y = rr[i].y - rmean.y;
		ri.z = rr[i].z - rmean.z;

		rgone += dot(&ri, &ri);
	}

	rgone = rgone/(double)(num_of_particles/2);
 
	return rgone;	
}

double getrgtwo()
{
	double rgtwo = 0;
	Vector rmeantwo, ri;
	Int_Type i,j;
	
	frenet();

	for (i = num_of_particles/2; i < num_of_particles; i++){
		rmeantwo.x += rr[i].x;
		rmeantwo.y += rr[i].y;
		rmeantwo.z += rr[i].z;
	}

	rmeantwo.x = rmeantwo.x/(double)(num_of_particles/2);
	rmeantwo.y = rmeantwo.y/(double)(num_of_particles/2);
	rmeantwo.z = rmeantwo.z/(double)(num_of_particles/2);

	for (j = num_of_particles/2; j < num_of_particles; j++){
		ri.x = rr[i].x - rmeantwo.x;
		ri.y = rr[i].y - rmeantwo.y;
		ri.z = rr[i].z - rmeantwo.z;

		rgtwo += dot(&ri, &ri);
	}

	rgtwo = rgtwo/(double)(num_of_particles/2);
 
	return rgtwo;	
}

double getrgthree()
{
	double rgthree = 0;
	Vector rmeanthree, ri;
	Int_Type i,j;
	
	frenet();

	for (i = 0; i < num_of_particles; i++){
		rmeanthree.x += rr[i].x;
		rmeanthree.y += rr[i].y;
		rmeanthree.z += rr[i].z;
	}

	rmeanthree.x = rmeanthree.x/(double)(num_of_particles);
	rmeanthree.y = rmeanthree.y/(double)(num_of_particles);
	rmeanthree.z = rmeanthree.z/(double)(num_of_particles);
	
	for (j = 0; j < num_of_particles; j++){
		ri.x = rr[i].x - rmeanthree.x;
		ri.y = rr[i].y - rmeanthree.y;
		ri.z = rr[i].z - rmeanthree.z;

		rgthree += dot(&ri, &ri);
	}

	rgthree = rgthree/(double)(num_of_particles);
 
	return rgthree;	
}


double getlinknew()
{
	Vector dummy,r1,r2,r3,r4,r13,r14,r23,r24,n1,n2,n3,n4;
	double Omega = 0;
	double norm1, norm2, norm3, norm4, ti, lknew = 0, Wr = 0, Tw = 0;
	Int_Type i,j;
	
	frenet();

	ti = dot(&newb[0], &newn[1]);
	Tw += ti;

	for (i = 2; i < num_of_particles/2-2; i++){
		ti = dot(&newb[i-1], &newn[i]);
		Tw += ti;

		for (j = 1; j < i; j++){
			if ((i-j) > 1){

				r1=rr[i];
				r2=rr[i+1];
				r3=rr[j];
				r4=rr[j+1];

				dummy = cproduct(&tt[j],&tt[i]);

				r13.x=r3.x-r1.x;
				r13.y=r3.y-r1.y;
				r13.z=r3.z-r1.z;

				r14.x=r4.x-r1.x;
				r14.y=r4.y-r1.y;
				r14.z=r4.z-r1.z;

				r23.x=r3.x-r2.x;
				r23.y=r3.y-r2.y;
				r23.z=r3.z-r2.z;

				r24.x=r4.x-r2.x;
				r24.y=r4.y-r2.y;
				r24.z=r4.z-r2.z;

				n1 = cproduct(&r13, &r14);
				n2 = cproduct(&r14, &r24);
				n3 = cproduct(&r24, &r23);
				n4 = cproduct(&r23, &r13);
				
				norm1 = sqrt(dot(&n1, &n1)) + 0.00001;
				norm2 = sqrt(dot(&n2, &n2)) + 0.00001;
				norm3 = sqrt(dot(&n3, &n3)) + 0.00001;
				norm4 = sqrt(dot(&n4, &n4)) + 0.00001;
			
				n1.x = n1.x/norm1;
				n1.y = n1.y/norm1;
				n1.z = n1.z/norm1;

				n2.x = n2.x/norm2;
				n2.y = n2.y/norm2;
				n2.z = n2.z/norm2;

				n3.x = n3.x/norm3;
				n3.y = n3.y/norm3;
				n3.z = n3.z/norm3;

				n4.x = n4.x/norm4;
				n4.y = n4.y/norm4;
				n4.z = n4.z/norm4;

				Omega += copysign(1.0, dot(&dummy, &r13))*(asin(dot(&n1, &n2))+ asin(dot(&n2, &n3))+ asin(dot(&n3, &n4))+ asin(dot(&n4, &n1)));
			}
		}	
	}
 
	Wr = Omega;
	lknew = Tw + Wr;
	lknew /= 2*PI;
 
	return lknew;	
}

double getlinknewsec()
{
	Vector dummy,r1,r2,r3,r4,r13,r14,r23,r24,n1,n2,n3,n4;
	double Omega = 0;
	double norm1, norm2, norm3, norm4, ti, lknewsec = 0, Wr = 0, Tw = 0;
	Int_Type i,j;
	
	frenet();

	ti = dot(&newb[num_of_particles/2], &newn[num_of_particles/2+1]);
	Tw += ti;

	for (i = num_of_particles/2+2; i < num_of_particles-2; i++){
		ti = dot(&newb[i-1], &newn[i]);
		Tw += ti;

		for (j = num_of_particles/2+1; j < i; j++){
			if ((i-j) > 1){

				r1=rr[i];
				r2=rr[i+1];
				r3=rr[j];
				r4=rr[j+1];

				dummy = cproduct(&tt[j],&tt[i]);

				r13.x=r3.x-r1.x;
				r13.y=r3.y-r1.y;
				r13.z=r3.z-r1.z;

				r14.x=r4.x-r1.x;
				r14.y=r4.y-r1.y;
				r14.z=r4.z-r1.z;

				r23.x=r3.x-r2.x;
				r23.y=r3.y-r2.y;
				r23.z=r3.z-r2.z;

				r24.x=r4.x-r2.x;
				r24.y=r4.y-r2.y;
				r24.z=r4.z-r2.z;

				n1 = cproduct(&r13, &r14);
				n2 = cproduct(&r14, &r24);
				n3 = cproduct(&r24, &r23);
				n4 = cproduct(&r23, &r13);
				
				norm1 = sqrt(dot(&n1, &n1)) + 0.00001;
				norm2 = sqrt(dot(&n2, &n2)) + 0.00001;
				norm3 = sqrt(dot(&n3, &n3)) + 0.00001;
				norm4 = sqrt(dot(&n4, &n4)) + 0.00001;
			
				n1.x = n1.x/norm1;
				n1.y = n1.y/norm1;
				n1.z = n1.z/norm1;

				n2.x = n2.x/norm2;
				n2.y = n2.y/norm2;
				n2.z = n2.z/norm2;

				n3.x = n3.x/norm3;
				n3.y = n3.y/norm3;
				n3.z = n3.z/norm3;

				n4.x = n4.x/norm4;
				n4.y = n4.y/norm4;
				n4.z = n4.z/norm4;

				Omega += copysign(1.0, dot(&dummy, &r13))*(asin(dot(&n1, &n2))+ asin(dot(&n2, &n3))+ asin(dot(&n3, &n4))+ asin(dot(&n4, &n1)));
			}
		}	
	}
 
	Wr = Omega;
	lknewsec = Tw + Wr;
	lknewsec /= 2*PI;
 
	return lknewsec;	
}


double get_link_n()
{
	Vector dummy,r1,r2,r3,r4,r13,r14,r23,r24,n1,n2,n3,n4;
	double Omega = 0;
	double norm1, norm2, norm3, norm4, ti, lk=0, Wr =0, Tw = 0;
	Int_Type i,j;
	
	frenet();

	ti = - (bb[1].x-bb[0].x)*nn[1].x - (bb[1].y-bb[0].y)*nn[1].y - (bb[1].z-bb[0].z)*nn[1].z;
	Tw += ti;

	for (i = 2; i < num_of_particles/2-2; i++){

		ti = - (bb[i].x-bb[i-1].x)*nn[i].x - (bb[i].y-bb[i-1].y)*nn[i].y - (bb[i].z-bb[i-1].z)*nn[i].z;
		Tw += ti;

		for (j = 1; j < i; j++){
			if ((i-j) > 1){

				r1=rr[i];
				r2=rr[i+1];
				r3=rr[j];
				r4=rr[j+1];

				dummy = cproduct(&tt[j],&tt[i]);

				r13.x=r3.x-r1.x;
				r13.y=r3.y-r1.y;
				r13.z=r3.z-r1.z;

				r14.x=r4.x-r1.x;
				r14.y=r4.y-r1.y;
				r14.z=r4.z-r1.z;

				r23.x=r3.x-r2.x;
				r23.y=r3.y-r2.y;
				r23.z=r3.z-r2.z;

				r24.x=r4.x-r2.x;
				r24.y=r4.y-r2.y;
				r24.z=r4.z-r2.z;

				n1 = cproduct(&r13, &r14);
				n2 = cproduct(&r14, &r24);
				n3 = cproduct(&r24, &r23);
				n4 = cproduct(&r23, &r13);
				
				norm1 = sqrt(dot(&n1, &n1)) + 0.00001;
				norm2 = sqrt(dot(&n2, &n2)) + 0.00001;
				norm3 = sqrt(dot(&n3, &n3)) + 0.00001;
				norm4 = sqrt(dot(&n4, &n4)) + 0.00001;
			
				n1.x = n1.x/norm1;
				n1.y = n1.y/norm1;
				n1.z = n1.z/norm1;

				n2.x = n2.x/norm2;
				n2.y = n2.y/norm2;
				n2.z = n2.z/norm2;

				n3.x = n3.x/norm3;
				n3.y = n3.y/norm3;
				n3.z = n3.z/norm3;

				n4.x = n4.x/norm4;
				n4.y = n4.y/norm4;
				n4.z = n4.z/norm4;

				Omega += copysign(1.0, dot(&dummy, &r13))*(asin(dot(&n1, &n2))+ asin(dot(&n2, &n3))+ asin(dot(&n3, &n4))+ asin(dot(&n4, &n1)));
			}
		}	
	}
 
	Wr = Omega;
	lk = Tw + Wr;
	lk /= 2*PI;
 
	return lk;	
}

double wrmid()
{
	Vector dummy, r1, r2, r3, r4, r13, r14, r23, r24, n1, n2, n3, n4, ttmidi, ttmidj;
	double Omega = 0;
	double norm1, norm2, norm3, norm4, normttmidi, normttmidj, Wr = 0;
	Int_Type i, j;
	
	frenet();

	for (i = 2; i < num_of_particles/2 - 1; i++){
		for (j = 0; j < i; j++){
			if ((i-j) > 1){
				r1.x = 0.5*(rr[i].x + rr[num_of_particles/2+i].x);
				r1.y = 0.5*(rr[i].y + rr[num_of_particles/2+i].y);
				r1.z = 0.5*(rr[i].z + rr[num_of_particles/2+i].z);

				r2.x = 0.5*(rr[i+1].x + rr[num_of_particles/2+i+1].x);
				r2.y = 0.5*(rr[i+1].y + rr[num_of_particles/2+i+1].y);
				r2.z = 0.5*(rr[i+1].z + rr[num_of_particles/2+i+1].z);

				r3.x = 0.5*(rr[j].x + rr[num_of_particles/2+j].x);
				r3.y = 0.5*(rr[j].y + rr[num_of_particles/2+j].y);
				r3.z = 0.5*(rr[j].z + rr[num_of_particles/2+j].z);

				r4.x = 0.5*(rr[j+1].x + rr[num_of_particles/2+j+1].x);
				r4.y = 0.5*(rr[j+1].y + rr[num_of_particles/2+j+1].y);
				r4.z = 0.5*(rr[j+1].z + rr[num_of_particles/2+j+1].z);

				ttmidi.x = r2.x-r1.x;
				ttmidi.x = r2.y-r1.y;
				ttmidi.x = r2.z-r1.z;

				normttmidi = sqrt(dot(&ttmidi, &ttmidi));

				ttmidi.x /= normttmidi;
				ttmidi.y /= normttmidi;
				ttmidi.z /= normttmidi;

				ttmidj.x = r4.x-r3.x;
				ttmidj.x = r4.y-r3.y;
				ttmidj.x = r4.z-r3.z;

				normttmidj = sqrt(dot(&ttmidj, &ttmidj));

				ttmidj.x /= normttmidj;
				ttmidj.y /= normttmidj;
				ttmidj.z /= normttmidj;

				dummy = cproduct(&ttmidj, &ttmidi);

				r13.x = r3.x-r1.x;
				r13.y = r3.y-r1.y;
				r13.z = r3.z-r1.z;

				r14.x = r4.x-r1.x;
				r14.y = r4.y-r1.y;
				r14.z = r4.z-r1.z;

				r23.x = r3.x-r2.x;
				r23.y = r3.y-r2.y;
				r23.z = r3.z-r2.z;

				r24.x = r4.x-r2.x;
				r24.y = r4.y-r2.y;
				r24.z = r4.z-r2.z;

				n1 = cproduct(&r13, &r14);
				n2 = cproduct(&r14, &r24);
				n3 = cproduct(&r24, &r23);
				n4 = cproduct(&r23, &r13);
				
				norm1 = sqrt(dot(&n1,&n1)) + 0.0001;
				norm2 = sqrt(dot(&n2,&n2)) + 0.0001;
				norm3 = sqrt(dot(&n3,&n3)) + 0.0001;
				norm4 = sqrt(dot(&n4,&n4)) + 0.0001;
			
				n1.x=n1.x/norm1;
				n1.y=n1.y/norm1;
				n1.z=n1.z/norm1;

				n2.x=n2.x/norm2;
				n2.y=n2.y/norm2;
				n2.z=n2.z/norm2;

				n3.x=n3.x/norm3;
				n3.y=n3.y/norm3;
				n3.z=n3.z/norm3;

				n4.x=n4.x/norm4;
				n4.y=n4.y/norm4;
				n4.z=n4.z/norm4;
			
				Omega += copysign(1.0, dot(&dummy, &r13))*(asin(dot(&n1,&n2)) + asin(dot(&n2,&n3)) + asin(dot(&n3,&n4)) + asin(dot(&n4,&n1)));
			}
		}	 
	}
 	
	Wr = Omega/(2*PI);

	return Wr;	
}

double get_link_nsec()
{
	Vector dummy,r1,r2,r3,r4,r13,r14,r23,r24,n1,n2,n3,n4;
	double Omega = 0;
	double norm1, norm2, norm3, norm4, ti, lksec = 0, Wr = 0, Tw = 0;
	Int_Type i,j;
	
	frenet();

	ti = - (bb[1].x-bb[0].x)*nn[1].x - (bb[1].y-bb[0].y)*nn[1].y - (bb[1].z-bb[0].z)*nn[1].z;
	Tw += ti;

	for (i = 2; i < num_of_particles-2; i++){

		ti = - (bb[i].x-bb[i-1].x)*nn[i].x - (bb[i].y-bb[i-1].y)*nn[i].y - (bb[i].z-bb[i-1].z)*nn[i].z;
		Tw += ti;

		for (j = 1; j < i; j++){
			if ((i-j) > 1){

				r1=rr[i];
				r2=rr[i+1];
				r3=rr[j];
				r4=rr[j+1];

				dummy = cproduct(&tt[j],&tt[i]);

				r13.x=r3.x-r1.x;
				r13.y=r3.y-r1.y;
				r13.z=r3.z-r1.z;

				r14.x=r4.x-r1.x;
				r14.y=r4.y-r1.y;
				r14.z=r4.z-r1.z;

				r23.x=r3.x-r2.x;
				r23.y=r3.y-r2.y;
				r23.z=r3.z-r2.z;

				r24.x=r4.x-r2.x;
				r24.y=r4.y-r2.y;
				r24.z=r4.z-r2.z;

				n1 = cproduct(&r13, &r14);
				n2 = cproduct(&r14, &r24);
				n3 = cproduct(&r24, &r23);
				n4 = cproduct(&r23, &r13);
				
				norm1 = sqrt(dot(&n1, &n1)) + 0.00001;
				norm2 = sqrt(dot(&n2, &n2)) + 0.00001;
				norm3 = sqrt(dot(&n3, &n3)) + 0.00001;
				norm4 = sqrt(dot(&n4, &n4)) + 0.00001;
			
				n1.x = n1.x/norm1;
				n1.y = n1.y/norm1;
				n1.z = n1.z/norm1;

				n2.x = n2.x/norm2;
				n2.y = n2.y/norm2;
				n2.z = n2.z/norm2;

				n3.x = n3.x/norm3;
				n3.y = n3.y/norm3;
				n3.z = n3.z/norm3;

				n4.x = n4.x/norm4;
				n4.y = n4.y/norm4;
				n4.z = n4.z/norm4;

				Omega += copysign(1.0, dot(&dummy, &r13))*(asin(dot(&n1, &n2))+ asin(dot(&n2, &n3))+ asin(dot(&n3, &n4))+ asin(dot(&n4, &n1)));
			}
		}	
	}
 
	Wr = Omega;
	lksec = Tw + Wr;
	lksec /= 2*PI;
 
	return lksec;	
}


double get_twist()
{
	double ti;
	double Tw = 0;
	Int_Type i;
	
	frenet();
	
	for (i = 1; i < num_of_particles/2-1; i++){

		ti = dot(&newb[i-1], &newn[i]);
		Tw += ti;
	}
	
	Tw = Tw / (2*PI);
	
	return Tw;
	
}


double get_writhe()
{
	Vector dummy,r1,r2,r3,r4,r13,r14,r23,r24,n1,n2,n3,n4;
	double Omega = 0;
	double norm1, norm2, norm3, norm4, Wr = 0 ;
	Int_Type i,j;
	
	frenet();

	for (i = 2; i < num_of_particles-2; i++){
		for (j = 1; j < i; j++){
			if ((i-j) > 1){

				r1=rr[i];
				r2=rr[i+1];
				r3=rr[j];
				r4=rr[j+1];

				dummy = cproduct(&tt[j],&tt[i]);

				r13.x=r3.x-r1.x;
				r13.y=r3.y-r1.y;
				r13.z=r3.z-r1.z;

				r14.x=r4.x-r1.x;
				r14.y=r4.y-r1.y;
				r14.z=r4.z-r1.z;

				r23.x=r3.x-r2.x;
				r23.y=r3.y-r2.y;
				r23.z=r3.z-r2.z;

				r24.x=r4.x-r2.x;
				r24.y=r4.y-r2.y;
				r24.z=r4.z-r2.z;

				n1 = cproduct(&r13, &r14);
				n2 = cproduct(&r14, &r24);
				n3 = cproduct(&r24, &r23);
				n4 = cproduct(&r23, &r13);
				
				norm1 = sqrt(dot(&n1, &n1)) + 0.00001;
				norm2 = sqrt(dot(&n2, &n2)) + 0.00001;
				norm3 = sqrt(dot(&n3, &n3)) + 0.00001;
				norm4 = sqrt(dot(&n4, &n4)) + 0.00001;
			
				n1.x = n1.x/norm1;
				n1.y = n1.y/norm1;
				n1.z = n1.z/norm1;

				n2.x = n2.x/norm2;
				n2.y = n2.y/norm2;
				n2.z = n2.z/norm2;

				n3.x = n3.x/norm3;
				n3.y = n3.y/norm3;
				n3.z = n3.z/norm3;

				n4.x = n4.x/norm4;
				n4.y = n4.y/norm4;
				n4.z = n4.z/norm4;

				Omega += copysign(1.0, dot(&dummy, &r13))*(asin(dot(&n1, &n2))+ asin(dot(&n2, &n3))+ asin(dot(&n3, &n4))+ asin(dot(&n4, &n1)));
			}
		}	
	}
 
	Wr = Omega;
	Wr = Wr/(2*PI);
	
	return Wr;
	
}

double get_end_to_end()
{
	Int_Type i;
	Vector pt;
	double r2;
	
	frenet();
	
	pt.x = 0;
	pt.y = 0;
	pt.z = 0;
	
	for (i = 0; i < num_of_particles-1; i++){
		pt.x += tt[i].x;
		pt.y += tt[i].y;
		pt.z += tt[i].z;
	}
	
	r2 = pt.x*pt.x+pt.y*pt.y+pt.z*pt.z;

	return r2;
}

double get_spring()
{
	double es = 0.0;
	Vector bonddir;
	Int_Type i;
	
	frenet();
	
	for (i=0; i < num_of_particles/2; i++){
		double sepdist;
		bonddir.x = rr[num_of_particles/2 + i].x - rr[i].x;
		bonddir.y = rr[num_of_particles/2 + i].y - rr[i].y;
		bonddir.z = rr[num_of_particles/2 + i].z - rr[i].z;

		sepdist = sqrt(dot(&bonddir, &bonddir));	
		es +=  kspring*(sepdist - 2.0)*(sepdist - 2.0);	
	}	

	return es;
}

double get_lj(){

	double parg, sqdistij, elj = 0.0;
	const int halfpatch = 4;
	Int_Type i,j;
	frenet();		
			
	for (i=2; i < num_of_particles; i++){
		if (i < num_of_particles/2 && (i-1)%halfpatch == 0){
			for(j = i+1; j < num_of_particles; j++){
				if (j < num_of_particles/2 && (j-1)%halfpatch == 0){
					sqdistij = (rr[i].x-rr[j].x)*(rr[i].x-rr[j].x) + (rr[i].y-rr[j].y)*(rr[i].y-rr[j].y) + (rr[i].z-rr[j].z)*(rr[i].z-rr[j].z);	
					if (sqdistij < pow(2, 1.0/3.0)*16.0){
						parg = 4.0/sqdistij; 
						elj += 4.0*mlj*(pow(parg,6) - pow(parg,3) + 0.25);
					}
				}
				
				else if (j > num_of_particles/2 && (j-num_of_particles/2-1)%halfpatch == 0){
					sqdistij = (rr[i].x-rr[j].x)*(rr[i].x-rr[j].x) + (rr[i].y-rr[j].y)*(rr[i].y-rr[j].y) + (rr[i].z-rr[j].z)*(rr[i].z-rr[j].z);	
					if (sqdistij < pow(2, 1.0/3.0)*16.0){
						parg = 4.0/sqdistij; 
						elj += 4.0*mlj*(pow(parg,6) - pow(parg,3) + 0.25);
					}
				}
			
				else {
					elj += 0.0;
				}

			}
		}	

		else if (i > num_of_particles/2+1 && (i-num_of_particles/2-1)%halfpatch == 0){
			for(j = i+1; j < num_of_particles; j++){
				if ((j-num_of_particles/2-1)%halfpatch == 0){
					sqdistij = (rr[i].x-rr[j].x)*(rr[i].x-rr[j].x) + (rr[i].y-rr[j].y)*(rr[i].y-rr[j].y) + (rr[i].z-rr[j].z)*(rr[i].z-rr[j].z);	
					if (sqdistij < pow(2, 1.0/3.0)*16.0){
						parg = 4.0/sqdistij; 
						elj += 4.0*mlj*(pow(parg,6) - pow(parg,3) + 0.25);
					}
				}
			
				else {
					elj += 0.0;
				}


			}
		}

		else {
			elj += 0.0;

		}	

	}
	
	return elj;

}

double get_planar()
{
	double epla = 0.0;
	Int_Type i;
	
	frenet();
	
	for (i = 1; i < num_of_particles/2; i++){
		double normttmidi, normmn;
		Vector mid, ttmidi, midnext, refv, cnext;

		mid.x = 0.5*(rr[num_of_particles/2 + i-1].x + rr[i-1].x);
        	mid.y = 0.5*(rr[num_of_particles/2 + i-1].y + rr[i-1].y);
        	mid.z = 0.5*(rr[num_of_particles/2 + i-1].z + rr[i-1].z);

		midnext.x = 0.5*(rr[num_of_particles/2 + i].x + rr[i].x);
        	midnext.y = 0.5*(rr[num_of_particles/2 + i].y + rr[i].y);
    		midnext.z = 0.5*(rr[num_of_particles/2 + i].z + rr[i].z);

		cnext.x = rr[num_of_particles/2 + i].x - rr[i].x;
		cnext.y = rr[num_of_particles/2 + i].y - rr[i].y;
		cnext.z = rr[num_of_particles/2 + i].z - rr[i].z;

		normmn = sqrt(dot(&cnext, &cnext));	

		cnext.x /= normmn;
		cnext.y /= normmn;
		cnext.z /= normmn;

		ttmidi.x = midnext.x-mid.x;
		ttmidi.y = midnext.y-mid.y;
		ttmidi.z = midnext.z-mid.z;

		normttmidi = sqrt(dot(&ttmidi, &ttmidi));

		ttmidi.x /= normttmidi;
		ttmidi.y /= normttmidi;
		ttmidi.z /= normttmidi;

		refv.x = 0.0;
		refv.y = 0.0;
		refv.z = 1.0;

		epla += 0.5*kpla*(acos(dot(&ttmidi, &cnext)) - 0.5*PI)*(acos(dot(&ttmidi, &cnext)) - 0.5*PI);		
	}
	
	return epla;
}

double get_dihedral()
{
	double edih = 0.0;
	Int_Type i;
	
	frenet();
	
	for (i = 1; i < num_of_particles/2; i++){
		double dihangle, normn1, normn2, normfb, normef, normea, normtestvec, signcos;
		Vector mid, ea, ef, midnext, fb, n1, n2, testvec;

		mid.x = 0.5*(rr[num_of_particles/2 + i-1].x + rr[i-1].x);
        	mid.y = 0.5*(rr[num_of_particles/2 + i-1].y + rr[i-1].y);
        	mid.z = 0.5*(rr[num_of_particles/2 + i-1].z + rr[i-1].z);

		ea.x = rr[num_of_particles/2 + i-1].x - mid.x;
		ea.y = rr[num_of_particles/2 + i-1].y - mid.y;
		ea.z = rr[num_of_particles/2 + i-1].z - mid.z;

		normea = sqrt(dot(&ea, &ea));

		ea.x /= normea;
		ea.y /= normea;
		ea.z /= normea;

		midnext.x = 0.5*(rr[num_of_particles/2 + i].x + rr[i].x);
        	midnext.y = 0.5*(rr[num_of_particles/2 + i].y + rr[i].y);
    		midnext.z = 0.5*(rr[num_of_particles/2 + i].z + rr[i].z);

		fb.x = rr[num_of_particles/2 + i].x - midnext.x;
		fb.y = rr[num_of_particles/2 + i].y - midnext.y;
		fb.z = rr[num_of_particles/2 + i].z - midnext.z;

		normfb = sqrt(dot(&fb, &fb));

		fb.x /= normfb;
		fb.y /= normfb;
		fb.z /= normfb;

		ef.x = midnext.x-mid.x;
		ef.y = midnext.y-mid.y;
		ef.z = midnext.z-mid.z;

		normef = sqrt(dot(&ef, &ef));

		ef.x /= normef;
		ef.y /= normef;
		ef.z /= normef;

		n1 = cproduct(&ea, &ef);
		normn1 = sqrt(dot(&n1, &n1));

		n1.x /= normn1;
		n1.y /= normn1;
		n1.z /= normn1;

		n2 = cproduct(&fb, &ef);
		normn2 = sqrt(dot(&n2, &n2));

		n2.x /= normn2;
		n2.y /= normn2;
		n2.z /= normn2;

		dihangle = acos(dot(&n1, &n2));

		testvec = cproduct(&n1, &n2);
	
		normtestvec = sqrt(dot(&testvec, &testvec));

		testvec.x /= normtestvec;
		testvec.y /= normtestvec;
		testvec.z /= normtestvec;

		signcos = copysign(1.0, dot(&ef, &testvec));
		
		if (signcos > 0){
			edih += mdih*(dihangle-36.0*PI/180.0)*(dihangle-36.0*PI/180.0);
		}

		else if (signcos < 0){
			edih += mdih*(2.0*PI-dihangle-36.0*PI/180.0)*(2.0*PI-dihangle-36.0*PI/180.0);
		}	

		ea.x = rr[i-1].x - mid.x;
		ea.y = rr[i-1].y - mid.y;
		ea.z = rr[i-1].z - mid.z;

		normea = sqrt(dot(&ea, &ea));

		ea.x /= normea;
		ea.y /= normea;
		ea.z /= normea;

		fb.x = rr[i].x - midnext.x;
		fb.y = rr[i].y - midnext.y;
		fb.z = rr[i].z - midnext.z;

		normfb = sqrt(dot(&fb, &fb));

		fb.x /= normfb;
		fb.y /= normfb;
		fb.z /= normfb;

		n1 = cproduct(&ef, &ea);
		normn1 = sqrt(dot(&n1, &n1));

		n1.x /= normn1;
		n1.y /= normn1;
		n1.z /= normn1;

		n2 = cproduct(&ef, &fb);
		normn2 = sqrt(dot(&n2, &n2));

		n2.x /= normn2;
		n2.y /= normn2;
		n2.z /= normn2;

		dihangle = acos(dot(&n1, &n2));
		edih += mdihsec*(dihangle-36.0*PI/180.0)*(dihangle-36.0*PI/180.0);
		
		}
	
	return edih;
}


double get_stackintdiag()
{
	double esidiag = 0.0;
	Int_Type i;
	
	frenet();
	
	for (i = 1; i < num_of_particles/2-1; i++){

		Vector diagone, diagtwo;
		double diagonedist, diagtwodist;

		diagone.x = rr[i].x - rr[num_of_particles/2 + i-1].x;
		diagone.y = rr[i].y - rr[num_of_particles/2 + i-1].y;
		diagone.z = rr[i].z - rr[num_of_particles/2 + i-1].z;

		diagtwo.x = rr[i].x - rr[num_of_particles/2 + i+1].x;
		diagtwo.y = rr[i].y - rr[num_of_particles/2 + i+1].y;
		diagtwo.z = rr[i].z - rr[num_of_particles/2 + i+1].z;

		diagonedist = sqrt(dot(&diagone, &diagone));
		diagtwodist = sqrt(dot(&diagtwo, &diagtwo));

		esidiag += mdiag*(diagonedist - 1.8)*(diagonedist - 1.8) + mdiag*(diagtwodist - 1.8)*(diagtwodist - 1.8);
	
		diagone.x = rr[num_of_particles/2 + i].x - rr[i-1].x;
		diagone.y = rr[num_of_particles/2 + i].y - rr[i-1].y;
		diagone.z = rr[num_of_particles/2 + i].z - rr[i-1].z;

		diagtwo.x = rr[num_of_particles/2 + i].x - rr[i+1].x;
		diagtwo.y = rr[num_of_particles/2 + i].y - rr[i+1].y;
		diagtwo.z = rr[num_of_particles/2 + i].z - rr[i+1].z;

		diagonedist = sqrt(dot(&diagone, &diagone));
		diagtwodist = sqrt(dot(&diagtwo, &diagtwo));

		esidiag += mdiag*(diagonedist - mratio*1.8)*(diagonedist - mratio*1.8) + mdiag*(diagtwodist - mratio*1.8)*(diagtwodist - mratio*1.8);
				
	}
	
	return esidiag;
}

double get_stackint()
{
	double esi = 0.0;
	Int_Type i;
	
	frenet();
	
	for (i = 1; i < num_of_particles/2; i++){

		Vector ellone, elltwo, crossdir, diffvec;
		double magcrossdirsq, shortestdist;

		ellone.x = rr[num_of_particles/2 + i-1].x - rr[i-1].x;
		ellone.y = rr[num_of_particles/2 + i-1].y - rr[i-1].y;
		ellone.z = rr[num_of_particles/2 + i-1].z - rr[i-1].z;

		elltwo.x = rr[num_of_particles/2 + i].x - rr[i].x;
		elltwo.y = rr[num_of_particles/2 + i].y - rr[i].y;
		elltwo.z = rr[num_of_particles/2 + i].z - rr[i].z;

		crossdir = cproduct(&ellone, &elltwo);	

		magcrossdirsq = dot(&crossdir, &crossdir);

		diffvec.x = rr[i].x - rr[i-1].x;
		diffvec.y = rr[i].y - rr[i-1].y;
		diffvec.z = rr[i].z - rr[i-1].z;

		shortestdist = fabs(dot(&crossdir, &diffvec)/sqrt(magcrossdirsq));

		esi += 0.5*mvert*(shortestdist - 0.34)*(shortestdist - 0.34);	
	}
	
	return esi;
}


double get_zeta()
{
	double z=0;
	Int_Type i;
	
	frenet();
	
	for (i=0; i<num_of_particles/2-1; i++){
		z += tt[i].z;
	}
	
	return z;
}

double get_zetasec()
{
	double ez = 0.0, ezsec = 0.0, etot = 0.0;
	Int_Type i;
	
	frenet();
	
	for (i=0; i<num_of_particles/2-1; i++){
		ez += force*tt[i].z;
	}

	for (i=num_of_particles/2; i<num_of_particles-1; i++){
		ezsec += force*tt[i].z;
	}

	etot = ez + ezsec;
	
	return etot;
}

double get_zetasecc()
{
	double ez = 0.0, normcvec, normcvecprev;
	Int_Type i;
	Vector cvec, cvecprev, vert, vertprev, deltacvec;
	frenet();
	
	for (i=1; i<num_of_particles/2; i++){

		cvec.x = rr[num_of_particles/2 + i].x - rr[i].x;
		cvec.y = rr[num_of_particles/2 + i].y - rr[i].y;
		cvec.z = rr[num_of_particles/2 + i].z - rr[i].z;

		normcvec = dot(&cvec, &cvec);
		
		cvec.x = cvec.x/sqrt(normcvec);
		cvec.y = cvec.y/sqrt(normcvec);
		cvec.z = cvec.z/sqrt(normcvec);

		cvecprev.x = rr[num_of_particles/2 + i - 1].x - rr[i-1].x;
		cvecprev.y = rr[num_of_particles/2 + i - 1].y - rr[i-1].y;
		cvecprev.z = rr[num_of_particles/2 + i - 1].z - rr[i-1].z;

		normcvecprev = dot(&cvecprev, &cvecprev);
		
		cvecprev.x = cvecprev.x/sqrt(normcvecprev);
		cvecprev.y = cvecprev.y/sqrt(normcvecprev);
		cvecprev.z = cvecprev.z/sqrt(normcvecprev);

		vert.x = rr[i].x + 0.5*sqrt(normcvec)*cvec.x;
		vert.y = rr[i].x + 0.5*sqrt(normcvec)*cvec.y;
		vert.z = rr[i].x + 0.5*sqrt(normcvec)*cvec.z;

		vertprev.x = rr[i-1].x + 0.5*sqrt(normcvecprev)*cvecprev.x;
		vertprev.y = rr[i-1].x + 0.5*sqrt(normcvecprev)*cvecprev.y;
		vertprev.z = rr[i-1].x + 0.5*sqrt(normcvecprev)*cvecprev.z;

		deltacvec.x = vert.x - vertprev.x;
		deltacvec.y = vert.y - vertprev.y;
		deltacvec.z = vert.z - vertprev.z;

		ez += force*sqrt(dot(&deltacvec, &deltacvec))*deltacvec.z;
	}
	
	return ez;
}



double get_curvature(){
	
	double curvature = 0;
	double ki;
	Int_Type i;
	
	frenet();
	
	for (i = 1; i< num_of_particles - 1; i++){
		curvature += 2*(1-cos(ribbon[i-1].psi));
	}
	
	return curvature/(double) num_of_particles;
	//averaged curvature squared//
}

void get_cf(double *f1, double *f2, double *f3, double *f4, double *f5, double *f6, double *f7, double *f8, double *f9, double *f10, double *f11, double *f12){
	
	int i, j, *nu, *mu;
	
	frenet();

	nu = calloc(num_of_particles/2, sizeof(int));
	mu = calloc(num_of_particles/2, sizeof(int));
	
	for (i  = 0; i < num_of_particles/2; i++){
		f1[i] = 0.;
		f2[i] = 0.;
		f3[i] = 0.;
		f4[i] = 0.;
		f5[i] = 0.;
		f6[i] = 0.;
		f7[i] = 0.;
		f8[i] = 0.;
		f9[i] = 0.;
		f10[i] = 0.;
		f11[i] = 0.;
		f12[i] = 0.;	
	}


	for (i = 0; i < num_of_particles/2-1; i++){
		for(j = 0; j < num_of_particles/2-1; j++){
			f1[abs(i-j)] += dot(&tt[i], &tt[j]);
			f2[abs(i-j)] += dot(&bb[i], &bb[j]);
			f3[abs(i-j)] += dot(&nn[i], &nn[j]);
			f4[abs(i-j)] += dot(&tt[i+num_of_particles/2], &tt[j+num_of_particles/2]);
			f5[abs(i-j)] += dot(&bb[i+num_of_particles/2], &bb[j+num_of_particles/2]);
			f6[abs(i-j)] += dot(&nn[i+num_of_particles/2], &nn[j+num_of_particles/2]);
			f7[abs(i-j)] += dot(&newn[i], &newn[j]);
			f8[abs(i-j)] += dot(&newb[i], &newb[j]);
			f9[abs(i-j)] += dot(&newn[i+num_of_particles/2], &newn[j+num_of_particles/2]);
			f10[abs(i-j)] += dot(&newb[i+num_of_particles/2], &newb[j+num_of_particles/2]);
			f12[abs(i-j)] += dot(&ttmid[i], &ttmid[j]);
			nu[abs(i-j)]++;
		}
	}
	
	for (i = 0; i < num_of_particles/2-1; i++){
		f1[i] /= (double) nu[i];
		f2[i] /= (double) nu[i];
		f3[i] /= (double) nu[i];
		f4[i] /= (double) nu[i];
		f5[i] /= (double) nu[i];
		f6[i] /= (double) nu[i];
		f7[i] /= (double) nu[i];
		f8[i] /= (double) nu[i];
		f9[i] /= (double) nu[i];
		f10[i] /= (double) nu[i];
		f12[i] /= (double) nu[i];
	}	

	free(nu);

	for (i = 0; i < num_of_particles/2; i++){
		for(j = 0; j < num_of_particles/2; j++){
			f11[abs(i-j)] += dot(&bv[i], &bv[j]);
			mu[abs(i-j)]++;
		}
	}

	for (i = 0; i < num_of_particles/2; i++){
		f11[i] /= (double) mu[i];
	}	

	free(mu);

}

void one_step(Int_Type i)
{
	Int_Type p, m, loop, elowestloc, elength;
	double edihold, eplanarold, esiold, esidiagold, espringold, e1, e2, boltzmann;
	double const alpha = 2.0;

	if (i == 0 || i == num_of_particles/2){
		e1 = 0.;
	}

	else if (i != 0 && i != num_of_particles/2){
		e1 = (ribbon[i-1].theta)*(ribbon[i-1].theta);
	}

	e2 = 0.;

	if (mdih > 0.0){
		edihold = get_dihedral();
		e1 += edihold/temperature;

	}

	if (kpla > 0.0){
		eplanarold = get_planar();
		e1 += eplanarold/temperature;
	}


	if (mvert > 0.0){
		esiold = get_stackint();
		e1 += esiold/temperature;
	}

	if (mdiag > 0.0){
		esidiagold = get_stackintdiag();
		e1 += esidiagold/temperature;
	}

	if (kspring > 0.0){
		espringold = get_spring();
		e1 += espringold/temperature;
	}

	for (p = 0; p < 64; p++){
		psinew[p] = PI*ran2(&seed);
    		taunew[p] = 2*PI*ran2(&seed);
		thetanew[p] = sqrt(-2*log(ran2(&seed))*temperature*0.64/alpha);
		energynew[p] = 0.;
		initdnew[p] = 2.0;

		anguonenew[p] = 2.0*PI*ran2(&seed);
		angutwonew[p] = PI*ran2(&seed);

		tempangonenew[p] = 2.0*PI*ran2(&seed);
		tempangtwonew[p] = PI*ran2(&seed);

		tempangonesecnew[p] = 2.0*PI*ran2(&seed);
		tempangtwosecnew[p] = PI*ran2(&seed);


	}
	


	#pragma omp parallel firstprivate(i)
	{
	Vector *mvecp, *dispvecp, *ttp, *rrp;
	Vector initvec;
	double *rpsip, *rtaup, *rthetap;
	Int_Type j, tid;
	
	#pragma omp critical
	rrp = (Vector *) malloc(num_of_particles*sizeof(Vector));
	ttp = (Vector *) malloc(num_of_particles*sizeof(Vector));
	mvecp = (Vector *) malloc(num_of_particles*sizeof(Vector));
	dispvecp = (Vector *) malloc(num_of_particles*sizeof(Vector));

	rpsip = (double *) calloc(num_of_particles, sizeof(double));
	rtaup = (double *) calloc(num_of_particles, sizeof(double));
	rthetap = (double *) calloc(num_of_particles, sizeof(double));
	// frame
	for (j = 1; j < num_of_particles - 1; j++){
		rpsip[j-1] = ribbon[j-1].psi;
		rtaup[j-1] = ribbon[j-1].tau;
		rthetap[j-1] = ribbon[j-1].theta;
	}

	tid = omp_get_thread_num(); 

	if (i != 0 && i != num_of_particles/2){
		rpsip[i-1] = psinew[tid];
		rtaup[i-1] = taunew[tid];
		rthetap[i-1] = thetanew[tid];; 
	}
	

	if (i == 0 || i == num_of_particles/2){
		initvec.x = cos(anguonenew[tid])*sin(angutwonew[tid]);
		initvec.y = sin(anguonenew[tid])*sin(angutwonew[tid]);	
		initvec.z = cos(angutwonew[tid]);

		rrp[0].x = -1.0*initdnew[tid]*initvec.x/2.0;
		rrp[0].y = -1.0*initdnew[tid]*initvec.y/2.0;
		rrp[0].z = -1.0*initdnew[tid]*initvec.z/2.0;

		ttp[0].x = cos(tempangonenew[tid])*sin(tempangtwonew[tid]);	
		ttp[0].y = sin(tempangonenew[tid])*sin(tempangtwonew[tid]);
		ttp[0].z = cos(tempangtwonew[tid]);
	}
	
	else if (i != 0 && i != num_of_particles/2){
		rrp[0].x = rr[0].x;
		rrp[0].y = rr[0].y;
		rrp[0].z = rr[0].z;

		ttp[0].x = tt[0].x;	
		ttp[0].y = tt[0].y;
		ttp[0].z = tt[0].z;
	}
	
	for (j = 1; j < num_of_particles/2 - 1; j++){
		double psi, tau, theta, sp, st, sth, cp, ct, cth, normdisp;
		psi = rpsip[j-1];
		tau = rtaup[j-1];
		theta = rthetap[j-1];

		sp = sin(psi);
		st = sin(tau);
		sth = sin(theta);
		cp = cos(psi);
		ct = cos(tau);
		cth = cos(theta);

		mvecp[j-1].x = ct*sp;
		mvecp[j-1].y = st*sp;
		mvecp[j-1].z = cp;
	
		dispvecp[j-1] = cproduct(&ttp[j-1], &mvecp[j-1]);
		normdisp = dot(&dispvecp[j-1], &dispvecp[j-1]);

		dispvecp[j-1].x = dispvecp[j-1].x/sqrt(normdisp);
		dispvecp[j-1].y = dispvecp[j-1].y/sqrt(normdisp);
		dispvecp[j-1].z = dispvecp[j-1].z/sqrt(normdisp);
		
		ttp[j].x = cth*ttp[j-1].x + sth*dispvecp[j-1].x;
		ttp[j].y = cth*ttp[j-1].y + sth*dispvecp[j-1].y;
		ttp[j].z = cth*ttp[j-1].z + sth*dispvecp[j-1].z;
			
		rrp[j].x  = rrp[j-1].x+len[j-1]*ttp[j-1].x;
		rrp[j].y  = rrp[j-1].y+len[j-1]*ttp[j-1].y;
		rrp[j].z  = rrp[j-1].z+len[j-1]*ttp[j-1].z;
		
	}	

	rrp[num_of_particles/2 - 1].x  = rrp[num_of_particles/2 - 2].x+len[num_of_particles/2 - 2]*ttp[num_of_particles/2 - 2].x;
	rrp[num_of_particles/2 - 1].y  = rrp[num_of_particles/2 - 2].y+len[num_of_particles/2 - 2]*ttp[num_of_particles/2 - 2].y;
	rrp[num_of_particles/2 - 1].z  = rrp[num_of_particles/2 - 2].z+len[num_of_particles/2 - 2]*ttp[num_of_particles/2 - 2].z;

	if (i == 0 || i == num_of_particles/2){
		initvec.x = cos(anguonenew[tid])*sin(angutwonew[tid]);
		initvec.y = sin(anguonenew[tid])*sin(angutwonew[tid]);	
		initvec.z = cos(angutwonew[tid]);

		rrp[num_of_particles/2].x = 1.0*initdnew[tid]*initvec.x/2.0;
		rrp[num_of_particles/2].y = 1.0*initdnew[tid]*initvec.y/2.0;
		rrp[num_of_particles/2].z = 1.0*initdnew[tid]*initvec.z/2.0;

		ttp[num_of_particles/2].x = cos(tempangonesecnew[tid])*sin(tempangtwosecnew[tid]);	
		ttp[num_of_particles/2].y = sin(tempangonesecnew[tid])*sin(tempangtwosecnew[tid]);
		ttp[num_of_particles/2].z = cos(tempangtwosecnew[tid]);
	}
	
	else if (i != 0 && i != num_of_particles/2){
		rrp[num_of_particles/2].x = rr[num_of_particles/2].x;
		rrp[num_of_particles/2].y = rr[num_of_particles/2].y;
		rrp[num_of_particles/2].z = rr[num_of_particles/2].z;

		ttp[num_of_particles/2].x = tt[num_of_particles/2].x;	
		ttp[num_of_particles/2].y = tt[num_of_particles/2].y;
		ttp[num_of_particles/2].z = tt[num_of_particles/2].z;
	}
	

	for (j = num_of_particles/2 + 1; j < num_of_particles - 1; j++){
		double psi, tau, theta, sp, st, sth, cp, ct, cth, normdisp;
		psi = rpsip[j-1];
		tau = rtaup[j-1];
		theta = rthetap[j-1];

		sp = sin(psi);
		st = sin(tau);
		sth = sin(theta);
		cp = cos(psi);
		ct = cos(tau);
		cth = cos(theta);

		mvecp[j-1].x = ct*sp;
		mvecp[j-1].y = st*sp;
		mvecp[j-1].z = cp;
	
		dispvecp[j-1] = cproduct(&ttp[j-1], &mvecp[j-1]);
		normdisp = dot(&dispvecp[j-1], &dispvecp[j-1]);

		dispvecp[j-1].x = dispvecp[j-1].x/sqrt(normdisp);
		dispvecp[j-1].y = dispvecp[j-1].y/sqrt(normdisp);
		dispvecp[j-1].z = dispvecp[j-1].z/sqrt(normdisp);
		
		ttp[j].x = cth*ttp[j-1].x + sth*dispvecp[j-1].x;
		ttp[j].y = cth*ttp[j-1].y + sth*dispvecp[j-1].y;
		ttp[j].z = cth*ttp[j-1].z + sth*dispvecp[j-1].z;
			
		rrp[j].x  = rrp[j-1].x+len[j-1]*ttp[j-1].x;
		rrp[j].y  = rrp[j-1].y+len[j-1]*ttp[j-1].y;
		rrp[j].z  = rrp[j-1].z+len[j-1]*ttp[j-1].z;

	}	
	
	rrp[num_of_particles - 1].x  = rrp[num_of_particles - 2].x+len[num_of_particles - 2]*ttp[num_of_particles - 2].x;
	rrp[num_of_particles - 1].y  = rrp[num_of_particles - 2].y+len[num_of_particles - 2]*ttp[num_of_particles - 2].y;
	rrp[num_of_particles - 1].z  = rrp[num_of_particles - 2].z+len[num_of_particles - 2]*ttp[num_of_particles - 2].z;

    //bending
	if (i != 0 && i != num_of_particles/2){
		energynew[tid] += thetanew[tid]*thetanew[tid]; 
	}
	
	else if (i == 0 || i == num_of_particles/2){
		energynew[tid] += 0.; 
	}
    
    // sb
	for (j = 1; j < num_of_particles/2; j++){
		Vector mid, ea, midnext, fb, ef, n1, n2, testvec;
		double normea, normfb, normef, normn1, normn2, dihangle, normtestvec, signcos;
		mid.x = 0.5*(rrp[num_of_particles/2 + j-1].x + rrp[j-1].x);
       		mid.y = 0.5*(rrp[num_of_particles/2 + j-1].y + rrp[j-1].y);
        	mid.z = 0.5*(rrp[num_of_particles/2 + j-1].z + rrp[j-1].z);

		ea.x = rrp[num_of_particles/2 + j-1].x - mid.x;
		ea.y = rrp[num_of_particles/2 + j-1].y - mid.y;
		ea.z = rrp[num_of_particles/2 + j-1].z - mid.z;

		normea = sqrt(dot(&ea, &ea));

		ea.x /= normea;
		ea.y /= normea;
		ea.z /= normea;

		midnext.x = 0.5*(rrp[num_of_particles/2 + j].x + rrp[j].x);
        	midnext.y = 0.5*(rrp[num_of_particles/2 + j].y + rrp[j].y);
    		midnext.z = 0.5*(rrp[num_of_particles/2 + j].z + rrp[j].z);

		fb.x = rrp[num_of_particles/2 + j].x - midnext.x;
		fb.y = rrp[num_of_particles/2 + j].y - midnext.y;
		fb.z = rrp[num_of_particles/2 + j].z - midnext.z;

		normfb = sqrt(dot(&fb, &fb));

		fb.x /= normfb;
		fb.y /= normfb;
		fb.z /= normfb;

		ef.x = midnext.x-mid.x;
		ef.y = midnext.y-mid.y;
		ef.z = midnext.z-mid.z;

		normef = sqrt(dot(&ef, &ef));

		ef.x /= normef;
		ef.y /= normef;
		ef.z /= normef;

		n1 = cproduct(&ea, &ef);
		normn1 = sqrt(dot(&n1, &n1));

		n1.x /= normn1;
		n1.y /= normn1;
		n1.z /= normn1;

		n2 = cproduct(&fb, &ef);
		normn2 = sqrt(dot(&n2, &n2));

		n2.x /= normn2;
		n2.y /= normn2;
		n2.z /= normn2;

		dihangle = acos(dot(&n1, &n2));

		testvec = cproduct(&n1, &n2);
	
		normtestvec = sqrt(dot(&testvec, &testvec));

		testvec.x /= normtestvec;
		testvec.y /= normtestvec;
		testvec.z /= normtestvec;

		signcos = copysign(1.0, dot(&ef, &testvec));
		
		if (signcos > 0){
			energynew[tid] += mdih*(dihangle-36.0*PI/180.0)*(dihangle-36.0*PI/180.0)/temperature;
		}

		else if (signcos < 0){
			energynew[tid] += mdih*(2.0*PI-dihangle-36.0*PI/180.0)*(2.0*PI-dihangle-36.0*PI/180.0)/temperature;
		}	
	
    }
    // Hbond
    for (j=0; j < num_of_particles/2; j++){
		Vector bonddir;
		double sepdist;

        bonddir.x = rrp[num_of_particles/2 + j].x - rrp[j].x;
        bonddir.y = rrp[num_of_particles/2 + j].y - rrp[j].y;
        bonddir.z = rrp[num_of_particles/2 + j].z - rrp[j].z;

        sepdist = sqrt(dot(&bonddir, &bonddir));	
        energynew[tid] +=  kspring*(sepdist - 2.0)*(sepdist - 2.0)/temperature;	
    }
    // diag
    for (j = 1; j < num_of_particles/2-1; j++){
		Vector diagone, diagtwo;
		double diagonedist, diagtwodist;

		diagone.x = rrp[j].x - rrp[num_of_particles/2 + j-1].x;
		diagone.y = rrp[j].y - rrp[num_of_particles/2 + j-1].y;
		diagone.z = rrp[j].z - rrp[num_of_particles/2 + j-1].z;

		diagtwo.x = rrp[j].x - rrp[num_of_particles/2 + j+1].x;
		diagtwo.y = rrp[j].y - rrp[num_of_particles/2 + j+1].y;
		diagtwo.z = rrp[j].z - rrp[num_of_particles/2 + j+1].z;

		diagonedist = sqrt(dot(&diagone, &diagone));
		diagtwodist = sqrt(dot(&diagtwo, &diagtwo));

		energynew[tid] += mdiag*(diagonedist - 1.8)*(diagonedist - 1.8)/temperature + mdiag*(diagtwodist - 1.8)*(diagtwodist - 1.8)/temperature;
	
		diagone.x = rrp[num_of_particles/2 + j].x - rrp[j-1].x;
		diagone.y = rrp[num_of_particles/2 + j].y - rrp[j-1].y;
		diagone.z = rrp[num_of_particles/2 + j].z - rrp[j-1].z;

		diagtwo.x = rrp[num_of_particles/2 + j].x - rrp[j+1].x;
		diagtwo.y = rrp[num_of_particles/2 + j].y - rrp[j+1].y;
		diagtwo.z = rrp[num_of_particles/2 + j].z - rrp[j+1].z;

		diagonedist = sqrt(dot(&diagone, &diagone));
		diagtwodist = sqrt(dot(&diagtwo, &diagtwo));

		energynew[tid] += mdiag*(diagonedist - mratio*1.8)*(diagonedist - mratio*1.8)/temperature + mdiag*(diagtwodist - mratio*1.8)*(diagtwodist - mratio*1.8)/temperature;		
    }

	
	free(rpsip);
	free(rtaup);
	free(rthetap);

	free(rrp);
	free(ttp);
	free(mvecp);
	free(dispvecp);

    #pragma omp barrier
}

// sort energynew... take the smallest element
    elength = 64;
    elowestloc = 0;
    for(loop = 0; loop < elength; loop++){
        if(energynew[loop] < energynew[elowestloc]){ 
            elowestloc = loop;
        }
    }
// total smallest energy
    e2 += energynew[elowestloc];
//Metropolis
 	boltzmann = exp(e1 - e2);	
	
	if (i != 0 && i != num_of_particles/2){
		if (e2 < e1){
			ribbon[i-1].psi = psinew[elowestloc];
			ribbon[i-1].tau = taunew[elowestloc];
			ribbon[i-1].theta = thetanew[elowestloc];
			stat.ar++;
			return;
		}

		else if (ran2(&seed) < boltzmann){
			ribbon[i-1].psi = psinew[elowestloc];
			ribbon[i-1].tau = taunew[elowestloc];
			ribbon[i-1].theta = thetanew[elowestloc];
			stat.ar++;	
		}
	}

	else if (i == 0 || i == num_of_particles/2){
		if (e2 < e1){
			initd = initdnew[elowestloc];
			anguone = anguonenew[elowestloc];
			angutwo = angutwonew[elowestloc];
			tempangone = tempangonenew[elowestloc];
			tempangtwo = tempangtwonew[elowestloc];
			tempangonesec = tempangonesecnew[elowestloc];
			tempangtwosec = tempangtwosecnew[elowestloc];
			stat.ar++;
			return;
		}

		else if (ran2(&seed) < boltzmann){
			initd = initdnew[elowestloc];
			anguone = anguonenew[elowestloc];
			angutwo = angutwonew[elowestloc];
			tempangone = tempangonenew[elowestloc];
			tempangtwo = tempangtwonew[elowestloc];
			tempangonesec = tempangonesecnew[elowestloc];
			tempangtwosec = tempangtwosecnew[elowestloc];
			stat.ar++;	
		}
	
	}

}



void measure()
{
	double w, psi, tau;	
	Int_Type i;
	
	get_cf(f_tt, f_bb, f_nn, f_ttsec, f_bbsec, f_nnsec, f_newb, f_newn, f_newbsec, f_newnsec, f_bv, f_ttmid);

	w = 2./num_of_iteration;

	psi = ribbon[0].psi;
	tau = ribbon[0].tau;

	tm[0][0] += w*cos(psi);
	tm[0][1] +=-w*sin(psi);
	tm[0][2] += 0;
	tm[1][0] += w*sin(psi)*cos(tau);
	tm[1][1] += w*cos(psi)*cos(tau);
	tm[1][2] +=-w*sin(tau);
	tm[2][0] += w*sin(psi)*sin(tau);
	tm[2][1] += w*cos(psi)*sin(tau);
	tm[2][2] += w*cos(tau);
	
	for (i = 0; i < num_of_particles/2-1; i++){
		stat.tt[i].x += w*tt[i].x;
		stat.tt[i].y += w*tt[i].y;
		stat.tt[i].z += w*tt[i].z;
		stat.bb[i].x += w*bb[i].x;
		stat.bb[i].y += w*bb[i].y;
		stat.bb[i].z += w*bb[i].z;
   		stat.nn[i].x += w*nn[i].x;
		stat.nn[i].y += w*nn[i].y;
		stat.nn[i].z += w*nn[i].z;
		stat.f_tt[i] += w*f_tt[i];
		stat.f_bb[i] += w*f_bb[i];
		stat.f_nn[i] += w*f_nn[i];
		stat.f_ttmid[i] += w*f_ttmid[i];
		stat.f_newb[i] += w*f_newb[i];
		stat.f_newn[i] += w*f_newn[i];
		stat.f_ttsec[i] += w*f_ttsec[i];
		stat.f_bbsec[i] += w*f_bbsec[i];
		stat.f_nnsec[i] += w*f_nnsec[i];
		stat.f_newbsec[i] += w*f_newbsec[i];
		stat.f_newnsec[i] += w*f_newnsec[i];

	}

	for (i = 0; i < num_of_particles/2; i++){
		stat.f_bv[i] += w*f_bv[i];

	}


	for (i = 0; i < num_of_particles-2; i++){

		cos_psi[i] += w*cos(ribbon[i].psi);
		cos_tau[i] += w*cos(ribbon[i].tau);
		sin_psi[i] += w*sin(ribbon[i].psi);
		sin_tau[i] += w*sin(ribbon[i].tau);

	}
	
	stat.end_to_end += w*get_end_to_end();	
	stat.curvature += w*get_curvature();
	stat.energy += w*get_energy();
	stat.zeta += w*get_zeta();
	stat.rg1 += w*getrgone();
	stat.rg2 += w*getrgtwo();
	stat.rg3 += w*getrgthree();
	stat.linkg += w*getlinkdbl();
	stat.distmid += w*getdistmid();
	stat.wrmiddle += w*wrmid();
	stat.link += w*get_link_n();
	stat.linksec += w*get_link_nsec();
	stat.linknew += w*getlinknew();
	stat.linknewsec += w*getlinknewsec();
	stat.twist += w*get_twist();
	stat.writhe += w*get_writhe();
}

void run()
{
	Int_Type t, i, j, period = 1000;
	Vector pt;
	FILE *f_hi;
	FILE *f_ou;

	f_hi = fopen("hi.dat", "w");
	f_ou = fopen("Confper50k.dat", "a");	
	
	for (t = 0; t < num_of_iteration; t++){
		for (i = 0; i < num_of_particles/2 - 1; i++){
			one_step(i);
			one_step(num_of_particles/2+i);
			
		}

		if (t>num_of_iteration/2) measure();

		if((t+1) % period == 0){
	     		printf("i%u ", t);
	     		write_hi(f_hi, t);
	     		fflush(stdout);
		}

		if ((t+1) % 50000 == 0){
			fprintf(f_ou, "st1%u = Graphics3D[{{GrayLevel[0],{Blue, Line[{{%g,%g,%g}", t, rr[0].x, rr[0].y, rr[0].z);	
	for (i = 1; i < num_of_particles/2; i++){
		pt.x = rr[0].x;
		pt.y = rr[0].y;
		pt.z = rr[0].z;
		for (j = 0; j < i; j++){
			pt.x += len[j]*tt[j].x;
			pt.y += len[j]*tt[j].y;
			pt.z += len[j]*tt[j].z;
		}
		fprintf(f_ou,",{%g,%g,%g}", pt.x, pt.y, pt.z);
	}
	
	fprintf(f_ou, "}], Point[{%g,%g,%g}]", rr[0].x, rr[0].y, rr[0].z);
	for (i = 1; i < num_of_particles/2; i++){
		pt.x = rr[0].x;
		pt.y = rr[0].y;
		pt.z = rr[0].z;
		for (j = 0; j < i; j++){
			pt.x += len[j]*tt[j].x;
			pt.y += len[j]*tt[j].y;
			pt.z += len[j]*tt[j].z;
		}
		fprintf(f_ou, ",Point[{%g,%g,%g}]", pt.x, pt.y, pt.z);
	}
	fprintf(f_ou, "}}},Boxed->False];\n");
		
   	fprintf(f_ou, "st2%u = Graphics3D[{{GrayLevel[0], {Red, Line[{{%g,%g,%g}", t, rr[num_of_particles/2].x, rr[num_of_particles/2].y, rr[num_of_particles/2].z);
	
	for (i = num_of_particles/2 + 1; i < num_of_particles; i++){
		pt.x = rr[num_of_particles/2].x;
		pt.y = rr[num_of_particles/2].y;
		pt.z = rr[num_of_particles/2].z;
		for (j = num_of_particles/2; j < i; j++){
			pt.x += len[j]*tt[j].x;
			pt.y += len[j]*tt[j].y;
			pt.z += len[j]*tt[j].z;
		}
		fprintf(f_ou,",{%g,%g,%g}", pt.x, pt.y, pt.z);
	}
	
	fprintf(f_ou, "}], Point[{%g,%g,%g}]", rr[num_of_particles/2].x, rr[num_of_particles/2].y, rr[num_of_particles/2].z);

	for (i = num_of_particles/2 + 1; i < num_of_particles; i++){
		pt.x = rr[num_of_particles/2].x;
		pt.y = rr[num_of_particles/2].y;
		pt.z = rr[num_of_particles/2].z;
		for (j = num_of_particles/2; j < i; j++){
			pt.x += len[j]*tt[j].x;
			pt.y += len[j]*tt[j].y;
			pt.z += len[j]*tt[j].z;
		}
		fprintf(f_ou, ",Point[{%g,%g,%g}]", pt.x, pt.y, pt.z);
	}
	fprintf(f_ou, "}}},Boxed->False];\n\n\n");

		}

	}
	fclose(f_ou);
	fclose(f_hi);
}


void write_hi(FILE *f_ou, Int_Type t)
{
  	double energydiag, energyspring, energydih, energyWLC, linkgauss, writhemid;

	energyWLC = get_energy();
	energydiag = get_stackintdiag();
	energyspring =  get_spring();
	energydih = get_dihedral();
	linkgauss = getlinkdbl();
	writhemid = wrmid();

	fprintf(f_ou,"%u\t%g\t%g\t%g\t%g\t%g\t%g\n", t, linkgauss, writhemid, energyWLC, energyspring, energydiag, energydih);
}

void write_conf(FILE *f_ou)
{
  	double psi, tau, theta;
	Int_Type i;
	
	for (i = 0; i < num_of_particles - 1; i++){
   		psi = fmod(ribbon[i].psi, 2*PI);
   		tau = fmod(ribbon[i].tau, 2*PI);
		theta = fmod(ribbon[i].theta, 2*PI);
	  	fprintf(f_ou,"%g\t%g\t%g\n", psi, tau, theta);
	}
}

void write_getlinkdbl(FILE *f_ou)
{
	Int_Type i,j, numbn;

	frenet();

	for (numbn = 2; numbn <= num_of_particles/2; numbn++){
		double lkdbl = 0.0;

		for (i = 1; i < numbn; i++){
			for (j = num_of_particles/2+1; j < num_of_particles/2+numbn; j++){
				Vector rrdiffstrand, rrdiffseg, rrsecdiffseg, cprr;
				double dprr, normrr;
		
				rrdiffstrand.x = rr[i].x-rr[j].x;
				rrdiffstrand.y = rr[i].y-rr[j].y;
				rrdiffstrand.z = rr[i].z-rr[j].z; 

				normrr = sqrt(dot(&rrdiffstrand, &rrdiffstrand));

				rrdiffseg.x = rr[i].x-rr[i-1].x;
				rrdiffseg.y = rr[i].y-rr[i-1].y;
				rrdiffseg.z = rr[i].z-rr[i-1].z;

				rrsecdiffseg.x = rr[j].x-rr[j-1].x;
				rrsecdiffseg.y = rr[j].y-rr[j-1].y;
				rrsecdiffseg.z = rr[j].z-rr[j-1].z;

				cprr = cproduct(&rrdiffseg, &rrsecdiffseg);
		
				dprr = dot(&rrdiffstrand, &cprr);

				lkdbl += dprr/(normrr*normrr*normrr);
			}
		}

	fprintf(f_ou,"%u\t%g\n", i, lkdbl/(4*PI));
	}

}

void write_wrmid(FILE *f_ou)
{
	Vector dummy, r1, r2, r3, r4, r13, r14, r23, r24, n1, n2, n3, n4, ttmidi, ttmidj;
	double Omega = 0;
	double norm1, norm2, norm3, norm4, normttmidi, normttmidj;
	Int_Type i, j;
	
	frenet();

	fprintf(f_ou,"%g\n", Omega);

	for (i = 2; i < num_of_particles/2 - 1; i++){
		for (j = 0; j < i; j++){
			if ((i-j) > 1){
				r1.x = 0.5*(rr[i].x + rr[num_of_particles/2+i].x);
				r1.y = 0.5*(rr[i].y + rr[num_of_particles/2+i].y);
				r1.z = 0.5*(rr[i].z + rr[num_of_particles/2+i].z);

				r2.x = 0.5*(rr[i+1].x + rr[num_of_particles/2+i+1].x);
				r2.y = 0.5*(rr[i+1].y + rr[num_of_particles/2+i+1].y);
				r2.z = 0.5*(rr[i+1].z + rr[num_of_particles/2+i+1].z);

				r3.x = 0.5*(rr[j].x + rr[num_of_particles/2+j].x);
				r3.y = 0.5*(rr[j].y + rr[num_of_particles/2+j].y);
				r3.z = 0.5*(rr[j].z + rr[num_of_particles/2+j].z);

				r4.x = 0.5*(rr[j+1].x + rr[num_of_particles/2+j+1].x);
				r4.y = 0.5*(rr[j+1].y + rr[num_of_particles/2+j+1].y);
				r4.z = 0.5*(rr[j+1].z + rr[num_of_particles/2+j+1].z);

				ttmidi.x = r2.x-r1.x;
				ttmidi.x = r2.y-r1.y;
				ttmidi.x = r2.z-r1.z;

				normttmidi = sqrt(dot(&ttmidi, &ttmidi));

				ttmidi.x /= normttmidi;
				ttmidi.y /= normttmidi;
				ttmidi.z /= normttmidi;

				ttmidj.x = r4.x-r3.x;
				ttmidj.x = r4.y-r3.y;
				ttmidj.x = r4.z-r3.z;

				normttmidj = sqrt(dot(&ttmidj, &ttmidj));

				ttmidj.x /= normttmidj;
				ttmidj.y /= normttmidj;
				ttmidj.z /= normttmidj;

				dummy = cproduct(&ttmidj, &ttmidi);

				r13.x = r3.x-r1.x;
				r13.y = r3.y-r1.y;
				r13.z = r3.z-r1.z;

				r14.x = r4.x-r1.x;
				r14.y = r4.y-r1.y;
				r14.z = r4.z-r1.z;

				r23.x = r3.x-r2.x;
				r23.y = r3.y-r2.y;
				r23.z = r3.z-r2.z;

				r24.x = r4.x-r2.x;
				r24.y = r4.y-r2.y;
				r24.z = r4.z-r2.z;

				n1 = cproduct(&r13,&r14);
				n2 = cproduct(&r14,&r24);
				n3 = cproduct(&r24,&r23);
				n4 =cproduct(&r23,&r13);
				
				norm1 = sqrt(dot(&n1,&n1))+0.0001;
				norm2 = sqrt(dot(&n2,&n2))+0.0001;
				norm3 = sqrt(dot(&n3,&n3))+0.0001;
				norm4 = sqrt(dot(&n4,&n4))+0.0001;
			
				n1.x = n1.x/norm1;
				n1.y = n1.y/norm1;
				n1.z = n1.z/norm1;

				n2.x = n2.x/norm2;
				n2.y = n2.y/norm2;
				n2.z = n2.z/norm2;

				n3.x = n3.x/norm3;
				n3.y = n3.y/norm3;
				n3.z = n3.z/norm3;

				n4.x = n4.x/norm4;
				n4.y = n4.y/norm4;
				n4.z = n4.z/norm4;
			
				Omega += copysign(1.0, dot(&dummy,&r13))*(asin(dot(&n1,&n2)) + asin(dot(&n2,&n3)) + asin(dot(&n3,&n4)) + asin(dot(&n4,&n1)));
			}
		}	
        
		fprintf(f_ou,"%g\n", Omega/(2*PI));
	}	
}

void end()
{
	char color[32];
	double *cf1;
	Vector pt;
  	FILE *f_ou;

	Int_Type i, j;
	
	frenet();
	
	printf("\n");
	printf("\tEnergy %g\n",stat.energy);
	printf("\tCurvature %g\n",stat.curvature);
	printf("\tEnd to end distance %g\n",stat.end_to_end);
	printf("\tExtension %g\n",stat.zeta);
 	printf("\tLinking %g\n",stat.link);
	printf("\tTwist %g\n", stat.twist);		
	printf("\tWrithe %g\n", stat.writhe);
	printf("\tAcceptance ratio %g\n", stat.ar/(double)(num_of_iteration*num_of_particles));
	printf("\n");

	f_ou = fopen("en.dat", "w");
	fprintf(f_ou,"Number of particles %u\n", num_of_particles);
	fprintf(f_ou,"Number of iteration %u\n", num_of_iteration);
	fprintf(f_ou,"(normalized) Temperature %g\n", temperature);
	fprintf(f_ou,"Hydrogen bond %g\n", kspring);
	fprintf(f_ou,"Diagonal %g\n", mdiag);
	fprintf(f_ou,"Ratio diagonal %g\n", mratio);
	fprintf(f_ou,"Vertical %g\n", mvert);
	fprintf(f_ou,"Dihedral 1 %g\n", mdih);
	fprintf(f_ou,"Dihedral 2 %g\n", mdihsec);
	fprintf(f_ou,"Planar %g\n", kpla);
	fprintf(f_ou,"Stretch modulus %g\n", smod);
	fprintf(f_ou,"L-J %g\n", mlj);
	fprintf(f_ou,"(normalized) P %g\n", bigp);
	fprintf(f_ou,"(normalized) Force %g\n", force);
	fprintf(f_ou,"(normalized) Torque %g\n", torque);
	fprintf(f_ou,"Energy %g\n", stat.energy);
	fprintf(f_ou,"End-to-end distance %g\n", stat.end_to_end);
	fprintf(f_ou,"Relative Extension %g\n", stat.zeta/(double)(num_of_particles));
	fprintf(f_ou,"Gauss link %g\n",stat.linkg);
	fprintf(f_ou,"mid dist %g\n",stat.distmid); 
	fprintf(f_ou,"Writhe (middle) %g\n", stat.wrmiddle); 
	fprintf(f_ou,"RgSq1 %g\n",stat.rg1);
	fprintf(f_ou,"RgSq2 %g\n",stat.rg2);
	fprintf(f_ou,"RgSq3 %g\n",stat.rg3);
	fprintf(f_ou,"Acceptance ratio %g\n", stat.ar/(double)(num_of_iteration*num_of_particles));
	fclose(f_ou);

	f_ou = fopen("bpdist.dat", "w");
	
   	for (i = 0; i < num_of_particles/2; i++){
		Vector diffr;

		diffr.x = rr[num_of_particles/2+i].x-rr[i].x;
		diffr.y = rr[num_of_particles/2+i].y-rr[i].y;
		diffr.z = rr[num_of_particles/2+i].z-rr[i].z;

	fprintf(f_ou, "%u\t%g\n", i, sqrt(dot(&diffr, &diffr)));
	}
	
	fclose(f_ou);

	f_ou = fopen("stackdist.dat", "w");
	
   	for (i = 1; i < num_of_particles/2; i++){
		Vector diffa, crossbv;
		double normcbv;

		diffa.x = rr[i].x-rr[i-1].x;
		diffa.y = rr[i].y-rr[i-1].y;
		diffa.z = rr[i].z-rr[i-1].z;

		crossbv = cproduct(&bv[i-1], &bv[i]);
		normcbv = sqrt(dot(&crossbv, &crossbv));

	fprintf(f_ou, "%u\t%g\n", i, fabs(dot(&diffa, &crossbv)/normcbv));
	}
	
	fclose(f_ou);

	f_ou = fopen("groove1.dat", "w");
	
   	for (i = num_of_particles/2+1 ; i < num_of_particles; i++){
		fprintf(f_ou, "%u\t%g\n", i, sqrt((rr[i].x - rr[1].x)*(rr[i].x - rr[1].x) + (rr[i].y - rr[1].y)*(rr[i].y - rr[1].y) + (rr[i].z - rr[1].z)*(rr[i].z - rr[1].z)));	
	}
	
	fclose(f_ou);

	f_ou = fopen("groove2.dat", "w");
	
   	for (i = 1 ; i < num_of_particles/2; i++){
		fprintf(f_ou, "%u\t%g\n", i, sqrt((rr[i].x - rr[num_of_particles/2+1].x)*(rr[i].x - rr[num_of_particles/2+1].x) + (rr[i].y - rr[num_of_particles/2+1].y)*(rr[i].y - rr[num_of_particles/2+1].y) + (rr[i].z - rr[num_of_particles/2+1].z)*(rr[i].z - rr[num_of_particles/2+1].z)));	
	}
	
	fclose(f_ou);


	f_ou = fopen("midConf.m", "w");	
   	fprintf(f_ou, "Graphics3D[{{GrayLevel[0],{Blue, Line[{{%g,%g,%g}", rrmid[0].x, rrmid[0].y, rrmid[0].z);	
	for (i = 1; i < num_of_particles/2; i++){
		fprintf(f_ou,",{%g,%g,%g}", rrmid[i].x, rrmid[i].y, rrmid[i].z);
	}
	
	fprintf(f_ou, "}], Point[{%g,%g,%g}]", rrmid[0].x, rrmid[0].y, rrmid[0].z);
	for (i = 1; i < num_of_particles/2; i++){
		fprintf(f_ou, ",Point[{%g,%g,%g}]", rrmid[i].x, rrmid[i].y, rrmid[i].z);
	}
	fprintf(f_ou, "}}},Boxed->False]\n");
	fclose(f_ou);

	f_ou = fopen("Conf.m", "w");	
   	fprintf(f_ou, "st1 = Graphics3D[{{GrayLevel[0],{Blue, Line[{{%g,%g,%g}", rr[0].x, rr[0].y, rr[0].z);	
	for (i = 1; i < num_of_particles/2; i++){
		pt.x = rr[0].x;
		pt.y = rr[0].y;
		pt.z = rr[0].z;
		for (j = 0; j < i; j++){
			pt.x += len[j]*tt[j].x;
			pt.y += len[j]*tt[j].y;
			pt.z += len[j]*tt[j].z;
		}
		fprintf(f_ou,",{%g,%g,%g}", pt.x, pt.y, pt.z);
	}
	
	fprintf(f_ou, "}], Point[{%g,%g,%g}]", rr[0].x, rr[0].y, rr[0].z);
	for (i = 1; i < num_of_particles/2; i++){
		pt.x = rr[0].x;
		pt.y = rr[0].y;
		pt.z = rr[0].z;
		for (j = 0; j < i; j++){
			pt.x += len[j]*tt[j].x;
			pt.y += len[j]*tt[j].y;
			pt.z += len[j]*tt[j].z;
		}
		fprintf(f_ou, ",Point[{%g,%g,%g}]", pt.x, pt.y, pt.z);
	}
	fprintf(f_ou, "}}},Boxed->False];\n");
		
   	fprintf(f_ou, "st2 = Graphics3D[{{GrayLevel[0], {Red, Line[{{%g,%g,%g}", rr[num_of_particles/2].x, rr[num_of_particles/2].y, rr[num_of_particles/2].z);
	
	for (i = num_of_particles/2 + 1; i < num_of_particles; i++){
		pt.x = rr[num_of_particles/2].x;
		pt.y = rr[num_of_particles/2].y;
		pt.z = rr[num_of_particles/2].z;
		for (j = num_of_particles/2; j < i; j++){
			pt.x += len[j]*tt[j].x;
			pt.y += len[j]*tt[j].y;
			pt.z += len[j]*tt[j].z;
		}
		fprintf(f_ou,",{%g,%g,%g}", pt.x, pt.y, pt.z);
	}
	
	fprintf(f_ou, "}], Point[{%g,%g,%g}]", rr[num_of_particles/2].x, rr[num_of_particles/2].y, rr[num_of_particles/2].z);

	for (i = num_of_particles/2 + 1; i < num_of_particles; i++){
		pt.x = rr[num_of_particles/2].x;
		pt.y = rr[num_of_particles/2].y;
		pt.z = rr[num_of_particles/2].z;
		for (j = num_of_particles/2; j < i; j++){
			pt.x += len[j]*tt[j].x;
			pt.y += len[j]*tt[j].y;
			pt.z += len[j]*tt[j].z;
		}
		fprintf(f_ou, ",Point[{%g,%g,%g}]", pt.x, pt.y, pt.z);
	}
	fprintf(f_ou, "}}},Boxed->False];\n");
	fprintf(f_ou, "Show[st1,st2]");
	fclose(f_ou);

	f_ou = fopen("Conf1.txt", "w");
	for (i = 0; i < num_of_particles/2; i++){
		fprintf(f_ou,"%g\t%g\t%g\n", rr[i].x, rr[i].y, rr[i].z);
	}
	fclose(f_ou);

	f_ou = fopen("Conf2.txt", "w");
	for (i = num_of_particles/2; i < num_of_particles; i++){
		fprintf(f_ou,"%g\t%g\t%g\n", rr[i].x, rr[i].y, rr[i].z);
	}
	fclose(f_ou);

	f_ou = fopen("cinit.txt", "w");
	fprintf(f_ou,"%g\t%g\t%g\n", initd*cos(anguone)*sin(angutwo)/2.0, initd*sin(anguone)*sin(angutwo)/2.0, initd*cos(angutwo)/2.0);
	fclose(f_ou);


	f_ou = fopen("corrbondvec.dat", "w");
	for (i = 0; i < num_of_particles/2; i++){
		fprintf(f_ou,"%u\t%g\n", i, stat.f_bv[i]);
	}
	fclose(f_ou);

	f_ou = fopen("corrttmid.dat", "w");
	for (i = 0; i < num_of_particles/2-1; i++){
		fprintf(f_ou,"%u\t%g\n", i, stat.f_ttmid[i]);
	}
	fclose(f_ou);

	f_ou = fopen("len.dat", "w");
	for (i = 0; i < num_of_particles-1; i++){
		fprintf(f_ou,"%u\t%g\n", i, len[i]);
	}
	fclose(f_ou);

	f_ou = fopen("checkpar.dat", "w");
	for (i = 0; i < 64; i++){
		fprintf(f_ou,"%g\n", energynew[i]);
	}
	fclose(f_ou);

	f_ou = fopen("cf.dat", "w");
	for (i = 0; i < num_of_particles/2-1; i++){
		fprintf(f_ou,"%u\t%g\t%g\t%g\n", i, stat.f_tt[i], stat.f_bb[i], stat.f_nn[i]);
	}
	fclose(f_ou);

	f_ou = fopen("cfnew.dat", "w");
	for (i = 0; i < num_of_particles/2-1; i++){
		fprintf(f_ou,"%u\t%g\t%g\n", i, stat.f_newb[i], stat.f_newn[i]);
	}
	fclose(f_ou);

	f_ou = fopen("cfnewsec.dat", "w");
	for (i = num_of_particles/2; i < num_of_particles-1; i++){
		fprintf(f_ou,"%u\t%g\t%g\n", i, stat.f_newb[i], stat.f_newn[i]);
	}
	fclose(f_ou);

	
	f_ou = fopen("cfsec.dat", "w");
	for (i = 0; i < num_of_particles/2-1; i++){
		fprintf(f_ou,"%u\t%g\t%g\t%g\n", i, stat.f_ttsec[i], stat.f_bbsec[i], stat.f_nnsec[i]);
	}
	fclose(f_ou);
	

	f_ou = fopen("ou.dat", "w");
	write_conf(f_ou);
	fclose(f_ou);
	
	f_ou = fopen("gausslink.dat", "w");
	write_getlinkdbl(f_ou);
	fclose(f_ou);

	
	f_ou = fopen("wrmiddle.dat", "w");
	write_wrmid(f_ou);
	fclose(f_ou);
	
	f_ou = fopen("gm.dat","w");
	for (i = 0; i < num_of_particles/2-1; i++){
		fprintf(f_ou,"%g\t%g\t%g\n", tt[i].x, tt[i].y, tt[i].z);	
	}
	fclose(f_ou);
	
 	f_ou = fopen("gmbb.dat","w");
	for (i = 0; i < num_of_particles/2-1; i++){
		fprintf(f_ou, "%g\t%g\t%g\n", bb[i].x, bb[i].y, bb[i].z);	
	}
	fclose(f_ou);
 
 	f_ou = fopen("gmnn.dat","w");

	for (i = 0; i < num_of_particles/2-1; i++){
		fprintf(f_ou, "%g\t%g\t%g\n", nn[i].x, nn[i].y, nn[i].z);	
	}
	fclose(f_ou);

	f_ou = fopen("gmsec.dat","w");
	for (i = num_of_particles/2; i < num_of_particles-1; i++){
		fprintf(f_ou,"%g\t%g\t%g\n", tt[i].x, tt[i].y, tt[i].z);	
	}
	fclose(f_ou);
	
 	f_ou = fopen("gmbbsec.dat","w");
	for (i = num_of_particles/2; i < num_of_particles-1; i++){
		fprintf(f_ou, "%g\t%g\t%g\n", bb[i].x, bb[i].y, bb[i].z);	
	}
	fclose(f_ou);
 
 	f_ou = fopen("gmnnsec.dat","w");

	for (i = num_of_particles/2; i < num_of_particles-1; i++){
		fprintf(f_ou, "%g\t%g\t%g\n", nn[i].x, nn[i].y, nn[i].z);	
	}
	fclose(f_ou);


}

void debug()
{
	double dot1, dot2, dot3, norm; 
	Int_Type i;
	
	Vector v1, v2, v3;
	
	free(tt);
	free(nn);
	free(bb);
	free(rr);
	free(ribbon);
	
	num_of_particles=5;
	allocate_memory();
	
	ribbon[0].psi = 1.;
	ribbon[0].tau = 1.;
	
	ribbon[1].psi = 2.;
	ribbon[1].tau = 2.;
	
	ribbon[2].psi = 3.;
	ribbon[2].tau = 3.;

	ribbon[3].psi = 4.;
	ribbon[3].tau = 4.;

	ribbon[4].psi = 5.;
	ribbon[4].tau = 5.;

	v1.x = 1; v2.x = 4;
	v1.y = 2; v2.y = 5;
	v1.z = 3; v2.z = 6;

	v3=cproduct(&v1,&v2);

	norm=sqrt(dot(&v3,&v3));

	printf("v3x %g\n",v3.x);
	printf("v3y %g\n",v3.y);
	printf("v3z %g\n",v3.z);
	printf("norm v3 %g\n",norm);
	printf("asin %g\n",asin(dot(&v3,&v3)/(norm*norm)));
	printf("sign of -3 is %g\n",copysign(1.0,-3));
	printf("dotproduct %g\n",dot(&v1,&v2));
	printf("energy %g\n",get_energy());

	frenet();

	printf("linkingnumber %g\n", get_link_n());

	for (i=0; i<num_of_particles; i++){
		printf("%i) t (%g,%g,%g) n (%g,%g,%g) b (%g,%g,%g) r (%g,%g,%g)\n",
		i,tt[i].x,tt[i].y,tt[i].z,nn[i].x,nn[i].y,nn[i].z,bb[i].x,bb[i].y,bb[i].z,rr[i].x,rr[i].y,rr[i].z);
	}

	for (i=0; i<num_of_particles; i++){
		dot1 = tt[i].x*tt[i].x+tt[i].y*tt[i].y+tt[i].z*tt[i].z;
		dot2 = nn[i].x*nn[i].x+nn[i].y*nn[i].y+nn[i].z*nn[i].z;
		dot3 = bb[i].x*bb[i].x+bb[i].y*bb[i].y+bb[i].z*bb[i].z;
		printf("%i) t.t %g n.n %g b.b %g\n",i,dot1,dot2,dot3);
	}

	
	for (i=0; i<num_of_particles; i++){
		dot1 = tt[i].x*nn[i].x+tt[i].y*nn[i].y+tt[i].z*nn[i].z;
		dot2 = tt[i].x*bb[i].x+tt[i].y*bb[i].y+tt[i].z*bb[i].z;
		dot3 = nn[i].x*bb[i].x+nn[i].y*bb[i].y+nn[i].z*bb[i].z;
		printf("%i) t.n %g t.b %g n.b %g\n",i,dot1,dot2,dot3);
	}

	for (i=0; i<num_of_particles-1; i++){
		dot1 = tt[i].x*tt[i+1].x+tt[i].y*tt[i+1].y+tt[i].z*tt[i+1].z;
		dot2 = bb[i].x*bb[i+1].x+bb[i].y*bb[i+1].y+bb[i].z*bb[i+1].z;
		printf("%i) cos(psi) %g t.t %g cos(tau) %g b.b %g\n",i,cos(ribbon[i].psi),dot1,cos(ribbon[i].tau),dot2);
	}
	
}










