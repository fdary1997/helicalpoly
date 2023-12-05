#include <math.h>
#include "chains.h"
#include "random.h"

double ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;
 
  if (*idum <= 0) { 
    if (-(*idum) < 1) *idum=1; 
    else *idum = -(*idum);
    idum2 = (*idum);
    for (j=NTAB+7; j>=0; j--) {
      k = (*idum)/IQ1;
      *idum = IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum)/IQ1; 
  *idum = IA1*(*idum-k*IQ1)-k*IR1; 
  if (*idum < 0) *idum += IM1; 
  k = idum2/IQ2;
  idum2 = IA2*(idum2-k*IQ2)-k*IR2; 
  if (idum2 < 0) idum2 += IM2;
  j = iy/NDIV; 
  iy = iv[j]-idum2;
  iv[j] = *idum; 
  if (iy < 1) iy += IMM1;
  if ((temp = AM*iy) > RNMX) return RNMX; 
  else return temp;
}

/*******************************************************************/

double gasdev(long *idum, double m, double sigma) 
{
  double fac, rsq, v1, v2, ran2(long *idum);
  static double gset1, gset2;
  static int iset = 0;

  if (*idum < 0) iset=0; 
  if (iset == 0) { 
    do {
      v1 = 2.0*ran2(idum)-1.0;
      v2 = 2.0*ran2(idum)-1.0; 
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0); 
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset1 = m+sigma*v1*fac;
    gset2 = m+sigma*v2*fac;
    iset = 1; 
    return gset2;
  } else { 
    iset = 0; 
    return gset1; 
  }
}
