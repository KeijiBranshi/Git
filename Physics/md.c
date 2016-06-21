#include "md.h"

#define nAtom 1000 //# of atoms
#define NMAX 100000 //Max # of atoms
double r[NMAX][3];
double Region[3];
double RegionH[3];

double SignR(double v, double x) {
	if (x > 0)
		return v;
	else
		return -v;
}

void ApplyBoundaryCond() {
	int n,k;
	for (n=0; n<nAtom; ++n)
		for (k=0; k<3; ++k)
			r[n][k] = r[n][k] 
				- SignR(RegionH[k],r[n][k]) 
				- SignR(RegionH[k],r[n][k]-Region[k]);
}