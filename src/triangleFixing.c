//
// Copyright (c) 2006-2019 of Toni Giorgino
//
// This file is part of the DTW package.
//
// DTW is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DTW is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with DTW.  If not, see <http://www.gnu.org/licenses/>.
//




/* 
 * Triangle fixing algorithm - implementation of algorithm 3.1
 * (Metric_Nearness_L2) in Brickell, J., Dhillon, I., Sra, S., and
 * Tropp, J. (2008). The Metric Nearness Problem. SIAM. J. Matrix
 * Anal. & Appl. 30, 375-396.
 *  
 */


#include <math.h>
#include <R.h>

static int n;

/* Index in 2d matrices:  0 <= i < j < n */
//#define ED(ii,jj) ((jj)*n+(ii))
static inline size_t ED(size_t ii, size_t jj) {
    if(ii<jj) {
	return ((jj)*n+(ii));
    } else {
	return ((ii)*n+(jj));
    }
}


/* Mostly from http://suvrit.de/work/soft/metricn.html metricL2.cc */
double fixOneTriangle(double *D, size_t i, size_t j, size_t k,
		      double *err) {

    double oab, alpha = 0.0;
    double del;

    double ab = D[ED(i,j)];
    double bc = D[ED(j,k)];
    double ca = D[ED(i,k)];
    
    // Save leading edge for abc
    oab = ab;   
    
    alpha = *err;
    del = ab - bc - ca + 3*alpha;
    
    if (del < 0) {
      ab = ab + alpha;
      bc = bc - alpha;
      ca = ca - alpha;
    } else {
      del = del / 3;
      ab  = ab + alpha - del;
      bc  = bc + del - alpha;
      ca  = ca + del - alpha;
    }

    D[ED(i,j)]=ab;
    D[ED(j,k)]=bc;
    D[ED(i,k)]=ca;
    
    /* gsl_matrix_set(d, i, j, ab); */
    /* gsl_matrix_set(d, j, k, bc); */
    /* gsl_matrix_set(d, i, k, ca); */

    *err = oab - ab + *err;
    double echange = fabs(ab - oab);
    return echange;
    
}




/* Index in the triangle inequalities:  0 <= i < j < k < n   following
 * http://machinelearning.wustl.edu/mlpapers/paper_files/NIPS2005_770.pdf
 * some code from http://suvrit.de/work/soft/metricn.html
 *
 * Algorithm L2B did NOT work.
 */

/* Only one triangle of the input matrix D is used (but a square
 * matrix is expected) */
void triangle_fixing_l2(
               /* IN+OUT */
	       double *D,		/* input matrix D, output M */
	       int *maxiter_p,		/* maximum iterations */
	       /* IN */
	       const int *n_p,		/* mtrx dimensions, int */
	       const double *kappa_p,	/* tolerance */
	       /* OUT */
	       double *delta_p		/* final sum of changes */
    ) {

    /* For convenience */
    n=*n_p;

    /* Initialize primal and dual */
    double *z = (double*) S_alloc(n*(n-1)*(n-2)/2,sizeof(double));

    *delta_p = 1.0 + *kappa_p;	/* first iteration */

    /* Convergence test */
    while( (*maxiter_p)-- && *delta_p > *kappa_p ) {

	size_t t=0;
	*delta_p=0.0;
	
	/* Foreach triangle inequality */
	for(size_t i=0; i<n; i++) {
	    for(size_t j=i+1; j<n; j++) {
		for(size_t k=j+1; k<n; k++) {
		    *delta_p += fixOneTriangle(D, i,j,k,z+t);
		    t++;
		    *delta_p += fixOneTriangle(D, j,k,i,z+t);
		    t++;
		    *delta_p += fixOneTriangle(D, k,i,j,z+t);
		    t++;
		}
	    }
	}
	/* delta = sum of changes in the e_ij values (?) */
    }

    for(size_t i=0; i<n; i++) {
	 for(size_t j=i+1; j<n; j++)  {
	      D[i*n+j]   =D[ED(i,j)]; /* symmetrize */
	 }
    }
    
    return;
}


