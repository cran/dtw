#include <Rinternals.h> // for SEXP

#ifndef _DTW_H
#define _DTW_H

// extern void computeCM(const int *s, const int *wm, const double *lm, const int *nstepsp, const double *dir, double *cm, int *sm);
extern SEXP computeCM_Call(SEXP wm, SEXP lm, SEXP cm, SEXP dir);
extern void triangle_fixing_l2(double *D, int *maxiter_p, const int *n_p, 
                               const double *kappa_p, double *delta_p);

#endif

