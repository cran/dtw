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



#ifndef _DTW_H
#define _DTW_H

#include <Rinternals.h> // for SEXP

// extern void computeCM(const int *s, const int *wm, const double *lm, const int *nstepsp, const double *dir, double *cm, int *sm);
extern SEXP computeCM_Call(SEXP wm, SEXP lm, SEXP cm, SEXP dir);
extern void triangle_fixing_l2(double *D, int *maxiter_p, const int *n_p, 
                               const double *kappa_p, double *delta_p);

#endif


