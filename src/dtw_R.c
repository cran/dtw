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



#include <stdlib.h> // for NULL
#include <Rinternals.h> // for SEXP
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "dtw_core.h"


SEXP computeCM_Call(SEXP wm, SEXP lm, SEXP cm, SEXP dir);


#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef R_CallDef[] = {
    CALLDEF(computeCM_Call, 4),
    {NULL, NULL, 0}
};

/* Old triangle fixing, disabled

static R_NativePrimitiveArgType triangle_fixing_l2_t[] = {
        REALSXP, INTSXP, // INOUT
        INTSXP, REALSXP, // IN
        REALSXP          // OUT
};

static R_NativeArgStyle triangle_fixing_l2_s[] = {
                R_ARG_IN_OUT,   R_ARG_IN_OUT,
                R_ARG_IN,       R_ARG_IN,
                R_ARG_OUT
};

static const R_CMethodDef R_CDef[] = {
        {"triangle_fixing_l2", (DL_FUNC) &triangle_fixing_l2,
         5, triangle_fixing_l2_t},
        {NULL, NULL, 0, NULL}
};

*/


void R_init_dtw(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
    R_RegisterCCallable("dtw","computeCM_Call",(DL_FUNC) computeCM_Call);
}




/* --------------------------------------------------
 *
 * Wrapper for .Call, avoids several copies. Returns a list with names
 * "costMatrix" and "directionMatrix"
 */

SEXP computeCM_Call(SEXP wm, 	/* logical */
                    SEXP lm,	/* double */
                    SEXP cm,	/* double */
                    SEXP dir)  	/* double */
{

    /* Get problem size */
    SEXP lm_dim;
    PROTECT(lm_dim = GET_DIM(lm)); /* ---- 1 */
    int *p_lm_dim = INTEGER_POINTER(lm_dim);

    /* Get pattern size */
    SEXP dir_dim;
    PROTECT(dir_dim = GET_DIM(dir)); /* ---- 2 */
    int nsteps=INTEGER_POINTER(dir_dim)[0];

    /* Cost matrix (input+output 1).  */
    SEXP cmo;
    PROTECT(cmo=duplicate(cm)); /* ---- 3 */

    /* Output 2: smo, INTEGER */
    SEXP smo;
    PROTECT(smo=allocMatrix(INTSXP,p_lm_dim[0],p_lm_dim[1]));  /* ---- 4 */

    /* Dispatch to C */
    computeCM(p_lm_dim,
              LOGICAL_POINTER(wm),
              NUMERIC_POINTER(lm),
              &nsteps,
              NUMERIC_POINTER(dir),
              NUMERIC_POINTER(cmo),
              INTEGER_POINTER(smo));

    /* cmo and smo are now set. Put them in a list. From S. Blay,
     http://www.sfu.ca/~sblay/R-C-interface.ppt */
    SEXP list_names;
    PROTECT(list_names = allocVector(STRSXP,2)); /* ---- 5 */
    SET_STRING_ELT(list_names,0,mkChar("costMatrix"));
    SET_STRING_ELT(list_names,1,mkChar("directionMatrix"));

// Creating a list with 2 vector elements:
    SEXP list;
    PROTECT(list = allocVector(VECSXP, 2)); /* ---- 6 */
    SET_VECTOR_ELT(list, 0, cmo);
    SET_VECTOR_ELT(list, 1, smo);
// and attaching the vector names:
    setAttrib(list, R_NamesSymbol, list_names);

    UNPROTECT(6);
    return list;
}

