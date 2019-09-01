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
#include <R_ext/Rdynload.h>
#include "dtw.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef R_CallDef[] = {
        CALLDEF(computeCM_Call, 4),
        {NULL, NULL, 0}
};

static R_NativePrimitiveArgType triangle_fixing_l2_t[] = {
        REALSXP, INTSXP, // INOUT
        INTSXP, REALSXP, // IN
        REALSXP          // OUT
};

/* For when it will be enabled
        static R_NativeArgStyle triangle_fixing_l2_s[] = {
                R_ARG_IN_OUT,   R_ARG_IN_OUT,
                R_ARG_IN,       R_ARG_IN,
                R_ARG_OUT
        };
 */

static const R_CMethodDef R_CDef[] = {
        {"triangle_fixing_l2", (DL_FUNC) &triangle_fixing_l2, 
         5, triangle_fixing_l2_t},
        {NULL, NULL, 0, NULL}
};



void R_init_dtw(DllInfo *dll)
{
        R_registerRoutines(dll, R_CDef, R_CallDef, NULL, NULL);
        R_useDynamicSymbols(dll, FALSE);
        R_forceSymbols(dll, TRUE);
        R_RegisterCCallable("dtw","computeCM_Call",(DL_FUNC) computeCM_Call);
}

