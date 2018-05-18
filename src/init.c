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
