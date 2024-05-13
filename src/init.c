#include "init.h"

// edit this file to register new model routines with R
// for each model, there must be
// one DECLARATIONS line and one METHODS line.

DECLARATIONS(LBDPwr2)
  
static const R_CallMethodDef callMethods[] = {
  METHODS(LBDPwr2),
  {NULL, NULL, 0}
};

void R_init_SimReassort (DllInfo *info) {
  // Register routines
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info,TRUE);
}
