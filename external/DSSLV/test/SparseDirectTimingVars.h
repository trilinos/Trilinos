#ifndef SPARSEDIRECTTIMINGVARS_H
#define SPARSEDIRECTTIMINGVARS_H
#ifdef SPARSE_DIRECT_TIMINGS
#include "SparseSolverResult.h" 
#include <fstream>

class SparseDirectTimingVars
{
 private:
 public:
  static SparseSolverResult SS_Result ; 
  static ofstream log_file ; 

} ;
#endif
#endif

