

#ifndef _MLIFPACKWRAP_
#define _MLIFPACKWRAP_

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

  int ML_Ifpack_Gen(ML *ml, int curr_level, int choice, int * options,
		  double * params, void ** Ifpack_Handle);
  int ML_Ifpack_Solve( void * Ifpack_Handle, double * x, double * rhs );
  void ML_Ifpack_Destroy(void * Ifpack_Handle);
  
#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
