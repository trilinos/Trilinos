/*!
 *  \file ml_ifpack_wrap.h
 *
 *  \brief Interface to the Trilinos package IFPACK.
 *
 *  IFPACK can be used to define ILU-type smoothers for ML. At present,
 *  this interface is still incomplete (but working). The choice of 
 *  factorizations and parameters is hardwired in the code. Plans are to move
 *  towards an Teuchos::ParameterList as input for all parameters.
 *  
 *  \date Last update do Doxygen: 22-Jul-04
 *
 */

#ifndef _MLIFPACKWRAP_
#define _MLIFPACKWRAP_

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

  /** Generates the IFPACK smoother. */
  int ML_Ifpack_Gen(ML *ml, int curr_level, int choice, int * options,
		  double * params, void ** Ifpack_Handle);
  
  /** Solves using IFPACK */
  int ML_Ifpack_Solve( void * Ifpack_Handle, double * x, double * rhs );

  /** Destroy all data associated to the IFPACK smoother. */
  void ML_Ifpack_Destroy(void * Ifpack_Handle);
  
#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
