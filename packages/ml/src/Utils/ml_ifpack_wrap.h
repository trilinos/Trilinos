/* Copyright (2004) Sandia Corportation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */ 

 /* NOTICE:  The United States Government is granted for itself and others 
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide 
 * license in ths data to reproduce, prepare derivative works, and 
 * perform publicly and display publicly.  Beginning five (5) years from 
 * July 25, 2001, the United States Government is granted for itself and 
 * others acting on its behalf a paid-up, nonexclusive, irrevocable 
 * worldwide license in this data to reproduce, prepare derivative works, 
 * distribute copies to the public, perform publicly and display 
 * publicly, and to permit others to do so.  
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT 
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES 
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR 
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY 
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS 
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */ 

#ifndef _MLIFPACKWRAP_
#define _MLIFPACKWRAP_

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

  /** Generates the IFPACK smoother. */
  int ML_Ifpack_Gen(ML *ml, int curr_level, int * options,
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
