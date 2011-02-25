/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_algsup_CudaCall_hpp
#define stk_algsup_CudaCall_hpp

#include <cuda.h>
#include <cuda_runtime.h>

//----------------------------------------------------------------
inline
void stk_cuda_call(cudaError err , const char* name )
{
  if ( err != cudaSuccess ) {
    fprintf(stderr, "%s error: %s\n",name, cudaGetErrorString(err) );
    exit(-1);
  }
}

#define CUDA_CALL( cuda_fn ) stk_cuda_call( cuda_fn , #cuda_fn )


#endif

