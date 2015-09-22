/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_algsup/CudaMemoryMgr.hpp>

namespace stk_classic {

#ifdef STK_HAVE_CUDA

CudaMemoryMgr& get_singleton()
{
  static CudaMemoryMgr cuda_memory_mgr;
  return cuda_memory_mgr;
}

#endif

CudaMemoryMgr::~CudaMemoryMgr()
{
#ifdef STK_HAVE_CUDA
  std::map<const void*,const void*>::iterator
    iter = device_to_host_map.begin(),
    iter_end = device_to_host_map.end();

  for(; iter!=iter_end; ++iter) {
    //cast away const so we can free the pointer:
    void* dev_ptr = const_cast<void*>(iter->first);
    CUDA_CALL( cudaFree(dev_ptr) );
  }
#endif
}

}//namespace stk_classic

