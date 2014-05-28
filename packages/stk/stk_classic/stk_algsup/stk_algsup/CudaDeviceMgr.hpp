/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_algsup_CudaDeviceMgr_hpp
#define stk_algsup_CudaDeviceMgr_hpp

#ifdef STK_HAVE_CUDA

#include <stk_algsup/CudaCall.hpp>

namespace stk {

class CudaDeviceMgr {
  public:
    CudaDeviceMgr(int device=0);

    virtual ~CudaDeviceMgr() {}

    int get_device() const { return m_device; }

    static CudaDeviceMgr& get_singleton();

 private:
  int m_device;
};//class CudaMemoryMgr

}//namespace stk

#endif

#endif

