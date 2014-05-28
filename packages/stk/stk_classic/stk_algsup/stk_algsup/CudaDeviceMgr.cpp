/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <iostream>

#ifdef STK_HAVE_CUDA

#include <stk_algsup/CudaDeviceMgr.hpp>

namespace stk {

CudaDeviceMgr& CudaDeviceMgr::get_singleton()
{
  static CudaDeviceMgr cuda_device_mgr;
  return cuda_device_mgr;
}

CudaDeviceMgr::CudaDeviceMgr(int device)
 : m_device(device)
{
  int deviceCount = 0;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount < 1) {
    std::cout << "CudaDeviceMgr: no devices detected." << std::endl;
    //what should we do here? Abort? Throw? Continue?
  }

  if (m_device >= deviceCount) {
    std::cout << "CudaDeviceMgr: specified device not valid, using device 0." << std::endl;
    m_device = 0;
  }

  //for now: if a cuda device is already in use, just use that one. In future we may
  //want to allow for using multiple different devices...

  int deviceAlreadyBeingUsed = -1;
  cudaGetDevice( &deviceAlreadyBeingUsed );
  if (deviceAlreadyBeingUsed >= 0 && deviceAlreadyBeingUsed < deviceCount) {
    m_device = deviceAlreadyBeingUsed;
  }
  else {
    cudaSetDevice(m_device);
  }

  cudaDeviceProp deviceProp;

  cudaGetDeviceProperties(&deviceProp, m_device);

  //TODO: make this output only occur in debug mode or verbose mode, or something:
  std::cout << "\nCudaDeviceMgr attached to device #"<<m_device<<" '"
    << deviceProp.name << "', compute capability " << deviceProp.major << "." << deviceProp.minor
    << std::endl;
}

}//namespace stk

#endif

