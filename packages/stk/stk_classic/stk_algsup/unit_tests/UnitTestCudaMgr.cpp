/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_algsup/CudaDeviceMgr.hpp>
#include <stk_algsup/CudaMemoryMgr.hpp>

STKUNIT_UNIT_TEST( UnitTestCudaMgr, UnitTest)
{
#ifdef STK_HAVE_CUDA
  stk_classic::CudaDeviceMgr cuda_device;
  stk_classic::CudaMemoryMgr cuda_memory;

  //create a vector (length 100) of ones:
  std::vector<int> ibuf(100, 1);

  //allocate a CUDA buffer and copy ibuf into it:
  int* d_ibuf = cuda_memory.get_buffer(&ibuf[0], ibuf.size());
  cuda_memory.copy_to_buffer(&ibuf[0], ibuf.size(), d_ibuf);

  //create a vector (length 100) of zeros:
  std::vector<int> ibuf_copy(100,0);

  //copy the CUDA buffer into ibuf_copy:
  cuda_memory.copy_from_buffer(&ibuf_copy[0], ibuf_copy.size(), d_ibuf);

  cuda_memory.destroy_buffer(d_ibuf);

  //test passes if ibuf_copy is now the same as ibuf:
  STKUNIT_ASSERT( ibuf == ibuf_copy );

#else
  std::cout << "\nSTK_HAVE_CUDA not enabled."<< std::endl;
#endif
}
