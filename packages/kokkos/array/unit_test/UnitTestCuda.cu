/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#include <iostream>

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>

#include <Kokkos_DeviceCuda.hpp>
#include <Kokkos_DeviceCuda_ValueView.hpp>
#include <Kokkos_DeviceCuda_MultiVectorView.hpp>
#include <Kokkos_DeviceCuda_MDArrayView.hpp>
#include <Kokkos_DeviceCuda_ParallelFor.hpp>
#include <Kokkos_DeviceCuda_ParallelReduce.hpp>


//----------------------------------------------------------------------------

#include <Kokkos_DeviceCuda_macros.hpp>
#include <UnitTestDeviceMemoryManagement.hpp>
#include <UnitTestValueView.hpp>
#include <UnitTestMultiVectorView.hpp>
#include <UnitTestMDArrayView.hpp>
#include <UnitTestMDArrayDeepCopy.hpp>
#include <UnitTestMDArrayIndexMap.hpp>
#include <UnitTestReduce.hpp>
#include <UnitTestMultiReduce.hpp>

namespace Test {

void test_device_cuda()
{
#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )
  try {
    Kokkos::DeviceCuda::initialize();

    UnitTestDeviceMemoryManagement< Kokkos::DeviceCuda >();
    UnitTestValueView<       Kokkos::DeviceCuda >();
    UnitTestMultiVectorView< Kokkos::DeviceCuda >();
    UnitTestMDArrayView<     Kokkos::DeviceCuda >();
    UnitTestMDArrayDeepCopy< Kokkos::DeviceCuda >();

    Test::UnitTestMDArrayIndexMap< Kokkos::DeviceCuda >();

    UnitTestReduce< long ,   Kokkos::DeviceCuda >( 1000000 );
    UnitTestReduce< double , Kokkos::DeviceCuda >( 1000000 );
    UnitTestReduceMulti< long , Kokkos::DeviceCuda >( 1000000 , 7 );

    std::cout << "PASSED : UnitTestCuda" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : UnitTestCuda : " << x.what() << std::endl ;
  }
#else
  std::cout << "PASSED : SKIPPED UnitTestCuda - NO DEVICE CUDA" << std::endl ;
#endif
}

}

#include <Kokkos_DeviceClear_macros.hpp>

//----------------------------------------------------------------------------

