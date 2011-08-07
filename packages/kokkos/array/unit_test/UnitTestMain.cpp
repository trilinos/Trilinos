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


#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>
#include <Kokkos_DeviceHost_ParallelFor.hpp>
#include <Kokkos_DeviceHost_ParallelReduce.hpp>

//----------------------------------------------------------------------------

#include <Kokkos_DeviceHost_macros.hpp>

#include <UnitTestDeviceMemoryManagement.hpp>
#include <UnitTestValueView.hpp>
#include <UnitTestMultiVectorView.hpp>
#include <UnitTestMDArrayView.hpp>
#include <UnitTestMDArrayDeepCopy.hpp>
#include <UnitTestMDArrayIndexMap.hpp>
#include <UnitTestReduce.hpp>
#include <UnitTestMultiReduce.hpp>

#include <Kokkos_DeviceClear_macros.hpp>

namespace Test {

void test_device_host()
{
  try {
    UnitTestDeviceMemoryManagement< Kokkos::DeviceHost >();
    UnitTestValueView<       Kokkos::DeviceHost >();
    UnitTestMultiVectorView< Kokkos::DeviceHost >();
    UnitTestMDArrayView<     Kokkos::DeviceHost >();
    UnitTestMDArrayDeepCopy< Kokkos::DeviceHost >();

    Test::UnitTestMDArrayIndexMap< Kokkos::DeviceHost >();

    UnitTestReduce< long ,   Kokkos::DeviceHost >( 1000000 );
    UnitTestReduce< double , Kokkos::DeviceHost >( 1000000 );
    UnitTestReduceMulti< long , Kokkos::DeviceHost >( 1000000 , 7 );

    std::cout << "PASSED : UnitTestHost" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : UnitTestHost : " << x.what() << std::endl ;
  }
}

}

//----------------------------------------------------------------------------

namespace Test {
void test_device_tpi();
void test_device_cuda();
void test_device_tbb();
void test_device_ferry();
}

//----------------------------------------------------------------------------

int main()
{
  Test::test_device_host();
  Test::test_device_tpi();
  Test::test_device_cuda();
  Test::test_device_tbb();
  Test::test_device_ferry();

  return 0 ;
}

