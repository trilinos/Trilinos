/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_checkMemoryType.hpp"
#include "Kokkos_Core.hpp"

namespace { // (anonymous)

  void
  testHostPointer (bool& success,
                   Teuchos::FancyOStream& out,
                   const void* ptr)
  {
    using Tpetra::Details::inMemorySpace;
    using Tpetra::Details::memorySpaceName;

    Kokkos::HostSpace hostSpace;
#ifdef HAVE_TPETRACORE_CUDA
    Kokkos::CudaSpace cudaSpace;
    Kokkos::CudaSpace cudaUVMSpace;
    Kokkos::CudaHostPinnedSpace cudaHostPinnedSpace;
#endif // HAVE_TPETRACORE_CUDA

    const std::string name = memorySpaceName (ptr);
    TEST_EQUALITY( name, "HostSpace" );
    const bool inHostSpace = inMemorySpace (ptr, hostSpace);
    TEST_ASSERT( inHostSpace );

#ifdef HAVE_TPETRACORE_CUDA
    const bool inCudaSpace = inMemorySpace (ptr, cudaSpace);
    TEST_ASSERT( ! inCudaSpace );
    const bool inCudaUVMSpace = inMemorySpace (ptr, cudaUVMSpace);
    TEST_ASSERT( ! inCudaUVMSpace );
    const bool inCudaHostPinnedSpace = inMemorySpace (ptr, cudaHostPinnedSpace);
    TEST_ASSERT( ! inCudaHostPinnedSpace );
#endif // HAVE_TPETRACORE_CUDA
  }

  TEUCHOS_UNIT_TEST( Utils, CheckMemoryType )
  {
    using Tpetra::Details::inMemorySpace;
    using Tpetra::Details::memorySpaceName;

    Kokkos::HostSpace hostSpace;
#ifdef HAVE_TPETRACORE_CUDA
    Kokkos::CudaSpace cudaSpace;
    Kokkos::CudaSpace cudaUVMSpace;
    Kokkos::CudaHostPinnedSpace cudaHostPinnedSpace;
#endif // HAVE_TPETRACORE_CUDA

    {
      int rawInt = 666;
      testHostPointer (success, out, &rawInt);
      TEST_ASSERT( rawInt == 666 );
    }
    {
      Kokkos::View<int, Kokkos::DefaultHostExecutionSpace> intView;
      intView() = 418;
      testHostPointer (success, out, intView.data ());
      TEST_ASSERT( intView() == 418 );
    }
#ifdef HAVE_TPETRACORE_CUDA
    {
      Kokkos::View<int, Kokkos::CudaSpace> cudaView;
      Kokkos::deep_copy (cudaView, 777);
      const bool inHostSpace =
        inMemorySpace (cudaView.data (), hostSpace);
      TEST_ASSERT( ! inHostSpace );
      const bool inCudaSpace =
        inMemorySpace (cudaView.data (), cudaSpace);
      TEST_ASSERT( inCudaSpace );
      const bool inCudaUVMSpace =
        inMemorySpace (cudaView.data (), cudaUVMSpace);
      TEST_ASSERT( ! inCudaUVMSpace );
      const bool inCudaHostPinnedSpace =
        inMemorySpace (cudaView.data (), cudaHostPinnedSpace);
      TEST_ASSERT( ! inCudaHostPinnedSpace );

      const std::string name = memorySpaceName (cudaView.data ());
      TEST_EQUALITY( name, "CudaSpace" );

      int finalValue = 0;
      Kokkos::deep_copy (finalValue, cudaView);
      TEST_ASSERT( finalValue == 777 );
    }
    {
      Kokkos::View<int, Kokkos::CudaUVMSpace> cudaUVMView;
      Kokkos::deep_copy (cudaUVMView, 31);
      const bool inHostSpace =
        inMemorySpace (cudaView.data (), hostSpace);
      // It's host accessible, but not in HostSpace.
      TEST_ASSERT( ! inHostSpace );
      const bool inCudaSpace =
        inMemorySpace (cudaUVMView.data (), cudaSpace);
      TEST_ASSERT( ! inCudaSpace );
      const bool inCudaUVMSpace =
        inMemorySpace (cudaUVMView.data (), cudaUVMSpace);
      TEST_ASSERT( inCudaUVMSpace );
      const bool inCudaHostPinnedSpace =
        inMemorySpace (cudaUVMView.data (), cudaHostPinnedSpace);
      TEST_ASSERT( ! inCudaHostPinnedSpace );

      const std::string name = memorySpaceName (cudaUVMView.data ());
      TEST_EQUALITY( name, "CudaUVMSpace" );

      int finalValue = 0;
      Kokkos::deep_copy (finalValue, cudaUVMView);
      TEST_ASSERT( finalValue == 31 );
    }
    {
      Kokkos::View<int, Kokkos::CudaHostPinnedSpace> cudaHostPinnedView;
      Kokkos::deep_copy (cudaHostPinnedView, 31);
      // The CUDA API should be able to distinguish between
      // host-pinned memory and plain old host memory.
      const bool inHostSpace =
        inMemorySpace (cudaView.data (), hostSpace);
      TEST_ASSERT( ! inHostSpace );
      const bool inCudaSpace =
        inMemorySpace (cudaHostPinnedView.data (), cudaSpace);
      TEST_ASSERT( ! inCudaSpace );
      const bool inCudaUVMSpace =
        inMemorySpace (cudaHostPinnedView.data (), cudaUVMSpace);
      TEST_ASSERT( ! inCudaUVMSpace );
      const bool inCudaHostPinnedSpace =
        inMemorySpace (cudaHostPinnedView.data (), cudaHostPinnedSpace);
      TEST_ASSERT( inCudaHostPinnedSpace );

      const std::string name = memorySpaceName (cudaHostPinnedView.data ());
      TEST_EQUALITY( name, "CudaHostPinnedSpace" );

      int finalValue = 0;
      Kokkos::deep_copy (finalValue, cudaHostPinnedView);
      TEST_ASSERT( finalValue == 31 );
    }
#endif // HAVE_TPETRACORE_CUDA
  }

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  // Initialize MPI (if enabled) before initializing Kokkos.  This
  // lets MPI control things like pinning processes to sockets.
  Teuchos::GlobalMPISession mpiSession (&argc, &argv);
  Kokkos::initialize (argc, argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  Kokkos::finalize ();
  return errCode;
}
