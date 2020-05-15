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
#include "Tpetra_Details_checkPointer.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Kokkos_Core.hpp"

namespace { // (anonymous)

  void
  testHostPointer (bool& success,
                   Teuchos::FancyOStream& out,
                   const void* ptr)
  {
    using Tpetra::Details::pointerAccessibleFromExecutionSpace;
    using Tpetra::Details::memorySpaceName;
    using std::endl;

    Kokkos::DefaultHostExecutionSpace hostExecSpace;
#ifdef HAVE_TPETRACORE_CUDA
    Kokkos::Cuda cudaExecSpace;
#endif // HAVE_TPETRACORE_CUDA

    const std::string name = memorySpaceName (ptr);
    out << "memorySpaceName returned \"" << name << "\"" << endl;
    TEST_EQUALITY( name, "HostSpace" );
    const bool canAccessFromHost =
      pointerAccessibleFromExecutionSpace (ptr, hostExecSpace);
    TEST_ASSERT( canAccessFromHost );

#ifdef HAVE_TPETRACORE_CUDA
    const bool canAccessFromCuda =
      pointerAccessibleFromExecutionSpace (ptr, cudaExecSpace);
    out << "pointerAccessibleFromExecutionSpace(ptr, Cuda()): "
        << (canAccessFromCuda ? "true" : "false") << endl;
    TEST_ASSERT( ! canAccessFromCuda );
#endif // HAVE_TPETRACORE_CUDA
  }

  void
  testCheckMemoryType (bool& success,
                       Teuchos::FancyOStream& out)
  {
    using Tpetra::Details::pointerAccessibleFromExecutionSpace;
    using Tpetra::Details::memorySpaceName;
    using std::endl;

    out << "Starting CheckMemoryType test" << endl;
    Teuchos::OSTab tab1 (out);

#ifdef HAVE_TPETRACORE_CUDA
    Kokkos::DefaultHostExecutionSpace hostExecSpace;
    Kokkos::Cuda cudaExecSpace;
#endif // HAVE_TPETRACORE_CUDA

    {
      out << "Test raw integer on host" << endl;
      Teuchos::OSTab tab2 (out);

      int rawInt = 666;
      testHostPointer (success, out, &rawInt);
      TEST_ASSERT( rawInt == 666 );

      out << "Result: " << (success ? "true" : "false") << endl;
    }
    {
      out << "Test host Kokkos::View" << endl;
      Teuchos::OSTab tab2 (out);

      using view_type =
        Kokkos::View<int, Kokkos::DefaultHostExecutionSpace>;
      view_type intView ("Host int");
      intView() = 418;
      testHostPointer (success, out, intView.data ());
      TEST_ASSERT( intView() == 418 );

      out << "Result: " << (success ? "true" : "false") << endl;
    }
#ifdef HAVE_TPETRACORE_CUDA
    {
      out << "Test raw cudaMalloc" << endl;
      Teuchos::OSTab tab2 (out);

      int* devPtr = nullptr;
      cudaError_t err = cudaMalloc (&devPtr, sizeof (int));

      if (err != cudaSuccess) {
        out << "Oh no!" << endl;
      }
      else {
        const bool canAccessFromHost =
          pointerAccessibleFromExecutionSpace (devPtr, hostExecSpace);
        TEST_ASSERT( ! canAccessFromHost );
        const bool canAccessFromCuda =
          pointerAccessibleFromExecutionSpace (devPtr, cudaExecSpace);
        TEST_ASSERT( canAccessFromCuda );

        const std::string name = memorySpaceName (devPtr);
        out << "memorySpaceName returned \"" << name << "\"" << endl;
        TEST_EQUALITY( name, "CudaSpace" );

        cudaFree (devPtr);
      }
      out << "Result: " << (success ? "true" : "false") << endl;
    }
    {
      out << "Test Kokkos::View of CudaSpace" << endl;
      Teuchos::OSTab tab2 (out);

      Kokkos::View<int, Kokkos::CudaSpace> cudaView ("Cuda");
      Kokkos::deep_copy (cudaView, 777);
      const bool canAccessFromHost =
        pointerAccessibleFromExecutionSpace (cudaView.data (), hostExecSpace);
      TEST_ASSERT( ! canAccessFromHost );
      const bool canAccessFromCuda =
        pointerAccessibleFromExecutionSpace (cudaView.data (), cudaExecSpace);
      TEST_ASSERT( canAccessFromCuda );

      const std::string name = memorySpaceName (cudaView.data ());
      out << "memorySpaceName returned \"" << name << "\"" << endl;
      TEST_EQUALITY( name, "CudaSpace" );

      int finalValue = 0;
      Kokkos::deep_copy (finalValue, cudaView);
      TEST_ASSERT( finalValue == 777 );

      out << "Result: " << (success ? "true" : "false") << endl;
    }
    {
      out << "Test Kokkos::View of CudaUVMSpace" << endl;
      Teuchos::OSTab tab2 (out);

      using view_type = Kokkos::View<int, Kokkos::CudaUVMSpace>;
      view_type view ("CudaUVMSpace");
      Kokkos::deep_copy (view, 31);
      const bool canAccessFromHost =
        pointerAccessibleFromExecutionSpace (view.data (), hostExecSpace);
      TEST_ASSERT( canAccessFromHost );
      const bool canAccessFromCuda =
        pointerAccessibleFromExecutionSpace (view.data (), cudaExecSpace);
      TEST_ASSERT( canAccessFromCuda );

      const std::string name = memorySpaceName (view.data ());
      out << "memorySpaceName returned \"" << name << "\"" << endl;
      TEST_EQUALITY( name, "CudaUVMSpace" );

      int finalValue = 0;
      Kokkos::deep_copy (finalValue, view);
      TEST_ASSERT( finalValue == 31 );

      out << "Result: " << (success ? "true" : "false") << endl;
    }
    {
      out << "Test Kokkos::View of CudaHostPinnedSpace" << endl;
      Teuchos::OSTab tab2 (out);

      using view_type = Kokkos::View<int, Kokkos::CudaHostPinnedSpace>;
      view_type view ("CudaHostPinned");
      Kokkos::deep_copy (view, 93);

      const bool canAccessFromHost =
        pointerAccessibleFromExecutionSpace (view.data (), hostExecSpace);
      TEST_ASSERT( canAccessFromHost );
      const bool canAccessFromCuda =
        pointerAccessibleFromExecutionSpace (view.data (), cudaExecSpace);
      TEST_ASSERT( canAccessFromCuda );

      // The CUDA API should be able to distinguish between
      // host-pinned memory and plain old host memory.
      const std::string name = memorySpaceName (view.data ());
      out << "memorySpaceName returned \"" << name << "\"" << endl;
      TEST_EQUALITY( name, "CudaHostPinnedSpace" );

      int finalValue = 0;
      Kokkos::deep_copy (finalValue, view);
      TEST_ASSERT( finalValue == 93 );

      out << "Result: " << (success ? "true" : "false") << endl;
    }
#endif // HAVE_TPETRACORE_CUDA
  }

  TEUCHOS_UNIT_TEST( Utils, CheckMemoryType )
  {
    // Replace 'out' with cerr in verbose mode.  This lets us diagnose
    // crashes, since the Teuchos unit test framework normally
    // captures output until the test finishes.
    using Teuchos::FancyOStream;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    RCP<FancyOStream> myOutPtr;
    const bool verbose = Tpetra::Details::Behavior::verbose ();
    if (verbose) {
      myOutPtr = Teuchos::getFancyOStream (rcpFromRef (std::cerr));
    }
    else {
      myOutPtr = rcpFromRef (out);
    }
    FancyOStream& myOut = *myOutPtr;
    testCheckMemoryType (success, myOut);
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
