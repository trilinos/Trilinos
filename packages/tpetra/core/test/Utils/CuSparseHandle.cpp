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
// ************************************************************************
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_CuSparseHandle_fwd.hpp"
#include "Tpetra_Details_CuSparseMatrix_fwd.hpp"
#include "Tpetra_Details_CuSparseVector_fwd.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Kokkos_Core.hpp"

namespace { // (anonymous)

  void
  testCuSparseHandle(bool& success, Teuchos::FancyOStream& out)
  {
#if ! defined(KOKKOS_ENABLE_CUDA) || ! defined(HAVE_TPETRACORE_CUSPARSE)
    out << "Running this test requires enabling CUDA in Kokkos, "
      "and enabling the CUSPARSE TPL in Tpetra." << std::endl;
    TEUCHOS_ASSERT( false );
#else
    using Tpetra::Details::CuSparseHandle;

    Kokkos::Cuda execSpace1;
    std::shared_ptr<CuSparseHandle> h1 =
      Tpetra::Details::getCuSparseHandle(execSpace1);
    TEST_ASSERT( h1.get() != nullptr );

    Kokkos::Cuda execSpace2;
    std::shared_ptr<CuSparseHandle> h2 =
      Tpetra::Details::getCuSparseHandle(execSpace2);
    TEST_ASSERT( h2.get() != nullptr );
    // We've created the singleton already with the same CUDA stream,
    // so we should just get the same handle back.
    TEST_EQUALITY( h1.get(), h2.get() );

    // Create a new CUDA stream, and get a cuSPARSE handle for it.
    // This handle should differ from the singleton, because the
    // stream differs from the singleton's stream.
    cudaStream_t stream = nullptr;
    auto status = cudaStreamCreate(&stream);
    TEST_EQUALITY( status, cudaSuccess );
    if (status == cudaSuccess) {
      {
        Kokkos::Cuda execSpace3(stream);
        std::shared_ptr<CuSparseHandle> h3 =
          Tpetra::Details::getCuSparseHandle(execSpace3);
        TEST_ASSERT( h3.get() != nullptr );
        TEST_ASSERT( h3.get() != h1.get() );
      }
      (void) cudaStreamDestroy(stream);
    }
#endif
  }

  void
  testCuSparseVector(bool& success, Teuchos::FancyOStream& out)
  {
#if ! defined(KOKKOS_ENABLE_CUDA) || ! defined(HAVE_TPETRACORE_CUSPARSE)
    out << "Running this test requires enabling CUDA in Kokkos, "
      "and enabling the CUSPARSE TPL in Tpetra." << std::endl;
    TEUCHOS_ASSERT( false );
#else
    using Tpetra::Details::CuSparseHandle;
    using Tpetra::Details::CuSparseVector;
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;

    // Don't do cuSPARSE things unless a cuSPARSE handle is active.
    Kokkos::Cuda execSpace1;
    std::shared_ptr<CuSparseHandle> h1 =
      Tpetra::Details::getCuSparseHandle(execSpace1);
    TEST_ASSERT( h1.get() != nullptr );
    if (! success) {
      return;
    }

    for (int64_t numRows : {0, 11}) {
      {
        Kokkos::View<float*, Kokkos::CudaSpace> x_c_f
          (view_alloc("x_c_f", WithoutInitializing), numRows);
        auto vec =
          Tpetra::Details::getCuSparseVector(x_c_f.data(), numRows);
        TEST_ASSERT( vec.get() != nullptr );
      }
      {
        Kokkos::View<float*, Kokkos::CudaUVMSpace> x_u_f
          (view_alloc("x_u_f", WithoutInitializing), numRows);
        auto vec =
          Tpetra::Details::getCuSparseVector(x_u_f.data(), numRows);
        TEST_ASSERT( vec.get() != nullptr );
      }
      {
        Kokkos::View<double*, Kokkos::CudaSpace> x_c_d
          (view_alloc("x_c_d", WithoutInitializing), numRows);
        auto vec =
          Tpetra::Details::getCuSparseVector(x_c_d.data(), numRows);
        TEST_ASSERT( vec.get() != nullptr );
      }
      {
        Kokkos::View<double*, Kokkos::CudaUVMSpace> x_u_d
          (view_alloc("x_u_d", WithoutInitializing), numRows);
        auto vec =
          Tpetra::Details::getCuSparseVector(x_u_d.data(), numRows);
        TEST_ASSERT( vec.get() != nullptr );
      }
    }
#endif
  }

  void
  testCuSparseMatrix(bool& success, Teuchos::FancyOStream& out)
  {
#if ! defined(KOKKOS_ENABLE_CUDA) || ! defined(HAVE_TPETRACORE_CUSPARSE)
    out << "Running this test requires enabling CUDA in Kokkos, "
      "and enabling the CUSPARSE TPL in Tpetra." << std::endl;
    TEUCHOS_ASSERT( false );
#else
    using Tpetra::Details::CuSparseHandle;
    using Tpetra::Details::CuSparseMatrix;
    using Tpetra::Details::getCuSparseMatrix;
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    using LO = Tpetra::Details::DefaultTypes::local_ordinal_type;

    // Don't do cuSPARSE things unless a cuSPARSE handle is active.
    Kokkos::Cuda execSpace1;
    std::shared_ptr<CuSparseHandle> h1 =
      Tpetra::Details::getCuSparseHandle(execSpace1);
    TEST_ASSERT( h1.get() != nullptr );
    if (! success) {
      return;
    }

    using memory_space = Kokkos::CudaUVMSpace;
    for (int64_t numRows : {0, 1, 11}) {
      Kokkos::View<LO*, memory_space> ptr("ptr", numRows+1);
      for(int64_t numEntries : {0, 32}) {
        Kokkos::View<LO*, memory_space> ind
          (view_alloc("ind", WithoutInitializing), numEntries);
        Kokkos::View<float*, memory_space> val_f
          (view_alloc("val_f", WithoutInitializing), numEntries);
        Kokkos::View<double*, memory_space> val_d
          (view_alloc("val_d", WithoutInitializing), numEntries);

        for (int64_t numCols : {0, 1, 5, 13}) {
          const int64_t numEnt = (numRows == 0 || numCols == 0) ?
            int64_t(0) : numEntries;
          auto mat_f = getCuSparseMatrix(numRows, numCols, numEnt,
                                         ptr.data(), ind.data(),
                                         val_f.data());
          TEST_ASSERT( mat_f.get() != nullptr );
          auto mat_d = getCuSparseMatrix(numRows, numCols, numEnt,
                                         ptr.data(), ind.data(),
                                         val_d.data());
          TEST_ASSERT( mat_d.get() != nullptr );
        }
      }
    }
#endif
  }

  void
  runTest(bool& success,
          Teuchos::FancyOStream& out,
          std::function<void(bool&, Teuchos::FancyOStream&)> testFunction)
  {
    // Replace 'out' with cerr in verbose mode.  This lets us diagnose
    // crashes, since the Teuchos unit test framework normally
    // captures output until the test finishes.
    using Teuchos::FancyOStream;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    RCP<FancyOStream> myOutPtr;
    const bool verbose = Tpetra::Details::Behavior::verbose();
    if (verbose) {
      myOutPtr = Teuchos::getFancyOStream(rcpFromRef(std::cerr));
    }
    else {
      myOutPtr = rcpFromRef(out);
    }
    FancyOStream& myOut = *myOutPtr;
    testFunction(success, myOut);
  }

  TEUCHOS_UNIT_TEST( Utils, CuSparseHandle )
  {
    runTest(success, out, testCuSparseHandle);
  }

  TEUCHOS_UNIT_TEST( Utils, CuSparseVector )
  {
    runTest(success, out, testCuSparseVector);
  }

  TEUCHOS_UNIT_TEST( Utils, CuSparseMatrix )
  {
    runTest(success, out, testCuSparseMatrix);
  }

} // namespace (anonymous)

int
main(int argc, char* argv[])
{
  int errCode = 0;
  {
    Kokkos::ScopeGuard kokkosScope(argc, argv);
    errCode = Teuchos::UnitTestRepository::
      runUnitTestsFromMain(argc, argv);
  }
  return errCode;
}
