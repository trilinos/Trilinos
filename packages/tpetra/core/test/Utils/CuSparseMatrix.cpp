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
// I'm actually using the class' methods in this test, so I need to
// include its declaration header.
#include "Tpetra_Details_CuSparseMatrix.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Kokkos_Core.hpp"

namespace { // (anonymous)

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
    using std::endl;
    using LO = Tpetra::Details::DefaultTypes::local_ordinal_type;

    out << "Test CuSparseMatrix" << endl;
    Teuchos::OSTab tab1(out);

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
      out << "numRows: " << numRows << endl;
      Teuchos::OSTab tab2(out);

      Kokkos::View<LO*, memory_space> ptr("ptr", numRows+1);
      for(int64_t numEntries : {0, 32}) {
        out << "numEntries: " << numEntries << endl;
        Teuchos::OSTab tab3(out);

        Kokkos::View<LO*, memory_space> ind
          (view_alloc("ind", WithoutInitializing), numEntries);
        Kokkos::View<float*, memory_space> val_f
          (view_alloc("val_f", WithoutInitializing), numEntries);
        Kokkos::View<double*, memory_space> val_d
          (view_alloc("val_d", WithoutInitializing), numEntries);

        for (int64_t numCols : {0, 1, 5, 13}) {
          const int64_t numEnt = (numRows == 0 || numCols == 0) ?
            int64_t(0) : numEntries;
          out << "numCols: " << numCols << ", numEnt: " << numEnt
              << endl;
          Teuchos::OSTab tab4(out);

          out << "Call getCuSparseMatrix (float)" << endl;
          auto mat_f = getCuSparseMatrix(numRows, numCols, numEnt,
                                         ptr.data(), ind.data(),
                                         val_f.data());
          out << "Test result of getCuSparseMatrix (float)" << endl;
          TEST_ASSERT( mat_f.get() != nullptr );
          if (mat_f.get() != nullptr) {
            cusparseMatDescr_t descr_f = mat_f->getDescr();
            out << "mat_f->getDescr() returned" << endl;
            TEST_ASSERT( cusparseGetMatType(descr_f) == CUSPARSE_MATRIX_TYPE_GENERAL );
            TEST_ASSERT( cusparseGetMatDiagType(descr_f) == CUSPARSE_DIAG_TYPE_NON_UNIT );
            TEST_ASSERT( cusparseGetMatIndexBase(descr_f) == CUSPARSE_INDEX_BASE_ZERO );
          }

          out << "Call getCuSparseMatrix (double)" << endl;
          auto mat_d = getCuSparseMatrix(numRows, numCols, numEnt,
                                         ptr.data(), ind.data(),
                                         val_d.data());
          out << "Test result of getCuSparseMatrix (double)" << endl;
          TEST_ASSERT( mat_d.get() != nullptr );
          if (mat_d.get() != nullptr) {
            cusparseMatDescr_t descr_d = mat_d->getDescr();
            out << "mat_d->getDescr() returned" << endl;
            TEST_ASSERT( cusparseGetMatType(descr_d) == CUSPARSE_MATRIX_TYPE_GENERAL );
            TEST_ASSERT( cusparseGetMatDiagType(descr_d) == CUSPARSE_DIAG_TYPE_NON_UNIT );
            TEST_ASSERT( cusparseGetMatIndexBase(descr_d) == CUSPARSE_INDEX_BASE_ZERO );
          }
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
