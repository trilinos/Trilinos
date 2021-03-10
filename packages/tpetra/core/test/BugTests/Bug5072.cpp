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

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Comm.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include <memory>
#include <sstream>

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;
  using Teuchos::null;
  using Tpetra::Import;
  using Tpetra::Map;
  using Tpetra::CrsGraph;
  using Tpetra::CrsMatrix;
  using Tpetra::Vector;
  using Tpetra::global_size_t;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }


  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST( readHBMatrix, Bug5072_ReadOneRowMPI )
  {
    // failure reading 1x4 matrix under MPI
    // typedef int                       LO;
    // typedef int                       GO;
    typedef Tpetra::Map<> map_type;
    typedef map_type::local_ordinal_type LO;
    typedef map_type::global_ordinal_type GO;
    typedef Tpetra::CrsMatrix<> crs_matrix_type;
    typedef typename crs_matrix_type::scalar_type SC;

    // create a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int myImageID = comm->getRank();
    RCP<const crs_matrix_type> readMatrix, testMatrix;

    {
      // this is what the file looks like: 1 3 4 9
      RCP<const map_type> rng = Tpetra::createUniformContigMap<LO, GO>(1,comm);
      RCP<const map_type> dom = Tpetra::createUniformContigMap<LO, GO>(4,comm);
      RCP<crs_matrix_type> A = Tpetra::createCrsMatrix<SC, LO, GO>(rng, 4);
      if (myImageID == 0) {
        A->insertGlobalValues( 0, Teuchos::tuple<GO>(0,1,2,3), Teuchos::tuple<SC>(1.0,3.0,4.0,9.0) );
      }
      A->fillComplete(dom,rng);
      testMatrix = A;
    }
    {
      RCP<crs_matrix_type> A;
      Tpetra::Utils::readHBMatrix("addA2.hb", comm, A);
      readMatrix = A;
    }
    // test that *readMatrix == *testMatrix
    TEST_EQUALITY( testMatrix->getNodeNumRows(), readMatrix->getNodeNumRows() );
    TEST_EQUALITY( testMatrix->getNodeNumCols(), readMatrix->getNodeNumCols() );
    TEST_EQUALITY( testMatrix->getNodeNumEntries(), readMatrix->getNodeNumEntries() );
    if (success) {
      Teuchos::ArrayView<const LO>    rowinds1, rowinds2;
      Teuchos::ArrayView<const SC> rowvals1, rowvals2;

      const LO lclNumRows = testMatrix->getNodeNumRows ();
      for (LO r = 0; r < lclNumRows; ++r) {
        testMatrix->getLocalRowView(r, rowinds1, rowvals1);
        readMatrix->getLocalRowView(r, rowinds2, rowvals2);
        TEST_COMPARE_ARRAYS( rowinds1, rowinds2 );
        TEST_COMPARE_ARRAYS( rowvals1, rowvals2 );
      }
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }
}
