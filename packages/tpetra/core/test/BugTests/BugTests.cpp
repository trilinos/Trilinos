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

#include <Tpetra_ConfigDefs.hpp>

#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_MatrixIO.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

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

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return Tpetra::getDefaultComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
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
    typedef Tpetra::CrsMatrix<double, LO, GO> crs_matrix_type;

    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    RCP<const crs_matrix_type> readMatrix, testMatrix;

    // readHBMatrix wants a Node instance, so we need to save it.
    RCP<map_type::node_type> node;
    {
      // this is what the file looks like: 1 3 4 9
      RCP<const map_type> rng = Tpetra::createUniformContigMap<LO, GO>(1,comm);
      node = rng->getNode (); // save the Node instance for readHBMatrix
      RCP<const map_type> dom = Tpetra::createUniformContigMap<LO, GO>(4,comm);
      RCP<crs_matrix_type> A = Tpetra::createCrsMatrix<double, LO, GO>(rng);
      if (myImageID == 0) {
        A->insertGlobalValues( 0, Teuchos::tuple<GO>(0,1,2,3), Teuchos::tuple<double>(1.0,3.0,4.0,9.0) );
      }
      A->fillComplete(dom,rng);
      testMatrix = A;
    }
    {
      RCP<crs_matrix_type> A;
      Tpetra::Utils::readHBMatrix("addA2.hb", comm, node, A);
      readMatrix = A;
    }
    // test that *readMatrix == *testMatrix
    TEST_EQUALITY( testMatrix->getNodeNumRows(), readMatrix->getNodeNumRows() );
    TEST_EQUALITY( testMatrix->getNodeNumCols(), readMatrix->getNodeNumCols() );
    TEST_EQUALITY( testMatrix->getNodeNumEntries(), readMatrix->getNodeNumEntries() );
    if (success) {
      Teuchos::ArrayView<const LO>    rowinds1, rowinds2;
      Teuchos::ArrayView<const double> rowvals1, rowvals2;

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

  ////
  // TEUCHOS_UNIT_TEST( Map, Bug4756_UnsignedGlobalOrdinal )
  // {
  //   // test bug where unsigned LO clashes with signed GO
  //   typedef unsigned int LO;
  //   // this still has to be bigger than LO
  //   typedef     long int GO;
  //   // this test assumes that global_size_t (default: size_t) is larger than GO=="long int", which should be true on 64-bit builds.
  //   // create a comm
  //   RCP<const Comm<int> > comm = getDefaultComm();
  //   const global_size_t numGlobal = comm->getSize();
  //   RCP<const Map<LO,GO> > map = rcp(new Map<LO,GO>(numGlobal,0,comm) );
  //   RCP<const CrsGraph<LO,GO> > graph = rcp(new CrsGraph<LO,GO>(map,0,Tpetra::DynamicProfile) );
  //   TEST_EQUALITY_CONST( graph != null, true );
  //   // All procs fail if any proc fails
  //   int globalSuccess_int = -1;
  //   reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
  //   TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  // }

  TEUCHOS_UNIT_TEST( DistObject, Bug5129_OverlyStrictMapComparison )
  {
    // test bug where map test checks that maps are the same object, instead of checking that maps are equivalent

    // mfh 06 Aug 2017: There was no obvious reason why this test
    // required LO=int and GO=int, nor did the original Bug 5129 bug
    // report depend on specific LO or GO types, so I relaxed this
    // requirement.  The test now uses the default LO and GO types.
#if 1
    typedef Tpetra::Map<>::local_ordinal_type LO;
    typedef Tpetra::Map<>::global_ordinal_type GO;
#else
    typedef int LO;
    // this still has to be bigger than LO
    typedef int GO;
#endif // 1

    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const global_size_t numGlobal = comm->getSize()*10;
    // create two separate, but identical, maps
    RCP<const Map<LO,GO> > mapImportIn  = Tpetra::createUniformContigMap<LO,GO>(numGlobal, comm);
    RCP<const Map<LO,GO> > mapImportOut = Tpetra::createUniformContigMap<LO,GO>(numGlobal, comm);
    RCP<const Map<LO,GO> > mapIn        = Tpetra::createUniformContigMap<LO,GO>(numGlobal, comm);
    RCP<const Map<LO,GO> > mapOut       = Tpetra::createUniformContigMap<LO,GO>(numGlobal, comm);
    TEST_EQUALITY_CONST( *mapImportIn  == *mapIn,  true );
    TEST_EQUALITY_CONST( *mapImportOut == *mapOut, true );
    TEST_EQUALITY_CONST( mapImportIn   == mapIn,  false );
    TEST_EQUALITY_CONST( mapImportOut  == mapOut, false );
    // create import, vectors from these maps
    RCP<const Import<LO,GO> > import = Tpetra::createImport(mapImportIn, mapImportOut);
    RCP<Vector<double,LO,GO> > vecIn = Tpetra::createVector<double>(mapIn);
    RCP<Vector<double,LO,GO> > vecOut = Tpetra::createVector<double>(mapOut);
    // do the import; under the bug, this should throw an exception
    TEST_NOTHROW( vecOut->doImport( *vecIn, *import, Tpetra::REPLACE ) )
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }
}


