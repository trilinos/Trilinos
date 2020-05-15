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
    typedef Tpetra::CrsMatrix<> crs_matrix_type;
    typedef typename crs_matrix_type::scalar_type SC;

    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
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
  //   RCP<const CrsGraph<LO,GO> > graph = rcp(new CrsGraph<LO,GO>(map,0,Tpetra::StaticProfile) );
  //   TEST_EQUALITY_CONST( graph != null, true );
  //   // All procs fail if any proc fails
  //   int globalSuccess_int = -1;
  //   reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
  //   TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  // }

  TEUCHOS_UNIT_TEST( DistObject, Bug5129_OverlyStrictMapComparison )
  {
    // test bug where map test checks that maps are the same object, instead of checking that maps are equivalent

    using std::endl;
    // mfh 06 Aug 2017: There was no obvious reason why this test
    // required LO=int and GO=int, nor did the original Bug 5129 bug
    // report depend on specific LO or GO types, so I relaxed this
    // requirement.  The test now uses the default LO and GO types.
#if 1
    typedef Tpetra::Map<>::local_ordinal_type LO;
    typedef Tpetra::Map<>::global_ordinal_type GO;
    typedef Tpetra::Vector<>::scalar_type SC;
#else
    typedef int LO;
    // this still has to be bigger than LO
    typedef int GO;
#endif // 1

    RCP<const Comm<int> > comm = getDefaultComm();
    const bool verbose = Tpetra::Details::Behavior::verbose ();
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      std::ostringstream os;
      os << "Proc " << comm->getRank () << ": Bug 5129 test: ";
      prefix = std::unique_ptr<std::string> (new std::string (os.str ()));
      os << endl;
      std::cerr << os.str ();
    }
    const global_size_t numGlobal = comm->getSize()*10;

    // create two separate, but identical, maps
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Create Maps" << endl;
      std::cerr << os.str ();
    }

    RCP<const Map<LO,GO> > mapImportIn  = Tpetra::createUniformContigMap<LO,GO>(numGlobal, comm);
    RCP<const Map<LO,GO> > mapImportOut = Tpetra::createUniformContigMap<LO,GO>(numGlobal, comm);
    RCP<const Map<LO,GO> > mapIn        = Tpetra::createUniformContigMap<LO,GO>(numGlobal, comm);
    RCP<const Map<LO,GO> > mapOut       = Tpetra::createUniformContigMap<LO,GO>(numGlobal, comm);
    TEST_EQUALITY_CONST( *mapImportIn  == *mapIn,  true );
    TEST_EQUALITY_CONST( *mapImportOut == *mapOut, true );
    TEST_EQUALITY_CONST( mapImportIn   == mapIn,  false );
    TEST_EQUALITY_CONST( mapImportOut  == mapOut, false );
    // create import, vectors from these maps

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Create Import" << endl;
      std::cerr << os.str ();
    }
    RCP<const Import<LO,GO> > import = Tpetra::createImport(mapImportIn, mapImportOut);

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Create Vectors" << endl;
      std::cerr << os.str ();
    }
    RCP<Vector<SC,LO,GO> > vecIn = Tpetra::createVector<SC>(mapIn);
    RCP<Vector<SC,LO,GO> > vecOut = Tpetra::createVector<SC>(mapOut);

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Call doImport" << endl;
      std::cerr << os.str ();
    }
    // do the import; under the bug, this should throw an exception
    TEST_NOTHROW( vecOut->doImport( *vecIn, *import, Tpetra::REPLACE ) )
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  TEUCHOS_UNIT_TEST( DistObject, BlockedViews_7234 )
  {
    // Test that replicates Trilinos issue #7234 
    // https://github.com/trilinos/Trilinos/issues/7234
    // On CUDA platforms, subviews of Kokkos::DualView did not perform properly
    // when pointing to the end of a DualView, as the shared Vector does
    // below.  See Kokkos issues #2981 and #2979 for more details
    // https://github.com/kokkos/kokkos/issues/2981
    // https://github.com/kokkos/kokkos/issues/2979
    using Teuchos::RCP;
    using Teuchos::rcp;
    using namespace Tpetra;
    using namespace Kokkos;
    using GST = Tpetra::global_size_t;
    using LO = int;
    using GO = Tpetra::Map<>::global_ordinal_type;
    using NodeType = Tpetra::Map<>::node_type;
    using TpetraMap = Tpetra::Map<LO,GO,NodeType>;
    using TpetraImport = Tpetra::Import<LO,GO,NodeType>;
    using TpetraExport = Tpetra::Export<LO,GO,NodeType>;
    using TpetraVec = Tpetra::Vector<double,LO,GO,NodeType>;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    const int rank = comm->getRank();
    const int size = comm->getSize();

    // Build the array: GIDs are BLOCKED - owned first, then shared
    Teuchos::Array<GO> owned_gids;
    Teuchos::Array<GO> shared_gids;
    Teuchos::Array<GO> owned_and_shared_gids;
    owned_gids.push_back(2*rank);
    owned_gids.push_back(2*rank+1);
    owned_and_shared_gids.push_back(2*rank);
    owned_and_shared_gids.push_back(2*rank+1);
    if (rank > 0) {
      shared_gids.push_back(owned_gids[0]-1);
      owned_and_shared_gids.push_back(shared_gids[0]);
    }

    out.setShowProcRank(true);
    out.setOutputToRootOnly(rank);
    out << "owned_map\n";
    for (int i=0; i < owned_gids.size(); ++i)
      out << i << " " << owned_gids[i] << std::endl;
    out << "owned_and_shared_map\n";
    for (int i=0; i < owned_and_shared_gids.size(); ++i)
      out << i << " " << owned_and_shared_gids[i] << std::endl;
    out.setOutputToRootOnly(0);

    RCP<TpetraMap> owned_map = rcp(new TpetraMap(Teuchos::OrdinalTraits<GST>::invalid(),owned_gids,GO(0),comm));
    RCP<TpetraMap> shared_map = rcp(new TpetraMap(Teuchos::OrdinalTraits<GST>::invalid(),shared_gids,GO(0),comm));
    RCP<TpetraMap> owned_and_shared_map = rcp(new TpetraMap(Teuchos::OrdinalTraits<GST>::invalid(),owned_and_shared_gids,GO(0),comm));

    bool zeroOut = true;
    RCP<TpetraVec> owned_and_shared = rcp(new TpetraVec(owned_and_shared_map,zeroOut));
    RCP<TpetraExport> exporter = rcp(new TpetraExport(shared_map,owned_map));
    RCP<TpetraImport> importer = rcp(new TpetraImport(owned_map,shared_map));

    // If we block the owned_and_shared, we can make the code much
    // more efficient with subviews.
    RCP<TpetraVec> owned = owned_and_shared->offsetViewNonConst(owned_map,0);
    RCP<TpetraVec> shared = owned_and_shared->offsetViewNonConst(shared_map,owned->getLocalLength());

    owned->putScalar(1.0);
    shared->doImport(*owned, *importer, Tpetra::REPLACE);

    shared->sync_host(); // sync device to host
    auto host_shared = shared->getLocalViewHost();
    const double tol = Teuchos::ScalarTraits<double>::eps() * 100.0;
    if (rank > 0) {
      TEST_FLOATING_EQUALITY(host_shared(0,0),1.0,tol);
    }

    owned->sync_host();
    auto host_owned_and_shared = owned_and_shared->getLocalViewHost();
    TEST_FLOATING_EQUALITY(host_owned_and_shared(0,0),1.0,tol);
    TEST_FLOATING_EQUALITY(host_owned_and_shared(1,0),1.0,tol);
    if (rank > 0) {
      TEST_FLOATING_EQUALITY(host_owned_and_shared(2,0),1.0,tol);
    }

    // now test the export
    owned->doExport(*shared, *exporter, Tpetra::ADD);

    owned->sync_host(); // sync device to host
    auto host_owned = owned->getLocalViewHost();
    TEST_FLOATING_EQUALITY(host_owned(0,0),1.0,tol);
    TEST_FLOATING_EQUALITY(host_owned_and_shared(0,0),1.0,tol); // check owned entries only
    if (rank != size-1) {
      TEST_FLOATING_EQUALITY(host_owned(1,0),2.0,tol);
      TEST_FLOATING_EQUALITY(host_owned_and_shared(1,0),2.0,tol); // check owned entries only
    }
    else {
      TEST_FLOATING_EQUALITY(host_owned(1,0),1.0,tol);
      TEST_FLOATING_EQUALITY(host_owned_and_shared(1,0),1.0,tol); // check owned entries only
    }
  }
}
