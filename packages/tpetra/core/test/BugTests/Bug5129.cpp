// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Comm.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
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

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
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
}
