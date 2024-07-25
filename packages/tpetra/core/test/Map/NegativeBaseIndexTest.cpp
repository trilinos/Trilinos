// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_OrdinalTraits.hpp"

namespace {

  using Teuchos::RCP;
  using Teuchos::Array;
  using Tpetra::global_size_t;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST( Map, Bug5401_NegativeBaseIndex )
  {
    using std::endl;
    using map_type = Tpetra::Map<>;
    using GO = Tpetra::Map<>::global_ordinal_type;
    using size_type = Teuchos::Array<GO>::size_type;

    out << "Bug 5401 (negative index base) test" << endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm ();
    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();

    out << "Comm has " << numProcs << " process"
        << (numProcs != 1 ? "es" : "") << endl;
    TEST_EQUALITY( numProcs, 2 );
    if (numProcs != 2) {
      out << "This test only works when running with exactly "
        "2 MPI processes." << endl;
      return;
    }

    TEST_ASSERT( std::is_signed<GO>::value );
    if (! std::is_signed<GO>::value) {
      out << "This test only works when the default GlobalOrdinal "
        "type is signed." << endl;
      return;
    }

    // failure reading 1x4 matrix under MPI
    const GO numElements = 78;
    const GO baseIndexIsNegOne = -1;
    const global_size_t GINV = Teuchos::OrdinalTraits<global_size_t>::invalid();
    Teuchos::Array<GO> elements (numElements);

    out << "Create array of global indices.  All processes have the same "
        << "global index.  The first global index on all processes is "
        << baseIndexIsNegOne << "." << endl;

    // first global element is -1
    for (size_type i = 0; i < elements.size (); ++i) {
      elements[i] = static_cast<GO> (i - 1);
    }

    //int localMapCtorSuccess = 0;
    RCP<map_type> map (new map_type (GINV, elements (),
                                     baseIndexIsNegOne, comm));
    out << "Process " << myRank << ":" << endl;
    {
      Teuchos::OSTab tab1 (out);
      out << "My number of global indices: " << map->getLocalNumElements ()
          << endl
          << "Global number of global indices: " << map->getGlobalNumElements ()
          << endl
          << "Index base: " << map->getIndexBase () << endl
          << "My min global index: " << map->getMinGlobalIndex () << endl
          << "Global min global index: " << map->getMinAllGlobalIndex () << endl;
    }

    TEST_EQUALITY( map->getLocalNumElements(),
                   static_cast<size_t> (numElements) );
    TEST_EQUALITY( map->getGlobalNumElements(),
                   static_cast<global_size_t> (numElements*numProcs) );
    TEST_EQUALITY( map->getIndexBase(), static_cast<GO> (-1) );
    TEST_EQUALITY( map->getMinGlobalIndex(),    static_cast<GO> (-1) );
    TEST_EQUALITY( map->getMinAllGlobalIndex(), static_cast<GO> (-1) );

    // All procs fail if any proc fails
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

}


