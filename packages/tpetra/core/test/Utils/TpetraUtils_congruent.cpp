// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Util.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#ifdef HAVE_MPI
#  include <Teuchos_DefaultMpiComm.hpp>
#endif // HAVE_MPI
#include <Teuchos_DefaultSerialComm.hpp>
#include <algorithm>
#include <iterator>

using Teuchos::Array;
using Teuchos::Comm;
//using Teuchos::as;
#ifdef HAVE_MPI
using Teuchos::MpiComm;
#endif // HAVE_MPI
using Teuchos::null;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::REDUCE_MIN;
using Teuchos::SerialComm;
using std::endl;

namespace {
  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP ();
    clp.addOutputSetupOptions (true);
  }

  //
  // UNIT TESTS
  //

#ifdef HAVE_MPI
  TEUCHOS_UNIT_TEST( TpetraUtils, Congruent2Mpi )
  {
    bool localResult = true;
    bool globalResult = true;
    int localResultInt = 1;
    int globalResultInt = 1;
    int err = MPI_SUCCESS;

    RCP<const Comm<int> > commWorld = Teuchos::DefaultComm<int>::getComm ();

    out << "Testing congruent with MPI_COMM_WORLD, MPI_COMM_WORLD" << endl;
    {
      RCP<const MpiComm<int> > comm1 = rcp (new MpiComm<int> (MPI_COMM_WORLD));
      RCP<const MpiComm<int> > comm2 = rcp (new MpiComm<int> (MPI_COMM_WORLD));

      // This should return true on all processes.  MPI_COMM_WORLD is
      // always identical (MPI_IDENT) to itself, which is included in
      // the definition of congruent.
      localResult = Tpetra::Details::congruent (*comm1, *comm2);
      localResultInt = localResult ? 1 : 0;
      reduceAll<int, int> (*commWorld, REDUCE_MIN, localResultInt,
                           outArg (globalResultInt));
      globalResult = (globalResultInt == 1);

      TEST_EQUALITY_CONST( globalResult, true );
    }

    out << "Testing congruent with MPI_COMM_SELF, MPI_COMM_SELF" << endl;
    {
      RCP<const MpiComm<int> > comm1 = rcp (new MpiComm<int> (MPI_COMM_SELF));
      RCP<const MpiComm<int> > comm2 = rcp (new MpiComm<int> (MPI_COMM_SELF));

      // This should return true on all processes.  MPI_COMM_SELF is
      // always identical (MPI_IDENT) to itself, which is included in
      // the definition of congruent.
      localResult = Tpetra::Details::congruent (*comm1, *comm2);
      localResultInt = localResult ? 1 : 0;
      reduceAll<int, int> (*commWorld, REDUCE_MIN, localResultInt,
                           outArg (globalResultInt));
      globalResult = (globalResultInt == 1);

      TEST_EQUALITY_CONST( globalResult, true );
    }

    out << "Testing congruent with MPI_COMM_WORLD, MPI_COMM_SELF" << endl;
    {
      // MPI_COMM_WORLD is only congruent to MPI_COMM_SELF if
      // MPI_COMM_WORLD has exactly one process.
      RCP<const MpiComm<int> > comm1 = rcp (new MpiComm<int> (MPI_COMM_WORLD));
      RCP<const MpiComm<int> > comm2 = rcp (new MpiComm<int> (MPI_COMM_SELF));
      const int numProcs = comm1->getSize ();

      // This should return the same result on all processes.  The
      // call to congruent() should return true if and only if
      // MPI_COMM_WORLD has one process.  We'll use REDUCE_SUM below
      // so that we can tell if any process didn't return the same
      // result as the others.
      const bool expectedGlobalResult = (numProcs == 1);
      const int expectedGlobalResultInt = expectedGlobalResult ? numProcs : 0;
      localResult = Tpetra::Details::congruent (*comm1, *comm2);
      localResultInt = localResult ? 1 : 0;
      reduceAll<int, int> (*commWorld, REDUCE_SUM, localResultInt,
                           outArg (globalResultInt));
      TEST_EQUALITY_CONST( globalResultInt, expectedGlobalResultInt );
    }

    out << "Testing congruent with MPI_COMM_WORLD, SerialComm" << endl;
    {
      // MPI_COMM_WORLD is only congruent to SerialComm if
      // MPI_COMM_WORLD has exactly one process.
      RCP<const MpiComm<int> > comm1 = rcp (new MpiComm<int> (MPI_COMM_WORLD));
      RCP<const SerialComm<int> > comm2 = rcp (new SerialComm<int>);
      const int numProcs = comm1->getSize ();

      // This should return the same result on all processes.  The
      // call to congruent() should return true if and only if
      // MPI_COMM_WORLD has one process.  We'll use REDUCE_SUM below
      // so that we can tell if any process didn't return the same
      // result as the others.
      const bool expectedGlobalResult = (numProcs == 1);
      const int expectedGlobalResultInt = expectedGlobalResult ? numProcs : 0;
      localResult = Tpetra::Details::congruent (*comm1, *comm2);
      localResultInt = localResult ? 1 : 0;
      reduceAll<int, int> (*commWorld, REDUCE_SUM, localResultInt,
                           outArg (globalResultInt));
      TEST_EQUALITY_CONST( globalResultInt, expectedGlobalResultInt );
    }

    out << "Testing congruent with MPI_COMM_WORLD, MPI_Comm_dup(MPI_COMM_WORLD)" << endl;
    {
      MPI_Comm rawComm1 = MPI_COMM_WORLD;
      MPI_Comm rawComm2;
      err = MPI_Comm_dup (rawComm1, &rawComm2);
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, "MPI_Comm_dup failed");

      RCP<const MpiComm<int> > comm1 = rcp (new MpiComm<int> (rawComm1));
      RCP<const MpiComm<int> > comm2 = rcp (new MpiComm<int> (rawComm2));

      // This should return true on all processes.  MPI_COMM_WORLD is
      // always MPI_CONGRUENT to the result of MPI_Comm_dup on itself.
      localResult = Tpetra::Details::congruent (*comm1, *comm2);
      localResultInt = localResult ? 1 : 0;
      reduceAll<int, int> (*commWorld, REDUCE_MIN, localResultInt,
                           outArg (globalResultInt));
      globalResult = (globalResultInt == 1);

      TEST_EQUALITY_CONST( globalResult, true );

      comm1 = null;
      comm2 = null;
      err = MPI_Comm_free (&rawComm2);
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, "MPI_Comm_free failed");
    }

    out << "Testing congruent with MPI_COMM_WORLD, reversed-ranks MPI_COMM_WORLD" << endl;
    {
      MPI_Comm rawComm1 = MPI_COMM_WORLD;

      MPI_Group grp1;
      err = MPI_Comm_group (rawComm1, &grp1);
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, "MPI_Comm_group failed");

      int numProcs;
      err = MPI_Group_size (grp1, &numProcs);
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, "MPI_Group_size failed");

      Array<int> ranks (numProcs); // ranks in reverse order
      for (int p = 0; p < numProcs; ++p) {
        ranks[p] = (numProcs - 1) - p;
      }

      MPI_Group grp2;
      err = MPI_Group_incl (grp1, numProcs, ranks.getRawPtr (), &grp2);
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, "MPI_Group_incl failed");

      MPI_Comm rawComm2;
      err = MPI_Comm_create (rawComm1, grp2, &rawComm2);
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, "MPI_Comm_create failed");

      RCP<const MpiComm<int> > comm1 = rcp (new MpiComm<int> (rawComm1));
      RCP<const MpiComm<int> > comm2 = rcp (new MpiComm<int> (rawComm2));

      // This should return the same result on all processes.  The
      // call to congruent() should return true if and only if
      // MPI_COMM_WORLD has one process.  We'll use REDUCE_SUM below
      // so that we can tell if any process didn't return the same
      // result as the others.
      const bool expectedGlobalResult = (numProcs == 1);
      const int expectedGlobalResultInt = expectedGlobalResult ? numProcs : 0;
      localResult = Tpetra::Details::congruent (*comm1, *comm2);
      localResultInt = localResult ? 1 : 0;
      reduceAll<int, int> (*commWorld, REDUCE_SUM, localResultInt,
                           outArg (globalResultInt));
      TEST_EQUALITY_CONST( globalResultInt, expectedGlobalResultInt );

      comm1 = null;
      comm2 = null;
      err = MPI_Comm_free (&rawComm2);
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, "MPI_Comm_free failed");

      err = MPI_Group_free (&grp2);
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, "MPI_Group_free failed");
    }
  }
#endif // HAVE_MPI

  TEUCHOS_UNIT_TEST( TpetraUtils, Congruent2Serial )
  {
    RCP<const Comm<int> > commWorld = Teuchos::DefaultComm<int>::getComm ();

    out << "Testing congruent with SerialComm, SerialComm" << endl;

    RCP<const SerialComm<int> > comm1 = rcp (new SerialComm<int>);
    RCP<const SerialComm<int> > comm2 = rcp (new SerialComm<int>);

    // This should return true on all processes.
    const bool localResult = Tpetra::Details::congruent (*comm1, *comm2);
    const int localResultInt = localResult ? 1 : 0;
    int globalResultInt = 1;
    reduceAll<int, int> (*commWorld, REDUCE_MIN, localResultInt,
                         outArg (globalResultInt));
    const bool globalResult = (globalResultInt == 1);

    TEST_EQUALITY_CONST( globalResult, true );
  }

} // namespace (anonymous)


