// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#ifdef HAVE_TEUCHOS_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_TEUCHOS_MPI
#include "Teuchos_UnitTestHarness.hpp"


TEUCHOS_UNIT_TEST( Comm, Issue1029 )
{
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_BOR;
  using Teuchos::reduceAll;
  using std::endl;

  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();

  const int numProcs = comm->getSize ();  
  TEST_ASSERT( numProcs > 1 );
  if (numProcs < 2) {
    out << "This test requires at least 2 MPI processes." << endl;
    return;
  }

  const int myRank = comm->getRank ();

  // Use a mix of positive, negative, and zero values, so that neither
  // min nor max would give the full picture.
  int lclResult;
  if (myRank == 0) {
    lclResult = 0;
  }
  else if (myRank == 1) {
    lclResult = -1;
  }
  else {
    lclResult = 1;
  }
  int gblResult = 0; // output argument
  reduceAll<int, int> (*comm, REDUCE_BOR, lclResult, outArg (gblResult));
  TEST_EQUALITY( gblResult, -1 );
}

