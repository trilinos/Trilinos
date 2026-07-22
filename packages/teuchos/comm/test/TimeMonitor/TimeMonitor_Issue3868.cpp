// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#ifdef HAVE_TEUCHOS_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_TEUCHOS_MPI

TEUCHOS_UNIT_TEST(TimeMonitor, Issue3868) {
  using Teuchos::TimeMonitor;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MAX;
  using Teuchos::reduceAll;

  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();

  const int numProcs = comm->getSize ();
  TEST_ASSERT( numProcs > 1 );
  if (numProcs < 4) {
    out << "This test requires at least 4 MPI processes." << std::endl;
    return;
  }

  const int myRank = comm->getRank();

  if (myRank == 3) {
    Teuchos::TimeMonitor timer(*Teuchos::TimeMonitor::getNewTimer(std::string("myTimer")));
  }

  std::ostringstream out1;
  Teuchos::TimeMonitor::summarize(out1,false,true,false,Teuchos::Union,"",true);
  int test = (out1.str().find("myTimer") != std::string::npos);

  int gblTest = false; // output argument
  reduceAll<int, int> (*comm, REDUCE_MAX, test, outArg (gblTest));
  TEST_EQUALITY(gblTest, 1);

}
