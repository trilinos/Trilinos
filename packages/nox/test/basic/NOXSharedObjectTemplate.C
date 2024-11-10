// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Tests the NOX's shared object template

#include "NOX.H"
#include "NOX_SharedObjectTemplate.H"
#include "NOX_MeritFunction_SumOfSquares.H"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

using namespace std;

int main(int argc, char *argv[])
{

  int status = 0;

  int myPID = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myPID);
#endif

  // Check verbosity level
  bool verbose = false;
  if (argc > 1)
    if (argv[1][0]=='-' && argv[1][1]=='v')
      verbose = true;

  // Create a print utility
  int outputInfo = NOX::Utils::Error;
  if (verbose)
    outputInfo += NOX::Utils::Debug;
  int printProc = 0;
  Teuchos::RCP<std::ostream> outputstream =
    Teuchos::rcp(&(std::cout), false);
  Teuchos::RCP<std::ostream> errorstream =
    Teuchos::rcp(&(std::cerr), false);
  Teuchos::RCP<NOX::Utils> u =
    Teuchos::rcp(new NOX::Utils(outputInfo, myPID, printProc, 6, outputstream,
                errorstream));

  // Let's share a merit function between status tests
  Teuchos::RCP<NOX::MeritFunction::Generic> mf =
    Teuchos::rcp(new NOX::MeritFunction::SumOfSquares(u));
  Teuchos::RCP<NOX::StatusTest::Generic> test1 =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND,
                        u.get()));
  Teuchos::RCP<NOX::StatusTest::Generic> test2 =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND,
                        u.get()));

  NOX::SharedObject<NOX::MeritFunction::Generic,NOX::StatusTest::Generic> sharedObject(mf);

  u->out(NOX::Utils::Debug)
    << "Testing isOwner() via non-const get(*newOwner)...";

  sharedObject.getObject(test1.get());

  if ( !(sharedObject.isOwner(test1.get())) )
    status +=1;

  if ( sharedObject.isOwner(test2.get()) )
    status += 1;

  if (status == 0)
    u->out(NOX::Utils ::Debug) << "passed!" << std::endl;
  else
    u->out(NOX::Utils ::Debug) << "failed!" << std::endl;

  u->out(NOX::Utils::Debug) << "Testing isOwner() via const get()...";

  sharedObject.getObject(test2.get());  // change ownership to test 2
  Teuchos::RCP<const NOX::MeritFunction::Generic> tmpPtr =
    sharedObject.getObject();

  if ( !(sharedObject.isOwner(test2.get())) )
    status +=1;

  if ( sharedObject.isOwner(test1.get()) )
    status += 1;

  if (status == 0)
    u->out(NOX::Utils ::Debug) << "passed!" << std::endl;
  else
    u->out(NOX::Utils ::Debug) << "failed!" << std::endl;



#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (myPID == printProc) {
    if (status == 0)
      std::cout << "\nTest passed!" << std::endl;
    else
      std::cout << "\nTest Failed!" << std::endl;
  }

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}
