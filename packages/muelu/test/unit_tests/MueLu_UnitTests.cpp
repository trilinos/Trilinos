// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file main.cpp

  \brief MueLu unit testing main program.

  This file is the main for the unit test executable.

NOTE: This file should *not* be built and included as part of the MueLu
library.  It is instead to be directly included in the build files for
specific unit test suites.

*/

#include <Teuchos_UnitTestRepository.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Kokkos_Core.hpp>

#include "MueLu_TestHelpers.hpp"

int main(int argc, char* argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Kokkos::initialize(argc, argv);

  bool success = false;
  bool verbose = true;
  int ierr     = -1;
  try {
    // Note: the command line parameter --linAlgebra= is take into account.
    // Xpetra parameters are added to the Teuchos::CommandLineProcessor of Teuchos::UnitTestRepository in MueLu_TestHelpers.cpp

#ifdef ParallelDebug
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    int mypid = comm->getRank();

    if (mypid == 0) std::cout << "Host and Process Ids for tasks" << std::endl;
    for (int i = 0; i < comm->getSize(); i++) {
      if (i == mypid) {
        char buf[80];
        char hostname[80];
        gethostname(hostname, sizeof(hostname));
        int pid = getpid();
        sprintf(buf, "Host: %s\tMPI rank: %d,\tPID: %d\n\tattach %d\n\tcontinue\n",
                hostname, mypid, pid, pid);
        printf("%s\n", buf);
        fflush(stdout);
        sleep(1);
      }
    }

    if (mypid == 0) {
      printf("** Enter a character to continue > ");
      fflush(stdout);
      char go = ' ';
      scanf("%c", &go);
    }
    comm->barrier();
#endif

    ierr = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  Kokkos::finalize();

  return (success ? ierr : EXIT_FAILURE);
}
