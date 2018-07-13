// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER


#include <Teuchos_UnitTestRepository.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#ifdef HAVE_MUELU_KOKKOSCORE
#include <Kokkos_Core.hpp>
#endif

#include "MueLu_TestHelpers_kokkos.hpp"

/*#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_UnitTestRepository.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Kokkos_Core.hpp>

#include "MueLu_TestHelpers_kokkos.hpp"

#include "MueLu_VerboseObject.hpp"*/

int main(int argc, char* argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Kokkos::initialize(argc, argv);

  bool success = false;
  bool verbose = true;
  int ierr = -1;
  try {
    // Note: the command line parameter --linAlgebra= is taken into account.
    // Xpetra parameters are added to the Teuchos::CommandLineProcessor of Teuchos::UnitTestRepository in MueLu_TestHelpers_kokkos.cpp

#ifdef ParallelDebug
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    int mypid = comm->getRank();

    if (mypid  == 0) std::cout << "Host and Process Ids for tasks" << std::endl;
    for (int i = 0; i <comm->getSize(); i++) {
      if (i == mypid ) {
        char buf[80];
        char hostname[80];
        gethostname(hostname, sizeof(hostname));
        int pid = getpid();
        sprintf(buf, "Host: %s\tMPI rank: %d,\tPID: %d\n\tattach %d\n\tcontinue\n",
                hostname, mypid, pid, pid);
        printf("%s\n",buf);
        fflush(stdout);
        sleep(1);
      }
    }

    if (mypid == 0) {
      printf( "** Enter a character to continue > "); fflush(stdout);
      char go = ' ';
      scanf("%c",&go);
    }
    comm->barrier();
#endif

    // Comment this line to get rid of MueLu output
    MueLu::VerboseObject::SetDefaultVerbLevel(MueLu::High);

    ierr = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
  Kokkos::finalize();
  return (success ? ierr : EXIT_FAILURE);
}
