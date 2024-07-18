// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <fstream>
#include <string>
#include <sys/types.h>
#include <unistd.h>

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"

#include "Kokkos_Core.hpp"

#include "Test_Utils.hpp"
#include "src/tLU2x2PreconditionerFactory_tpetra.hpp"
#include "src/tLSCStablePreconditionerFactory_tpetra.hpp"
#include "src/tSIMPLEPreconditionerFactory_tpetra.hpp"
#include "src/tLSCStabilized_tpetra.hpp"
#include "src/tJacobi2x2PreconditionerFactory_tpetra.hpp"
#include "src/tBlockUpperTriInverseOp_tpetra.hpp"
#include "src/tBlockLowerTriInverseOp_tpetra.hpp"
#include "src/tLSCIntegrationTest_tpetra.hpp"
#include "src/tLSCHIntegrationTest_tpetra.hpp"
#include "src/tGraphLaplacian_tpetra.hpp"
#include "src/tParallelInverse_tpetra.hpp"

#include "src/tLumping_tpetra.hpp"
#include "src/tAbsRowSum_tpetra.hpp"
#include "src/tNeumannSeries_tpetra.hpp"
#include "src/tPCDStrategy_tpetra.hpp"

#ifdef TEKO_HAVE_EPETRA
#include "src/tDiagonalPreconditionerFactory_tpetra.hpp"
#include "src/tBlockJacobiPreconditionerFactory_tpetra.hpp"
#include "src/tExplicitOps_tpetra.hpp"
#include "src/Tpetra/tTpetraOperatorWrapper.hpp"
#include "src/Tpetra/tStridedTpetraOperator.hpp"
#include "src/Tpetra/tBlockedTpetraOperator.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#endif

#include "src/Tpetra/tInterlacedTpetra.hpp"
#include "src/Tpetra/tBlockingTpetra.hpp"

#include "src/Tpetra/tTpetraThyraConverter.hpp"

#include "Tpetra_Core.hpp"

void gdbIn() {
  int i = 0;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  printf("PID %d on %s ready for attach\n", getpid(), hostname);
  fflush(stdout);
  while (0 == i) sleep(5);
}

int main(int argc, char* argv[]) {
  bool status = false;

  {
    // need to protect kokkos and MPI
    // calls MPI_Init and MPI_Finalize
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    Kokkos::initialize(argc, argv);

#ifdef TEKO_HAVE_EPETRA
    // build MPI/Serial communicators
#ifdef HAVE_MPI
    Epetra_MpiComm Comm_epetra(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm_epetra;
#endif  // HAVE_MPI
    Teko::Test::UnitTest::SetComm(Teuchos::rcpFromRef(Comm_epetra));
#endif  // TEKO_HAVE_EPETRA

    Teuchos::RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();
    Teko::Test::UnitTest::SetComm_tpetra(Comm);

    Teuchos::CommandLineProcessor clp;

    int verbosity       = 1;
    std::string faillog = "failure.log";
    bool isfast         = false;

    clp.setOption("verb", &verbosity, "How verbose is the output? 1 is normal 10 is a lot.");
    clp.setOption("log", &faillog,
                  "File for failure information to go to (also high verbosity text)");
    clp.setOption("fast", "notfast", &isfast, "Run only fast tests");
    clp.parse(argc, argv);

    Teuchos::RCP<Teuchos::FancyOStream> termout =
        Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::RCP<Teuchos::FancyOStream> failout;
    std::ofstream failure;

    if (faillog == "stdout") {
      failout = termout;
    } else {
      failure.open(faillog.c_str());
      failout = Teuchos::getFancyOStream(Teuchos::rcpFromRef(failure));
    }

    termout->setOutputToRootOnly(0);
    failout->setOutputToRootOnly(0);

    // gdbIn();
    Teko_ADD_UNIT_TEST(Teko::Test::tSIMPLEPreconditionerFactory_tpetra,
                       SIMPLEPreconditionerFactory_tpetra);
#ifdef TEKO_HAVE_EPETRA
    Teko_ADD_UNIT_TEST(Teko::Test::tDiagonalPreconditionerFactory_tpetra,
                       DiagonalPreconditionerFactory_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tBlockJacobiPreconditionerFactory_tpetra,
                       BlockJacobiPreconditionerFactory_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tExplicitOps_tpetra, tExplicitOps_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tTpetraOperatorWrapper, tTpetraOperatorWrapper);
#endif
    Teko_ADD_UNIT_TEST(Teko::Test::tLU2x2PreconditionerFactory_tpetra,
                       LU2x2PreconditionerFactory_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tLSCStablePreconditionerFactory_tpetra,
                       LSCStablePreconditionerFactory_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tLSCStabilized_tpetra, LSCStabilized_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tJacobi2x2PreconditionerFactory_tpetra,
                       Jacobi2x2PreconditionerFactory_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tBlockUpperTriInverseOp_tpetra, BlockUpperTriInverseOp_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tBlockLowerTriInverseOp_tpetra, BlockLowerTriInverseOp_tpetra);

    Teko_ADD_UNIT_TEST(Teko::Test::tInterlacedTpetra, InterlacedTpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tBlockingTpetra, BlockingTpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tTpetraThyraConverter, TpetraThyraConverter);
    Teko_ADD_UNIT_TEST(Teko::Test::tGraphLaplacian_tpetra, tGraphLaplacian_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tParallelInverse_tpetra, tParallelInverse_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tLSCHIntegrationTest_tpetra, LSCHIntegrationTest_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tLumping_tpetra, Lumping_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tAbsRowSum_tpetra, AbsRowSum_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tNeumannSeries_tpetra, NeumannSeries_tpetra);
    Teko_ADD_UNIT_TEST(Teko::Test::tPCDStrategy_tpetra, PCDStrategy_tpetra);
    if (not isfast) {
      Teko_ADD_UNIT_TEST(Teko::Test::tLSCIntegrationTest_tpetra, LSCIntegrationTest_tpetra);
#ifdef TEKO_HAVE_EPETRA
      Teko_ADD_UNIT_TEST(Teko::Test::tStridedTpetraOperator, tStridedTpetraOperator);
      Teko_ADD_UNIT_TEST(Teko::Test::tBlockedTpetraOperator, tBlockedTpetraOperator);
#endif
    }

    status = Teko::Test::UnitTest::RunTests_tpetra(verbosity, *termout, *failout);

    if (not status) *termout << "Teko tests failed" << std::endl;

    // release any stored Kokkos memory
    Teko::Test::UnitTest::ClearTests();
  }

  if (!Kokkos::is_finalized()) Kokkos::finalize();

  return status ? 0 : -1;
}
