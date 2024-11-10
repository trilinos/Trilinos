// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <unistd.h>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

#include <MueLu_Details_DefaultTypes.hpp>

#ifdef _MSC_VER
#define NOMINMAX
#include <windows.h>  // for Sleep
#endif

/*
   This driver simply generates a Tpetra matrix, prints it to screen, and exits.

   THIS EXAMPLE DOES NOT USE XPETRA.

   Use the "--help" option to get verbose help.
   */

int main(int argc, char** argv) {
  using Teuchos::RCP;

  typedef typename Tpetra::Map<>::local_ordinal_type LO;   // LocalOrdinal
  typedef typename Tpetra::Map<>::global_ordinal_type GO;  // GlobalOrdinal
  typedef MueLu::DefaultScalar SC;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    /**********************************************************************************/
    /* SET TEST PARAMETERS                                                            */
    /**********************************************************************************/
    // Note: use --help to list available options.
    Teuchos::CommandLineProcessor clp(false);

    Galeri::Xpetra::Parameters<GO> matrixParameters(clp);  // manage parameters of the test case

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    matrixParameters.check();
    std::cout << matrixParameters;

    /**********************************************************************************/
    /* CREATE INITAL MATRIX                                                           */
    /**********************************************************************************/
    RCP<const Tpetra::Map<LO, GO> > map = Teuchos::rcp(new Tpetra::Map<LO, GO>(matrixParameters.GetNumGlobalElements(), 0, comm));
    RCP<Galeri::Xpetra::Problem<Tpetra::Map<LO, GO>, Tpetra::CrsMatrix<SC, LO, GO>, Tpetra::MultiVector<SC, LO, GO> > > problem =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Tpetra::Map<LO, GO>, Tpetra::CrsMatrix<SC, LO, GO>, Tpetra::MultiVector<SC, LO, GO> >(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
    RCP<Tpetra::CrsMatrix<SC, LO, GO> > A = problem->BuildMatrix();

    /**********************************************************************************/
    /*                                                                                */
    /**********************************************************************************/

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    if (comm->getRank() == 0)
      std::cout << "\n================ MAP =====================================================\n"
                << std::endl;
    map->describe(*out, Teuchos::VERB_EXTREME);
    comm->barrier();
#ifdef _MSC_VER
    Sleep(1000);
#else
    sleep(1);
#endif
    if (comm->getRank() == 0)
      std::cout << "\n================ MATRIX ==================================================\n"
                << std::endl;
    A->describe(*out, Teuchos::VERB_EXTREME);

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}
