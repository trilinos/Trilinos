// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>

// Gallery
#define XPETRA_ENABLED  // == Gallery have to be build with the support of Xpetra matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

#include "Xpetra_ConfigDefs.hpp"

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>

#include <Xpetra_Parameters.hpp>

#include "Xpetra_UseDefaultTypes.hpp"
#include "Xpetra_UseShortNames.hpp"

/*
  This driver simply generates a Xpetra matrix, prints it to screen, and exits.

  Use the "--help" option to get verbose help.
*/

int main(int argc, char* argv[]) {
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);

  MueLu::Gallery::Parameters<GO> matrixParameters(clp);  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);              // manage parameters of xpetra

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  matrixParameters.check();
  xpetraParameters.check();

  if (comm->getRank() == 0) {
    std::cout << xpetraParameters << matrixParameters;
  }

  //   std::cout << "#threads = " << numThreads << std::endl;
  //   std::cout << "problem size = " << nx*ny << std::endl;
  //   std::cout << "matrix type = " << matrixType << std::endl;

  //   Teuchos::ParameterList pl;
  //   pl.set("Num Threads",numThreads);

  // typedef Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;

  const RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);

  {
    RCP<Matrix> A = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, Matrix>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    if (comm->getRank() == 0)
      std::cout << "\n================ MAP =====================================================\n"
                << std::endl;
    map->describe(*out, Teuchos::VERB_EXTREME);
    comm->barrier();
    //    sleep(1);

    if (comm->getRank() == 0)
      std::cout << "\n================ MATRIX ==================================================\n"
                << std::endl;
    A->describe(*out, Teuchos::VERB_EXTREME);
  }

  return EXIT_SUCCESS;
}
