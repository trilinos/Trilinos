// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_Utilities.hpp"

/**********************************************************************************/
/* CREATE INITAL MATRIX                                                           */
/**********************************************************************************/
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraFactory.hpp>

#include "MueLu_UseDefaultTypes.hpp"
#include <unistd.h>
/**********************************************************************************/

int main(int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);

  // Default is Laplace1D with nx = 8748.
  // It's a nice size for 1D and perfect aggregation. (6561=3^8)
  // Nice size for 1D and perfect aggregation on small numbers of processors. (8748=4*3^7)
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748);  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                    // manage parameters of xpetra

  // custom parameters
  LO maxLevels = 3;
  LO its       = 10;
  clp.setOption("maxLevels", &maxLevels, "maximum number of levels allowed");
  clp.setOption("its", &its, "number of multigrid cycles");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  matrixParameters.check();
  xpetraParameters.check();
  // TODO: check custom parameters

  if (comm->getRank() == 0) {
    matrixParameters.print();
    xpetraParameters.print();
    // TODO: print custom parameters
  }

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  const RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());  // TODO: Matrix vs. CrsMatrixWrap
  RCP<Matrix> Op = Pr->BuildMatrix();
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  MueLu::AggregationOptions aggOptions;
  aggOptions.SetPrintFlag(6);
  aggOptions.SetMinNodesPerAggregate(2);
  aggOptions.SetMaxNeighAlreadySelected(0);
  aggOptions.SetOrdering("natural");
  UncoupledAggregationFactory aggFact(aggOptions);
  RCP<Graph> graph           = rcp(new Graph(Op->getCrsGraph(), "someGraphLabel"));
  double t0                  = MPI_Wtime();
  RCP<Aggregates> aggregates = aggFact.Build(*graph);
  double t1                  = MPI_Wtime() - t0;
  printf("Aggregation time only: %g\n", t1);

  return EXIT_SUCCESS;
}
