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
#include <Teuchos_StandardCatchMacros.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Parameters.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"

#include "MueLu_Utilities.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"

#include "MueLu_Exceptions.hpp"

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

#include <unistd.h>
/**********************************************************************************/

// Belos
#ifdef HAVE_MUELU_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosMueLuAdapter.hpp"  // this header defines Belos::MueLuOp()
#endif

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;

  bool success = false;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    RCP<Teuchos::FancyOStream> out      = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);

#ifndef HAVE_XPETRA_INT_LONG_LONG
    *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
#endif

    // Default is Laplace1D with nx = 8748.
    // It's a nice size for 1D and perfect aggregation. (6561=3^8)
    // Nice size for 1D and perfect aggregation on small numbers of processors. (8748=4*3^7)
    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748);  // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);                    // manage parameters of xpetra

    // custom parameters
    std::string aggOrdering   = "natural";
    int minPerAgg             = 2;
    int maxNbrAlreadySelected = 0;

    clp.setOption("aggOrdering", &aggOrdering, "aggregation ordering strategy (natural,random,graph)");
    clp.setOption("minPerAgg", &minPerAgg, "minimum #DOFs per aggregate");
    clp.setOption("maxNbrSel", &maxNbrAlreadySelected, "maximum # of nbrs allowed to be in other aggregates");

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
      std::cout << matrixParameters << xpetraParameters << std::endl;
      // TODO: print custom parameters
    }

    /**********************************************************************************/
    /* CREATE INITIAL MATRIX                                                          */
    /**********************************************************************************/
    const RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
    Teuchos::RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());  // TODO: Matrix vs. CrsMatrixWrap
    RCP<Matrix> A = Pr->BuildMatrix();

    //  return EXIT_SUCCESS;
    /**********************************************************************************/
    /*                                                                                */
    /**********************************************************************************/

    Level Finest;
    Finest.SetLevelID(0);  // must be level 0 for NullspaceFactory
    Finest.Set("A", A);

    RCP<FactoryManager> factMngr = rcp(new FactoryManager());
    factMngr->SetKokkosRefactor(false);
    Finest.SetFactoryManager(factMngr);

    UncoupledAggregationFactory UncoupledAggFact;
    Finest.Request(UncoupledAggFact);
    *out << "========================= Aggregate option summary  =========================" << std::endl;
    *out << "min DOFs per aggregate :                " << minPerAgg << std::endl;
    *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
    UncoupledAggFact.SetMinNodesPerAggregate(minPerAgg);  // TODO should increase if run anything other than 1D
    UncoupledAggFact.SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
    std::transform(aggOrdering.begin(), aggOrdering.end(), aggOrdering.begin(), ::tolower);
    if (aggOrdering == "natural" || aggOrdering == "graph" || aggOrdering == "random") {
      *out << "aggregate ordering :                    " << aggOrdering << std::endl;
      UncoupledAggFact.SetOrdering(aggOrdering);
    } else {
      std::string msg =
          "main: bad aggregation option "
          "" +
          aggOrdering +
          ""
          ".";
      throw(MueLu::Exceptions::RuntimeError(msg));
    }
    *out << "=============================================================================" << std::endl;

    UncoupledAggFact.Build(Finest);

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
