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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_Parameters.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"

#include "MueLu_Utilities.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_UCAggregationFactory.hpp"

#include "MueLu_Exceptions.hpp"

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraMatrixFactory.hpp>


typedef double Scalar;
typedef int    LocalOrdinal;
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
//typedef long long int    GlobalOrdinal;
typedef int    GlobalOrdinal;
#else
typedef int GlobalOrdinal;
#warning Teuchos support for long long not enabled.
#endif
typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;

#include "MueLu_UseShortNames.hpp"
#include <unistd.h>
/**********************************************************************************/

// Belos
#ifdef HAVE_MUELU_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosMueLuAdapter.hpp" // this header defines Belos::MueLuOp()
#endif

int main(int argc, char *argv[]) {
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);

#ifndef HAVE_TEUCHOS_LONG_LONG_INT
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
#endif

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);

  // Default is Laplace1D with nx = 8748.
  // It's a nice size for 1D and perfect aggregation. (6561=3^8)
  //Nice size for 1D and perfect aggregation on small numbers of processors. (8748=4*3^7)
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);             // manage parameters of xpetra

  // custom parameters
  std::string aggOrdering = "natural";
  int minPerAgg=2;
  int maxNbrAlreadySelected=0;

  clp.setOption("aggOrdering",&aggOrdering,"aggregation ordering strategy (natural,random,graph)");
  clp.setOption("minPerAgg",&minPerAgg,"minimum #DOFs per aggregate");
  clp.setOption("maxNbrSel",&maxNbrAlreadySelected,"maximum # of nbrs allowed to be in other aggregates");
  
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
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
  RCP<Operator> A = Galeri::Xpetra::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Operator vs. CrsOperator

  //  return EXIT_SUCCESS;
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  Level Finest;
  Finest.SetLevelID(0);  // must be level 0 for NullspaceFactory
  Finest.Set("A", A);
  
  Finest.SetFactoryManager( rcp( new FactoryManager() ));

  UCAggregationFactory UCAggFact;
  Finest.Request(UCAggFact);
  *out << "========================= Aggregate option summary  =========================" << std::endl;
  *out << "min DOFs per aggregate :                " << minPerAgg << std::endl;
  *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
  UCAggFact.SetMinNodesPerAggregate(minPerAgg);  //TODO should increase if run anything other than 1D
  UCAggFact.SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  std::transform(aggOrdering.begin(), aggOrdering.end(), aggOrdering.begin(), ::tolower);
  if (aggOrdering == "natural") {
       *out << "aggregate ordering :                    NATURAL" << std::endl;
       UCAggFact.SetOrdering(MueLu::AggOptions::NATURAL);
  } else if (aggOrdering == "random") {
       *out << "aggregate ordering :                    RANDOM" << std::endl;
       UCAggFact.SetOrdering(MueLu::AggOptions::RANDOM);
  } else if (aggOrdering == "graph") {
       *out << "aggregate ordering :                    GRAPH" << std::endl;
       UCAggFact.SetOrdering(MueLu::AggOptions::GRAPH);
  } else {
    std::string msg = "main: bad aggregation option """ + aggOrdering + """.";
    throw(MueLu::Exceptions::RuntimeError(msg));
  }
  UCAggFact.SetPhase3AggCreation(0.5);
  *out << "=============================================================================" << std::endl;

  UCAggFact.Build(Finest);

  return EXIT_SUCCESS;
}
