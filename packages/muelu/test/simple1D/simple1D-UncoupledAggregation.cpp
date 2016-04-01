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
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include "MueLu_Hierarchy.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_Utilities.hpp"

/**********************************************************************************/
/* CREATE INITAL MATRIX                                                           */
/**********************************************************************************/
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Galeri
// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <unistd.h>
/**********************************************************************************/

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  Teuchos::oblackholestream blackhole;

  bool success = false;
  bool verbose = true;
  try {
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    /**********************************************************************************/
    /* SET TEST PARAMETERS                                                            */
    /**********************************************************************************/

    // Default is Laplace1D with nx = 8748.
    // It's a nice size for 1D and perfect aggregation. (6561=3^8)
    //Nice size for 1D and perfect aggregation on small numbers of processors. (8748=4*3^7)
    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);             // manage parameters of xpetra

    // custom parameters
    int pauseForDebugger=0;
    //std::string aggOrdering = "natural";
    int minPerAgg=2; //was 3 in simple
    int maxNbrAlreadySelected=0;
    int printTimings=0;

    //clp.setOption("aggOrdering",&aggOrdering,"aggregation ordering strategy (natural,graph)");
    clp.setOption("debug",&pauseForDebugger,"pause to attach debugger");
    clp.setOption("maxNbrSel",&maxNbrAlreadySelected,"maximum # of nbrs allowed to be in other aggregates");
    clp.setOption("minPerAgg",&minPerAgg,"minimum #DOFs per aggregate");
    clp.setOption("timings",&printTimings,"print timings to screen");

    switch (clp.parse(argc,argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    Teuchos::RCP<Teuchos::TimeMonitor> globalTimeMonitor = Teuchos::rcp (new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Timings: Global Time")));

    if (pauseForDebugger) {
      Utilities::PauseForDebugger();
    }

    matrixParameters.check();
    xpetraParameters.check();

    if (comm->getRank() == 0) {
      std::cout << xpetraParameters << matrixParameters;
    }

    /**********************************************************************************/
    /* CREATE INITIAL MATRIX                                                          */
    /**********************************************************************************/
    Teuchos::RCP<const Map> map;
    Teuchos::RCP<Matrix> A;

    {
      Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("Timings: Matrix Build"));

      map = MapFactory::Build(lib, matrixParameters.GetNumGlobalElements(), 0, comm);
      Teuchos::RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Matrix vs. CrsMatrixWrap
      A = Pr->BuildMatrix();

    }
    /**********************************************************************************/
    /*                                                                                */
    /**********************************************************************************/


    Teuchos::RCP<Hierarchy> hierarchy = Teuchos::rcp(new Hierarchy(A));
    hierarchy->SetDefaultVerbLevel(MueLu::toMueLuVerbLevel(Teuchos::VERB_EXTREME));
    hierarchy->SetMaxCoarseSize(100);

    ///////////////////////////////////////////////////////////

    Teuchos::RCP<MueLu::Level> Finest = hierarchy->GetLevel();  // get finest level

    Finest->Set("A",A);

    Teuchos::RCP<AmalgamationFactory> amalgFact = Teuchos::rcp(new AmalgamationFactory());
    Teuchos::RCP<CoalesceDropFactory> dropFact = Teuchos::rcp(new CoalesceDropFactory());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    // aggregation factory
    Teuchos::RCP<UncoupledAggregationFactory> UnCoupledAggFact = Teuchos::rcp(new UncoupledAggregationFactory(/*dropFact*/));
    UnCoupledAggFact->SetFactory("Graph", dropFact);
    UnCoupledAggFact->SetFactory("DofsPerNode", dropFact);
    //UnCoupledAggFact->SetFactory("Graph", dropFact); // UnCoupledAggFact not changed to new factory handling
    //UnCoupledAggFact->SetMinNodesPerAggregate(minPerAgg);
    //UnCoupledAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
    //UnCoupledAggFact->SetOrdering(MueLu::AggOptions::GRAPH);
    UnCoupledAggFact->SetParameter("aggregation: max selected neighbors", Teuchos::ParameterEntry(maxNbrAlreadySelected));
    UnCoupledAggFact->SetParameter("aggregation: min agg size",           Teuchos::ParameterEntry(minPerAgg));
    UnCoupledAggFact->SetParameter("aggregation: ordering",               Teuchos::ParameterEntry(std::string("graph")));

    Finest->Request("Graph",dropFact.get());
    Finest->Request("DofsPerNode",dropFact.get());
    Finest->Request("UnAmalgamationInfo",amalgFact.get());
    UnCoupledAggFact->Build(*Finest);

    // Timer final summaries
    globalTimeMonitor = Teuchos::null; // stop this timer before summary

    if (printTimings)
      Teuchos::TimeMonitor::summarize();

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

int main(int argc, char* argv[]) {
  bool success = false;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  try {
    const bool throwExceptions     = false;
    const bool recogniseAllOptions = false;

    Teuchos::CommandLineProcessor clp(throwExceptions, recogniseAllOptions);
    Xpetra::Parameters xpetraParameters(clp);

    std::string node = "";  clp.setOption("node", &node, "node type (serial | openmp | cuda)");

    switch (clp.parse(argc, argv, NULL)) {
      case Teuchos::CommandLineProcessor::PARSE_ERROR:               return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      return main_<double,int,int,Xpetra::EpetraNode>(clp, lib, argc, argv);
#else
      throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
    }

    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      if (node == "") {
        typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#else
#    if   defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double,int,int,Node> (clp, lib, argc, argv);
#  elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#  elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double,int,long long,Node>(clp, lib, argc, argv);
#  else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#  endif
#endif
      } else if (node == "serial") {
#ifdef KOKKOS_HAVE_SERIAL
        typedef Kokkos::Compat::KokkosSerialWrapperNode Node;

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#  else
#    if   defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double,int,int,Node> (clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double,int,long long,Node>(clp, lib, argc, argv);
#    else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#    endif
#  endif
#else
        throw MueLu::Exceptions::RuntimeError("Serial node type is disabled");
#endif
      } else if (node == "openmp") {
#ifdef KOKKOS_HAVE_OPENMP
        typedef Kokkos::Compat::KokkosOpenMPWrapperNode Node;

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, argc, argv);
#  else
#    if   defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double,int,int,Node> (clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double,int,long long,Node>(clp, lib, argc, argv);
#    else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#    endif
#  endif
#else
        throw MueLu::Exceptions::RuntimeError("OpenMP node type is disabled");
#endif
      } else if (node == "cuda") {
#ifdef KOKKOS_HAVE_CUDA
        typedef Kokkos::Compat::KokkosCudaWrapperNode Node;

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, argc, argv);
#  else
#    if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double,int,int,Node> (clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double,int,long long,Node>(clp, lib, argc, argv);
#    else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#    endif
#  endif
#else
        throw MueLu::Exceptions::RuntimeError("CUDA node type is disabled");
#endif
      } else {
        throw MueLu::Exceptions::RuntimeError("Unrecognized node type");
      }
#else
      throw MueLu::Exceptions::RuntimeError("Tpetra is not available");
#endif
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

