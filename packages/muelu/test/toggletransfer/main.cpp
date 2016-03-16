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
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <MueLu_Exceptions.hpp>
#include <MueLu_ParameterListInterpreter.hpp>

#include "MueLu_SemiCoarsenPFactory.hpp" // for semi-coarsening constants

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
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

#include <unistd.h>
/**********************************************************************************/

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  Teuchos::oblackholestream blackhole;

  bool success = true;
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
    std::string xmlFile="parameters.xml";

    //clp.setOption("aggOrdering",&aggOrdering,"aggregation ordering strategy (natural,graph)");
    clp.setOption("debug",&pauseForDebugger,"pause to attach debugger");
    clp.setOption("maxNbrSel",&maxNbrAlreadySelected,"maximum # of nbrs allowed to be in other aggregates");
    clp.setOption("minPerAgg",&minPerAgg,"minimum #DOFs per aggregate");
    clp.setOption("timings",&printTimings,"print timings to screen");
    clp.setOption("xmlFile",&xmlFile,"file name containing MueLu multigrid parameters in XML format");

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
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFile, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *comm);

    // create parameter list interpreter
    Teuchos::RCP<HierarchyManager> mueluFactory = Teuchos::rcp(new ParameterListInterpreter(paramList));

    Teuchos::RCP<Hierarchy> H = mueluFactory->CreateHierarchy();

    H->GetLevel(0)->template Set< Teuchos::RCP<Matrix> >("A", A);

    Teuchos::RCP<MultiVector> nullspace = MultiVectorFactory::Build(A->getRowMap(), 1);
    nullspace->putScalar(1.0);
    H->GetLevel(0)->Set("Nullspace", nullspace);

    // set minimal information about number of layers for semicoarsening...
    // This information can also be provided as a user parameter in the xml file using the
    // parameter: "semicoarsen: num layers"
    // TAW: 3/16/2016: NumZLayers has to be provided as LO
    GO zLayers = matrixParameters.GetParameterList().template get<GO>("nz");
    H->GetLevel(0)->Set("NumZLayers",Teuchos::as<LO>(zLayers));


    mueluFactory->SetupHierarchy(*H);

    for (int l=0; l<H->GetNumLevels(); l++) {
      Teuchos::RCP<MueLu::Level> level = H->GetLevel(l);
      if(level->IsAvailable("A", MueLu::NoFactory::get()) == false) { success = false; H->GetLevel(l)->print(std::cout, MueLu::Debug);}
      if(level->IsAvailable("P", MueLu::NoFactory::get()) == false && l>0) { success = false; H->GetLevel(l)->print(std::cout, MueLu::Debug);}
      if(level->IsAvailable("R", MueLu::NoFactory::get()) == false && l>0) { success = false; H->GetLevel(l)->print(std::cout, MueLu::Debug);}
      if(level->IsAvailable("PreSmoother",  MueLu::NoFactory::get()) == false) { success = false; H->GetLevel(l)->print(std::cout, MueLu::Debug);}
      if(level->IsAvailable("PostSmoother", MueLu::NoFactory::get()) == false && l<H->GetNumLevels()-1) { success = false; H->GetLevel(l)->print(std::cout, MueLu::Debug);}
      if(level->IsAvailable("NumZLayers",   MueLu::NoFactory::get()) == true && l>0) {  success = false; H->GetLevel(l)->print(std::cout, MueLu::Debug);}
      H->GetLevel(l)->print(std::cout, MueLu::Debug);
    }
    ///////////////////////////////////////////////////////////

    // =========================================================================
    // System solution (Ax = b)
    // =========================================================================
    comm->barrier();
    typedef Teuchos::ScalarTraits<SC> STS;
    SC zero = STS::zero(), one = STS::one();

    Teuchos::RCP<Vector> X = VectorFactory::Build(A->getRowMap());
    Teuchos::RCP<Vector> B = VectorFactory::Build(A->getRowMap());

    {
      // we set seed for reproducibility
      Utilities::SetRandomSeed(*comm);
      X->randomize();
      A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

      Teuchos::Array<typename STS::magnitudeType> norms(1);
      B->norm2(norms);
      B->scale(one/norms[0]);
      X->putScalar(zero);
    }

    comm->barrier();

    H->IsPreconditioner(false);
    H->Iterate(*B, *X, 20);

    // Timer final summaries
    globalTimeMonitor = Teuchos::null; // stop this timer before summary

    if (printTimings)
      Teuchos::TimeMonitor::summarize();
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
