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
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <unistd.h>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_YamlParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Xpetra
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_IO.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

// MueLu
#include <MueLu.hpp>

#include <MueLu_BaseClass.hpp>
#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif
#include <MueLu_Level.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>

// Belos
#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosBiCGStabSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#ifdef HAVE_MUELU_TPETRA
#include <BelosTpetraAdapter.hpp>    // => This header defines Belos::TpetraOp
#endif
#endif


#ifdef HAVE_MUELU_CUDA
#include "cuda_profiler_api.h"
#endif

// MueLu and Xpetra Tpetra stack
#ifdef HAVE_MUELU_TPETRA
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include "Xpetra_TpetraMultiVector.hpp"
#include <KokkosBlas1_abs.hpp>
#include <Tpetra_leftAndOrRightScaleCrsMatrix.hpp>
#include <Tpetra_computeRowAndColumnOneNorms.hpp>
#endif

// Xpetra Epetra stack
#ifdef HAVE_MUELU_EPETRA
#include "Xpetra_EpetraMultiVector.hpp"
#endif

#include <MueLu_CreateXpetraPreconditioner.hpp>

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)
#include <Amesos2_config.h>
#include <Amesos2.hpp>
#endif

// Region MG headers
#include "SetupRegionHierarchy_def.hpp"


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib& lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;
  using Teuchos::TimeMonitor;
  using Teuchos::ParameterList;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  const int myRank = comm->getRank();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  using STS = Teuchos::ScalarTraits<SC>;
  SC zero = STS::zero(), one = STS::one();
  using real_type = typename STS::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type,LO,GO,NO>;
  const real_type realOne = Teuchos::ScalarTraits<real_type>::one();

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  GO nx = 10, ny = 10, nz = 10;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName       = "";                clp.setOption("xml",                   &xmlFileName,       "read parameters from an xml file");
  std::string yamlFileName      = "";                clp.setOption("yaml",                  &yamlFileName,      "read parameters from a yaml file");
  int         maxIts            = 200;               clp.setOption("its",                   &maxIts,            "maximum number of solver iterations");
  int         smootherIts       =  20;               clp.setOption("smootherIts",           &smootherIts,       "number of smoother iterations");
  double      smootherDamp      = 0.67;              clp.setOption("smootherDamp",          &smootherDamp,      "damping parameter for the level smoother");
  double      tol               = 1e-12;             clp.setOption("tol",                   &tol,               "solver convergence tolerance");
  bool        scaleResidualHist = true;              clp.setOption("scale", "noscale",      &scaleResidualHist, "scaled Krylov residual history");
  bool        solvePreconditioned = true;            clp.setOption("solve-preconditioned","no-solve-preconditioned", &solvePreconditioned, "use MueLu preconditioner in solve");
#ifdef HAVE_MUELU_TPETRA
  std::string equilibrate = "no" ;                   clp.setOption("equilibrate",           &equilibrate,       "equilibrate the system (no | diag | 1-norm)");
#endif
#ifdef HAVE_MUELU_CUDA
  bool profileSetup = false;                         clp.setOption("cuda-profile-setup", "no-cuda-profile-setup", &profileSetup, "enable CUDA profiling for setup");
  bool profileSolve = false;                         clp.setOption("cuda-profile-solve", "no-cuda-profile-solve", &profileSolve, "enable CUDA profiling for solve");
#endif
  int  cacheSize = 0;                                clp.setOption("cachesize",               &cacheSize,       "cache size (in KB)");

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(xmlFileName != "" && yamlFileName != "", std::runtime_error,
                             "Cannot provide both xml and yaml input files");

  // Instead of checking each time for rank, create a rank 0 stream
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out = *fancy;
  out.setOutputToRootOnly(0);

  ParameterList paramList;
  auto inst = xpetraParameters.GetInstantiation();

  if (yamlFileName != "") {
    Teuchos::updateParametersFromYamlFileAndBroadcast(yamlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);
  } else {
    if (inst == Xpetra::COMPLEX_INT_INT)
      xmlFileName = (xmlFileName != "" ? xmlFileName : "structured_1dof-complex.xml");
    else
      xmlFileName = (xmlFileName != "" ? xmlFileName : "structured_1dof.xml");
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);
  }


  // Retrieve matrix parameters (they may have been changed on the command line)
  // [for instance, if we changed matrix type from 2D to 3D we need to update nz]
  ParameterList galeriList = galeriParameters.GetParameterList();

  // =========================================================================
  // Problem construction
  // =========================================================================
  std::ostringstream galeriStream;
#ifdef HAVE_MUELU_OPENMP
  std::string node_name = Node::name();
  if(!comm->getRank() && !node_name.compare("OpenMP/Wrapper"))
    galeriStream<<"OpenMP Max Threads = "<<omp_get_max_threads()<<std::endl;
#endif


  comm->barrier();
  Teuchos::TimeMonitor::setStackedTimer(Teuchos::null);
  RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: S - Global Time")));
  RCP<TimeMonitor> tm                = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1 - Build Composite Matrix")));


  RCP<Matrix>      A;
  RCP<Xpetra::Map<LO,GO,NO> > nodeMap, dofMap;
  RCP<RealValuedMultiVector> coordinates;
  RCP<Vector> X, B;

  galeriStream << "========================================================\n" << xpetraParameters << galeriParameters;

  // Galeri will attempt to create a square-as-possible distribution of subdomains di, e.g.,
  //                                 d1  d2  d3
  //                                 d4  d5  d6
  //                                 d7  d8  d9
  //                                 d10 d11 d12
  // A perfect distribution is only possible when the #processors is a perfect square.
  // This *will* result in "strip" distribution if the #processors is a prime number or if the factors are very different in
  // size. For example, np=14 will give a 7-by-2 distribution.
  // If you don't want Galeri to do this, specify mx or my on the galeriList.
  std::string matrixType = galeriParameters.GetMatrixType();
  int numDimensions  = 0;
  int numDofsPerNode = 1;
  Teuchos::Array<GO> procsPerDim(3);
  Teuchos::Array<GO> gNodesPerDim(3);
  Teuchos::Array<LO> lNodesPerDim(3);

  // Create map and coordinates
  // In the future, we hope to be able to first create a Galeri problem, and then request map and coordinates from it
  // At the moment, however, things are fragile as we hope that the Problem uses same map and coordinates inside
  if (matrixType == "Laplace1D") {
    numDimensions = 1;
    nodeMap = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double,LO,GO,Map,RealValuedMultiVector>("1D", nodeMap, galeriList);

  } else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
             matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
    numDimensions = 2;
    nodeMap = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double,LO,GO,Map,RealValuedMultiVector>("2D", nodeMap, galeriList);

  } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
    numDimensions = 3;
    nodeMap = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double,LO,GO,Map,RealValuedMultiVector>("3D", nodeMap, galeriList);
  }

  // Expand map to do multiple DOF per node for block problems
  if (matrixType == "Elasticity2D") {
    dofMap = Xpetra::MapFactory<LO,GO,Node>::Build(nodeMap, 2);
  } else if (matrixType == "Elasticity3D") {
    dofMap = Xpetra::MapFactory<LO,GO,Node>::Build(nodeMap, 3);
  } else {
    dofMap = Xpetra::MapFactory<LO,GO,Node>::Build(nodeMap, 1);
  }

  galeriStream << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
               << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
               << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
               << "========================================================" << std::endl;

  if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
    // Our default test case for elasticity: all boundaries of a square/cube have Neumann b.c. except left which has Dirichlet
    galeriList.set("right boundary" , "Neumann");
    galeriList.set("bottom boundary", "Neumann");
    galeriList.set("top boundary"   , "Neumann");
    galeriList.set("front boundary" , "Neumann");
    galeriList.set("back boundary"  , "Neumann");
  }

  RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
    Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(galeriParameters.GetMatrixType(), dofMap, galeriList);
  A = Pr->BuildMatrix();

  if(matrixType == "Elasticity2D") {
    numDofsPerNode = 2;
  } else if(matrixType == "Elasticity3D") {
    numDofsPerNode = 3;
  }

  A->SetFixedBlockSize(numDofsPerNode);

  X = VectorFactory::Build(dofMap);
  B = VectorFactory::Build(dofMap);

  // we set seed for reproducibility
  Utilities::SetRandomSeed(*comm);
  X->randomize();
  A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

  Teuchos::Array<typename STS::magnitudeType> norms(1);
  B->norm2(norms);
  B->scale(one/norms[0]);
  galeriStream << "Galeri complete.\n========================================================" << std::endl;

  out << galeriStream.str();

  comm->barrier();
  tm = Teuchos::null;

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 2 - Compute region data")));

  // Loading geometric info from galeri
  if(numDimensions == 1) {
    gNodesPerDim[0] = galeriList.get<GO>("nx");
    gNodesPerDim[1] = 1;
    gNodesPerDim[2] = 1;

    lNodesPerDim[0] = galeriList.get<LO>("lnx");
    lNodesPerDim[1] = 1;
    lNodesPerDim[2] = 1;

    procsPerDim[0] = galeriList.get<GO>("mx");
    procsPerDim[1] = 1;
    procsPerDim[2] = 1;
  } else if(numDimensions == 2) {
    gNodesPerDim[0] = galeriList.get<GO>("nx");
    gNodesPerDim[1] = galeriList.get<GO>("ny");
    gNodesPerDim[2] = 1;

    lNodesPerDim[0] = galeriList.get<LO>("lnx");
    lNodesPerDim[1] = galeriList.get<LO>("lny");
    lNodesPerDim[2] = 1;

    procsPerDim[0] = galeriList.get<GO>("mx");
    procsPerDim[1] = galeriList.get<GO>("my");
    procsPerDim[2] = 1;
  } else if(numDimensions == 3) {
    gNodesPerDim[0] = galeriList.get<GO>("nx");
    gNodesPerDim[1] = galeriList.get<GO>("ny");
    gNodesPerDim[2] = galeriList.get<GO>("nz");

    lNodesPerDim[0] = galeriList.get<LO>("lnx");
    lNodesPerDim[1] = galeriList.get<LO>("lny");
    lNodesPerDim[2] = galeriList.get<LO>("lnz");

    procsPerDim[0] = galeriList.get<GO>("mx");
    procsPerDim[1] = galeriList.get<GO>("my");
    procsPerDim[2] = galeriList.get<GO>("mz");
  }

  Teuchos::Array<GO> startIndices(3);
  Teuchos::Array<GO> endIndices(3);
  const GO startGID = dofMap->getMinGlobalIndex();
  startIndices[2] = startGID / (gNodesPerDim[1]*gNodesPerDim[0]);
  const GO rem    = startGID % (gNodesPerDim[1]*gNodesPerDim[0]);
  startIndices[1] = rem / gNodesPerDim[0];
  startIndices[0] = rem % gNodesPerDim[0];
  endIndices[0] = startIndices[0] + lNodesPerDim[0] - 1;
  endIndices[1] = startIndices[1] + lNodesPerDim[1] - 1;
  endIndices[2] = startIndices[2] + lNodesPerDim[2] - 1;

  const LO numLocalCompositeNodes = lNodesPerDim[0]*lNodesPerDim[1]*lNodesPerDim[2];

  int leftBC = 0, rightBC = 0, frontBC = 0, backBC = 0, bottomBC = 0, topBC = 0;
  if(startIndices[0] == 0) {leftBC = 1;}
  if(startIndices[1] == 0) {frontBC = 1;}
  if(startIndices[2] == 0) {bottomBC = 1;}

  if(endIndices[0] == gNodesPerDim[0] - 1) {rightBC = 1;}
  if(endIndices[1] == gNodesPerDim[1] - 1) {backBC = 1;}
  if(endIndices[2] == gNodesPerDim[2] - 1) {topBC = 1;}

  // std::cout << "p=" << myRank << " | startGID= " << startGID
  //           << ", startIndices: " << startIndices
  //           << ", endIndices: " << endIndices
  //           << ", gNodesPerDim: " << gNodesPerDim
  //           << ", BCs={" << leftBC << ", " << rightBC << ", "
  //           << frontBC << ", " << backBC << ", "
  //           << bottomBC << ", " << topBC << "}" << std::endl;

  // Rule for boundary duplication
  // For any two ranks that share an interface:
  // the lowest ranks owns the interface and the highest rank gets extra nodes

  // First we count how many nodes the region needs to send and receive
  // and allocate arrays accordingly
  const int maxRegPerProc = 1;
  int maxRegPerGID = 0;
  LO numReceive = 0, numSend = 0;
  Teuchos::Array<GO>  receiveGIDs;
  Teuchos::Array<int> receivePIDs;
  Teuchos::Array<LO>  receiveLIDs;
  Teuchos::Array<GO>  sendGIDs;
  Teuchos::Array<int> sendPIDs;
  Teuchos::Array<LO>  sendLIDs;
  if(numDimensions == 1) {
    maxRegPerGID = 2;
    if(leftBC == 0) {
      numReceive = 1;
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      receiveGIDs[0] = startIndices[0] - 1;
      receivePIDs[0] = myRank - 1;
    }
    if(rightBC == 0) {
      numSend = 1;
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      sendGIDs[0] = endIndices[0];
      sendGIDs[0] = myRank + 1;
      sendLIDs[0] = lNodesPerDim[0] - 1;
    }
  } else if(numDimensions == 2) {
    maxRegPerGID = 4;
    // Received nodes
    if(frontBC == 0 && leftBC == 0) {
      numReceive = lNodesPerDim[0] + lNodesPerDim[1] + 1;
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front-left corner node
      receiveGIDs[countIDs] = startGID - gNodesPerDim[0] - 1;
      receivePIDs[countIDs] = myRank - procsPerDim[0] - 1;
      ++countIDs;
      // Receive front edge nodes
      for(LO i = 0; i < lNodesPerDim[0]; ++i) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] + i;
        receivePIDs[countIDs] = myRank - procsPerDim[0];
        ++countIDs;
      }
      // Receive left edge nodes
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        receiveGIDs[countIDs] = startGID - 1 + j*gNodesPerDim[0];
        receivePIDs[countIDs] = myRank - 1;
        ++countIDs;
      }
    } else if(frontBC == 0) {
      numReceive = lNodesPerDim[0];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front edge nodes
      for(LO i = 0; i < lNodesPerDim[0]; ++i) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] + i;
        receivePIDs[countIDs] = myRank - procsPerDim[0];
        ++countIDs;
      }
    } else if(leftBC == 0) {
      numReceive = lNodesPerDim[1];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive left edge nodes
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        receiveGIDs[countIDs] = startGID - 1 + j*gNodesPerDim[0];
        receivePIDs[countIDs] = myRank - 1;
        ++countIDs;
      }
    }

    // Sent nodes
    if(rightBC == 0 && backBC == 0) {
      numSend = lNodesPerDim[0] + lNodesPerDim[1] + 1;
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right edge
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        sendGIDs[countIDs] = j*gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + 1;
        sendLIDs[countIDs] = (j + 1)*lNodesPerDim[0] - 1;
        ++countIDs;
      }
      // Send nodes of back edge
      for(LO i = 0; i < lNodesPerDim[0]; ++i) {
        sendGIDs[countIDs] = i + startGID + (lNodesPerDim[1] - 1)*gNodesPerDim[0];
        sendPIDs[countIDs] = myRank + procsPerDim[0];
        sendLIDs[countIDs] = (lNodesPerDim[1] - 1)*lNodesPerDim[0] + i;
        ++countIDs;
      }
      // Send node of back-right corner
      sendGIDs[countIDs] = startGID + (lNodesPerDim[1] - 1)*gNodesPerDim[0] + lNodesPerDim[0] - 1;
      sendPIDs[countIDs] = myRank + procsPerDim[1] + 1;
      sendLIDs[countIDs] = lNodesPerDim[1]*lNodesPerDim[0] - 1;
      ++countIDs;
    } else if(backBC == 0) {
      numSend = lNodesPerDim[0];

      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of back edge
      for(LO i = 0; i < lNodesPerDim[0]; ++i) {
        sendGIDs[countIDs] = i + startGID + (lNodesPerDim[1] - 1)*gNodesPerDim[0];
        sendPIDs[countIDs] = myRank + procsPerDim[0];
        sendLIDs[countIDs] = (lNodesPerDim[1] - 1)*lNodesPerDim[0] + i;
        ++countIDs;
      }
    } else if(rightBC == 0) {
      numSend = lNodesPerDim[1];
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right edge
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        sendGIDs[countIDs] = j*gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + 1;
        sendLIDs[countIDs] = (j + 1)*lNodesPerDim[0] - 1;
        ++countIDs;
      }
    }
  } else if(numDimensions == 3) {
    maxRegPerGID = 8;
    // Received nodes
    if( (bottomBC == 0) && (frontBC == 0) && (leftBC == 0) ) {
      numReceive = lNodesPerDim[0]*lNodesPerDim[1]     // bottom face
        + (lNodesPerDim[0] + 1)*lNodesPerDim[2]        // front face
        + (lNodesPerDim[1] + 1)*(lNodesPerDim[2] + 1); // left face
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front-left-bottom corner node
      receiveGIDs[countIDs] = startGID - gNodesPerDim[0] - 1
          - gNodesPerDim[1]*gNodesPerDim[0];
      receivePIDs[countIDs] = myRank - procsPerDim[0] - 1
          - procsPerDim[1]*procsPerDim[0];
      ++countIDs;
      // Receive front-bottom edge nodes
      for(LO i = 0; i < lNodesPerDim[0]; ++i) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0]*gNodesPerDim[1]
              - gNodesPerDim[0] + i;
        receivePIDs[countIDs] = myRank - procsPerDim[0]*procsPerDim[1] - procsPerDim[0];
        ++countIDs;
      }
      // Receive left-bottom edge nodes
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0]*gNodesPerDim[1]
              - 1 + j*gNodesPerDim[0];
        receivePIDs[countIDs] = myRank - procsPerDim[0]*procsPerDim[1] - 1;
        ++countIDs;
      }
      // Receive bottom face nodes
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = startGID - gNodesPerDim[0]*gNodesPerDim[1]
              + i
              + j*gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[0]*procsPerDim[1];
          ++countIDs;
        }
      }
      // Receive front-left edge nodes
      for(LO k = 0; k < lNodesPerDim[1]; ++k) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0]
              - 1 + k*gNodesPerDim[0]*gNodesPerDim[1];
        receivePIDs[countIDs] = myRank - procsPerDim[0] - 1;
        ++countIDs;
      }
      // Receive front face nodes
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = startGID - gNodesPerDim[0] + i
              + k*(gNodesPerDim[1]*gNodesPerDim[0]);
          receivePIDs[countIDs] = myRank - procsPerDim[0];
          ++countIDs;
        }
      }
      // Receive left face nodes
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO j = 0; j < lNodesPerDim[1]; ++j) {
          receiveGIDs[countIDs] = startGID - 1
              + j*gNodesPerDim[0]
              + k*(gNodesPerDim[1]*gNodesPerDim[0]);
          receivePIDs[countIDs] = myRank - 1;
          ++countIDs;
        }
      }

    // Two faces received
    } else if( (bottomBC == 0) && (frontBC == 0) ) {
      numReceive = lNodesPerDim[0]*lNodesPerDim[1]     // bottom face
        + lNodesPerDim[0]*(lNodesPerDim[2] + 1);       // front face;
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front-bottom edge nodes
      for(LO i = 0; i < lNodesPerDim[0]; ++i) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0]*gNodesPerDim[1]
              - gNodesPerDim[0] + i;
        receivePIDs[countIDs] = myRank - procsPerDim[0]*procsPerDim[1] - procsPerDim[0];
        ++countIDs;
      }
      // Receive bottom face nodes
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = startGID - gNodesPerDim[0]*gNodesPerDim[1]
              + i
              + j*gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[0]*procsPerDim[1];
          ++countIDs;
        }
      }
      // Receive front face nodes
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = startGID - gNodesPerDim[0] + i
              + k*(gNodesPerDim[1]*gNodesPerDim[0]);
          receivePIDs[countIDs] = myRank - procsPerDim[0];
          ++countIDs;
        }
      }

    } else if( (bottomBC == 0) && (leftBC == 0) ) {
      numReceive = lNodesPerDim[1]*(lNodesPerDim[0] + lNodesPerDim[2] + 1);
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive left-bottom edge nodes
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        receiveGIDs[countIDs] =  j*gNodesPerDim[0]
          + startGID - gNodesPerDim[1]*gNodesPerDim[0] - 1;
        receivePIDs[countIDs] =  myRank - procsPerDim[1]*procsPerDim[0] - 1;
        ++countIDs;
      }
      // Receive bottom face nodes
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = j*gNodesPerDim[0] + i
            + startGID - gNodesPerDim[1]*gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[1]*procsPerDim[0];
          ++countIDs;
        }
      }
      // Receive left face nodes
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO j = 0; j < lNodesPerDim[1]; ++j) {
          receiveGIDs[countIDs] = k*gNodesPerDim[1]*gNodesPerDim[0] + j*gNodesPerDim[0]
            + startGID - 1;
          receivePIDs[countIDs] = myRank - 1;
          ++countIDs;
        }
      }

    } else if( (frontBC == 0) && (leftBC == 0) ) {
      numReceive = lNodesPerDim[2]*(lNodesPerDim[1] + lNodesPerDim[0] + 1);
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front-left edge nodes
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        receiveGIDs[countIDs] =  k*gNodesPerDim[1]*gNodesPerDim[0]
          + startGID - gNodesPerDim[0] - 1;
        receivePIDs[countIDs] =  myRank - procsPerDim[0] - 1;
        ++countIDs;
      }
      // Receive front face nodes
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = k*gNodesPerDim[1]*gNodesPerDim[0] + i
            + startGID - gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[0];
          ++countIDs;
        }
      }
      // Receive left face nodes
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO j = 0; j < lNodesPerDim[1]; ++j) {
          receiveGIDs[countIDs] = k*gNodesPerDim[1]*gNodesPerDim[0] + j*gNodesPerDim[0]
            + startGID - 1;
          receivePIDs[countIDs] = myRank - 1;
          ++countIDs;
        }
      }

    // Single face received
    } else if(bottomBC == 0) {
      numReceive = lNodesPerDim[0]*lNodesPerDim[1];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive bottom face nodes
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = j*gNodesPerDim[0] + i
            + startGID - gNodesPerDim[1]*gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[1]*procsPerDim[0];
          ++countIDs;
        }
      }

    } else if(frontBC == 0) {
      numReceive = lNodesPerDim[0]*lNodesPerDim[2];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front face nodes
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = k*gNodesPerDim[1]*gNodesPerDim[0] + i
            + startGID - gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[0];
          ++countIDs;
        }
      }

    } else if(leftBC == 0) {
      numReceive = lNodesPerDim[1]*lNodesPerDim[2];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Recive left face nodes
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO j = 0; j < lNodesPerDim[1]; ++j) {
          receiveGIDs[countIDs] = k*gNodesPerDim[1]*gNodesPerDim[0]
            + j*gNodesPerDim[0] + startGID - 1;
          receivePIDs[countIDs] = myRank - 1;
          ++countIDs;
        }
      }

    }

    // Sent nodes
    if( (topBC == 0) && (backBC == 0) && (rightBC == 0) ) {
      numSend = (lNodesPerDim[0])*(lNodesPerDim[1])
        + (lNodesPerDim[0])*(lNodesPerDim[2])
        + (lNodesPerDim[1])*(lNodesPerDim[2])
        + lNodesPerDim[0]
        + lNodesPerDim[1]
        + lNodesPerDim[2]
        + 1;
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right face
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO j = 0; j < lNodesPerDim[1]; ++j) {
          sendGIDs[countIDs] = k*(gNodesPerDim[1]*gNodesPerDim[0])
            + j*gNodesPerDim[0]
            + startGID + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + 1;
          ++countIDs;
        }
      }
      // Send nodes of back face
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = k*(gNodesPerDim[1]*gNodesPerDim[0]) + i
            + startGID + (lNodesPerDim[1] - 1)*gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of right-back edge
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        sendGIDs[countIDs] = k*(gNodesPerDim[1]*gNodesPerDim[0])
          + startGID + (lNodesPerDim[1] - 1)*gNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + procsPerDim[0] + 1;
        ++countIDs;
      }
      // Send nodes of top face
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = j*gNodesPerDim[0] + i
            + startGID + (lNodesPerDim[2] - 1)*gNodesPerDim[1]*gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[1]*procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of top-right edge
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        sendGIDs[countIDs] = j*gNodesPerDim[0]
          + startGID + (lNodesPerDim[2] - 1)*gNodesPerDim[1]*gNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + procsPerDim[1]*procsPerDim[0] + 1;
        ++countIDs;
      }
      // Send nodes of top-back edge
      for(LO i = 0; i < lNodesPerDim[0]; ++i) {
        sendGIDs[countIDs] = i
          + startGID + (lNodesPerDim[2] - 1)*gNodesPerDim[1]*gNodesPerDim[0]
          + (lNodesPerDim[0] - 1)*gNodesPerDim[0];
        sendPIDs[countIDs] = myRank + procsPerDim[1]*procsPerDim[0] + procsPerDim[1];
        ++countIDs;
      }
      // Send node of top-back-right corner
      sendGIDs[countIDs] = startGID + (lNodesPerDim[2] - 1)*gNodesPerDim[1]*gNodesPerDim[0]
        + (lNodesPerDim[0] - 1)*gNodesPerDim[0] + lNodesPerDim[0] - 1;
      sendPIDs[countIDs] = myRank + procsPerDim[1]*procsPerDim[0] + procsPerDim[1] + 1;
      ++countIDs;

    } else if( (topBC == 0) && (backBC == 0) ) {
      numSend = (lNodesPerDim[0]*lNodesPerDim[2]) // back face
              + (lNodesPerDim[0]*lNodesPerDim[1]) // Top face
              + (lNodesPerDim[0]); // top-back edge
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of back face
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = k*(gNodesPerDim[1]*gNodesPerDim[0]) + i
            + startGID + (lNodesPerDim[1] - 1)*gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of top face
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = j*gNodesPerDim[0] + i
            + startGID + (lNodesPerDim[2] - 1)*gNodesPerDim[1]*gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[1]*procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of top-back edge
      for(LO i = 0; i < lNodesPerDim[0]; ++i) {
        sendGIDs[countIDs] = i
          + startGID + (lNodesPerDim[2] - 1)*gNodesPerDim[1]*gNodesPerDim[0]
          + (lNodesPerDim[0] - 1)*gNodesPerDim[0];
        sendPIDs[countIDs] = myRank + procsPerDim[1]*procsPerDim[0] + procsPerDim[1];
        ++countIDs;
      }

    } else if( (topBC == 0) && (rightBC == 0) ) {
      numSend = (lNodesPerDim[1]*lNodesPerDim[2]) // right face
              + (lNodesPerDim[0]*lNodesPerDim[1]) // Top face
              + (lNodesPerDim[1]); // top-right edge
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right face
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO j = 0; j < lNodesPerDim[1]; ++j) {
          sendGIDs[countIDs] = k*(gNodesPerDim[1]*gNodesPerDim[0])
            + j*gNodesPerDim[0]
            + startGID + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + 1;
          ++countIDs;
        }
      }
      // Send nodes of top face
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = j*gNodesPerDim[0] + i
            + startGID + (lNodesPerDim[2] - 1)*gNodesPerDim[1]*gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[1]*procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of top-right edge
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        sendGIDs[countIDs] = j*gNodesPerDim[0]
          + startGID + (lNodesPerDim[2] - 1)*gNodesPerDim[1]*gNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + procsPerDim[1]*procsPerDim[0] + 1;
        ++countIDs;
      }

    } else if( (backBC == 0) && (rightBC == 0) ) {
      numSend = lNodesPerDim[2]*(lNodesPerDim[0] + lNodesPerDim[1] + 1);
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right face
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO j = 0; j < lNodesPerDim[1]; ++j) {
          sendGIDs[countIDs] = k*gNodesPerDim[1]*gNodesPerDim[0] + j*gNodesPerDim[0]
            + startGID + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + 1;
          ++countIDs;
        }
      }
      // Send nodes of back face
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = k*gNodesPerDim[1]*gNodesPerDim[0] + i
            + startGID + (lNodesPerDim[1] - 1)*gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of back-right edge
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
          sendGIDs[countIDs] = k*gNodesPerDim[1]*gNodesPerDim[0]
            + startGID + (lNodesPerDim[1] - 1)*gNodesPerDim[0] + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + procsPerDim[0] + 1;
          ++countIDs;
      }

    } else if(topBC == 0) {
      numSend = lNodesPerDim[0]*lNodesPerDim[1];
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of top face
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = j*gNodesPerDim[0] + i
            + startGID + (lNodesPerDim[2] - 1)*gNodesPerDim[0]*gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[1]*procsPerDim[0];
          ++countIDs;
        }
      }

    } else if(backBC == 0) {
      numSend = lNodesPerDim[0]*lNodesPerDim[2];
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of back face
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = k*gNodesPerDim[1]*gNodesPerDim[0] + i
            + startGID + (lNodesPerDim[1] - 1)*gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[0];
          ++countIDs;
        }
      }

    } else if(rightBC == 0) {
      numSend = lNodesPerDim[1]*lNodesPerDim[2];
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right face
      for(LO k = 0; k < lNodesPerDim[2]; ++k) {
        for(LO j = 0; j < lNodesPerDim[1]; ++j) {
          sendGIDs[countIDs] = k*gNodesPerDim[1]*gNodesPerDim[0]
            + j*gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + 1;
          ++countIDs;
        }
      }

    }
  }

  // Second we actually fill the send and receive arrays with appropriate data
  // which will allow us to compute the region and composite maps.

  // std::cout << "p=" << myRank << " | numReceive=" << numReceive
  //           << ", numSend=" << numSend << std::endl;
  // std::cout << "p=" << myRank << " | receiveGIDs: " << receiveGIDs << std::endl;
  // std::cout << "p=" << myRank << " | receivePIDs: " << receivePIDs << std::endl;
  // std::cout << "p=" << myRank << " | sendGIDs: " << sendGIDs << std::endl;
  // std::cout << "p=" << myRank << " | sendPIDs: " << sendPIDs << std::endl;

  // Now we can construct a list of GIDs that corresponds to rowMapPerGrp
  const LO numLocalRegionNodes = numLocalCompositeNodes + numReceive;
  Array<LO> rNodesPerDim(3);
  Array<GO> quasiRegionGIDs(numLocalRegionNodes*numDofsPerNode);
  Array<GO> regionGIDs(numLocalRegionNodes*numDofsPerNode);

  rNodesPerDim[0] = lNodesPerDim[0];
  rNodesPerDim[1] = lNodesPerDim[1];
  rNodesPerDim[2] = lNodesPerDim[2];
  if(leftBC   == 0) {rNodesPerDim[0] += 1;}
  if(frontBC  == 0) {rNodesPerDim[1] += 1;}
  if(bottomBC == 0) {rNodesPerDim[2] += 1;}

  // std::cout << "p=" << myRank << " | numLocalRegionNodes=" << numLocalRegionNodes
  //           << ", rNodesPerDim: " << rNodesPerDim << std::endl;

  Array<LO> compositeToRegionLIDs(nodeMap->getNodeNumElements());

  // Using receiveGIDs, rNodesPerDim and numLocalRegionNodes, build quasi-region row map
  // This will potentially be done by the application or in a MueLu interface but for now
  // doing it in the driver seem to avoid design hassle...
  LO interfaceCount = 0, compositeIdx = 0;
  Array<LO> regionIJK(3);
  for(LO nodeRegionIdx = 0; nodeRegionIdx < numLocalRegionNodes; ++nodeRegionIdx) {
    regionIJK[2] = nodeRegionIdx / (rNodesPerDim[1]*rNodesPerDim[0]);
    LO tmp       = nodeRegionIdx % (rNodesPerDim[1]*rNodesPerDim[0]);
    regionIJK[1] = tmp / rNodesPerDim[0];
    regionIJK[0] = tmp % rNodesPerDim[0];

    // std::cout << "p=" << myRank << " | regionIJK=" << regionIJK << std::endl;

    if( (regionIJK[0] == 0 && leftBC   == 0) ||
        (regionIJK[1] == 0 && frontBC  == 0) ||
        (regionIJK[2] == 0 && bottomBC == 0) ) {
      for(int dof = 0; dof < numDofsPerNode; ++dof) {
        quasiRegionGIDs[nodeRegionIdx*numDofsPerNode + dof] =
          receiveGIDs[interfaceCount]*numDofsPerNode + dof;
      }
      receiveLIDs[interfaceCount] = nodeRegionIdx;
      ++interfaceCount;
    } else {
      compositeIdx = (regionIJK[2] + bottomBC - 1)*lNodesPerDim[1]*lNodesPerDim[0]
        + (regionIJK[1] + frontBC  - 1)*lNodesPerDim[0]
        + (regionIJK[0] + leftBC - 1);
      for(int dof = 0; dof < numDofsPerNode; ++dof) {
        quasiRegionGIDs[nodeRegionIdx*numDofsPerNode + dof]
          = dofMap->getGlobalElement(compositeIdx*numDofsPerNode + dof);
      }
      compositeToRegionLIDs[compositeIdx] = nodeRegionIdx;
    }
  }

  // std::cout << "p=" << myRank << " | receiveLIDs: " << receiveLIDs() << std::endl;
  // std::cout << "p=" << myRank << " | sendLIDs: " << sendLIDs() << std::endl;
  // std::cout << "p=" << myRank << " | compositeToRegionLIDs: " << compositeToRegionLIDs() << std::endl;
  // std::cout << "p=" << myRank << " | quasiRegionGIDs: " << quasiRegionGIDs << std::endl;

  // In our very particular case we know that a node is at most shared by 4 regions.
  // Other geometries will certainly have different constrains and a parallel reduction using MAX
  // would be appropriate.
  RCP<Xpetra::MultiVector<LO, LO, GO, NO> > regionsPerGID
    = Xpetra::MultiVectorFactory<LO, LO, GO, NO>::Build(dofMap, maxRegPerGID, false);

  { // Scope for regionsPerGIDView
    Array<ArrayRCP<LO> > regionsPerGIDView(maxRegPerGID);
    for(int regionIdx = 0; regionIdx < maxRegPerGID; ++regionIdx) {
      regionsPerGIDView[regionIdx] = regionsPerGID->getDataNonConst(regionIdx);
    }

    // Initialize all entries to myRank in first column and to -1 in other columns
    for(LO dofIdx = 0; dofIdx < numLocalCompositeNodes*numDofsPerNode; ++dofIdx) {
      regionsPerGIDView[0][dofIdx] = myRank;
      for(int regionIdx = 1; regionIdx < maxRegPerGID; ++regionIdx) {
        regionsPerGIDView[regionIdx][dofIdx] = -1;
      }
    }

    // Now loop over the sendGIDs array to fill entries with values in sendPIDs
    LO nodeIdx = 0;
    for(LO sendIdx = 0; sendIdx < numSend; ++sendIdx) {
      nodeIdx = nodeMap->getLocalElement(sendGIDs[sendIdx]);
      for(int dof = 0; dof < numDofsPerNode; ++dof) {
        LO dofIdx = nodeIdx*numDofsPerNode + dof;
        for(int regionIdx = 1; regionIdx < maxRegPerGID; ++regionIdx) {
          if(regionsPerGIDView[regionIdx][dofIdx] == -1) {
            regionsPerGIDView[regionIdx][dofIdx] = sendPIDs[sendIdx];
            break;
          }
        }
      }
    }
  } // end of regionsPerGIDView's scope

  // if(myRank == 0) std::cout << "regionsPerGID:" << std::endl;
  // regionsPerGID->describe(*Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)), Teuchos::VERB_EXTREME);

  comm->barrier();
  tm = Teuchos::null;

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3 - Build Region Matrix")));

  std::vector<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > rowMapPerGrp(maxRegPerProc), colMapPerGrp(maxRegPerProc);
  std::vector<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > revisedRowMapPerGrp(maxRegPerProc), revisedColMapPerGrp(maxRegPerProc);
  rowMapPerGrp[0] = Xpetra::MapFactory<LO,GO,Node>::Build(dofMap->lib(),
                                                          Teuchos::OrdinalTraits<GO>::invalid(),
                                                          quasiRegionGIDs(),
                                                          dofMap->getIndexBase(),
                                                          dofMap->getComm());
  colMapPerGrp[0] = rowMapPerGrp[0];
  revisedRowMapPerGrp[0] = Xpetra::MapFactory<LO,GO,Node>::Build(dofMap->lib(),
                                                                 Teuchos::OrdinalTraits<GO>::invalid(),
                                                                 numLocalRegionNodes,
                                                                 dofMap->getIndexBase(),
                                                                 dofMap->getComm());
  revisedColMapPerGrp[0] = revisedRowMapPerGrp[0];

  // Setup importers
  std::vector<RCP<Import> > rowImportPerGrp(maxRegPerProc);
  std::vector<RCP<Import> > colImportPerGrp(maxRegPerProc);
  rowImportPerGrp[0] = ImportFactory::Build(dofMap, rowMapPerGrp[0]);
  colImportPerGrp[0] = ImportFactory::Build(dofMap, colMapPerGrp[0]);

  RCP<Xpetra::MultiVector<LO, LO, GO, NO> > regionsPerGIDWithGhosts
    = Xpetra::MultiVectorFactory<LO, LO, GO, NO>::Build(A->getColMap(), maxRegPerGID, false);
  RCP<Import> regionsPerGIDImport = ImportFactory::Build(A->getRowMap(), A->getColMap());
  regionsPerGIDWithGhosts->doImport(*regionsPerGID, *regionsPerGIDImport, Xpetra::INSERT);

  // if(myRank == 0) std::cout << "regionsPerGIDWithGhosts:" << std::endl;
  // regionsPerGIDWithGhosts->describe(*Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)), Teuchos::VERB_EXTREME);

  std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > quasiRegionGrpMats(1);
  MakeQuasiregionMatrices(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A), maxRegPerProc,
                          regionsPerGIDWithGhosts, rowMapPerGrp, colMapPerGrp, rowImportPerGrp,
                          quasiRegionGrpMats);

  std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionGrpMats(1);
  MakeRegionMatrices(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A), A->getRowMap(), rowMapPerGrp,
                     revisedRowMapPerGrp, revisedColMapPerGrp,
                     rowImportPerGrp, maxRegPerProc, quasiRegionGrpMats, regionGrpMats);

  // if(myRank == 0) std::cout << "composite A:" << std::endl;
  // A->describe(*Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)), Teuchos::VERB_EXTREME);

  // if(myRank == 0) std::cout << "quasi-region A:" << std::endl;
  // quasiRegionGrpMats[0]->describe(*Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)), Teuchos::VERB_EXTREME);

  // if(myRank == 0) std::cout << "region A:" << std::endl;
  // regionGrpMats[0]->describe(*Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)), Teuchos::VERB_EXTREME);

  comm->barrier();
  tm = Teuchos::null;

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 4 - Build Region Hierarchy")));

  // Setting up parameters before hierarchy construction
  // These need to stay in the driver as they would be provide by an app
  Array<Array<int> > regionNodesPerDim(maxRegPerProc);
  Array<std::string> aggregationRegionType(maxRegPerProc);
  Array<RCP<MultiVector> > nullspace(maxRegPerProc);
  Array<RCP<RealValuedMultiVector> > regionCoordinates(maxRegPerProc);

  // Set mesh structure data
  regionNodesPerDim[0] = rNodesPerDim;

  // Set aggregation type for each region
  aggregationRegionType[0] = "structured";

  // create nullspace vector
  nullspace[0] = MultiVectorFactory::Build(revisedRowMapPerGrp[0], 1);
  nullspace[0]->putScalar(one);

  // create region coordinates vector
  regionCoordinates[0] = Xpetra::MultiVectorFactory<real_type,LO,GO,NO>::Build(revisedRowMapPerGrp[0],
                                                                               numDimensions);
  using Tpetra_CrsMatrix = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Tpetra_MultiVector = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  /* Stuff for multi-level algorithm
   *
   * To allow for multi-level schemes with more than two levels, we need to store
   * maps, matrices, vectors, and stuff like that on each level. Since we call the
   * multi-level scheme recursively, this should be reflected in the design of
   * variables.
   *
   * We use Teuchos::Array<T> to store each quantity on each level.
   */
  Array<RCP<Xpetra::Map<LO,GO,NO> > > compRowMaps; // composite row maps on each level
  Array<RCP<Xpetra::Map<LO,GO,NO> > > compColMaps; // composite columns maps on each level
  Array<std::vector<RCP<Xpetra::Map<LO,GO,NO> > > > regRowMaps; // regional row maps on each level
  Array<std::vector<RCP<Xpetra::Map<LO,GO,NO> > > > regColMaps; // regional column maps on each level
  Array<std::vector<RCP<Xpetra::Map<LO,GO,NO> > > > quasiRegRowMaps; // quasiRegional row maps on each level
  Array<std::vector<RCP<Xpetra::Map<LO,GO,NO> > > > quasiRegColMaps; // quasiRegional column maps on each level
  Array<std::vector<RCP<Matrix> > > regMatrices; // regional matrices on each level
  Array<std::vector<RCP<Matrix> > > regProlong; // regional prolongators on each level
  Array<std::vector<RCP<Import> > > regRowImporters; // regional row importers on each level
  Array<std::vector<RCP<Vector> > > regInterfaceScalings; // regional interface scaling factors on each level
  Teuchos::RCP<Matrix> coarseCompOp = Teuchos::null;
  RCP<Amesos2::Solver<Tpetra_CrsMatrix, Tpetra_MultiVector> > coarseCompositeDirectSolver = Teuchos::null;


  // Create multigrid hierarchy
  createRegionHierarchy(maxRegPerProc,
                        numDimensions,
                        regionNodesPerDim,
                        aggregationRegionType,
                        xmlFileName,
                        nullspace,
                        regionCoordinates,
                        regionGrpMats,
                        dofMap,
                        rowMapPerGrp,
                        colMapPerGrp,
                        revisedRowMapPerGrp,
                        revisedColMapPerGrp,
                        rowImportPerGrp,
                        compRowMaps,
                        compColMaps,
                        regRowMaps,
                        regColMaps,
                        quasiRegRowMaps,
                        quasiRegColMaps,
                        regMatrices,
                        regProlong,
                        regRowImporters,
                        regInterfaceScalings,
                        coarseCompOp,
                        maxRegPerGID,
                        compositeToRegionLIDs(),
                        coarseCompositeDirectSolver);


  comm->barrier();
  tm = Teuchos::null;

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 5 - Solve with V-cycle")));

  {
    std::cout << myRank << " | Running V-cycle ..." << std::endl;

    // Extract the number of levels from the prolongator data structure
    int numLevels = regProlong.size();

    TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels>0), "We require numLevel > 0. Probably, numLevel has not been set, yet.");

    /* We first use the non-level container variables to setup the fine grid problem.
     * This is ok since the initial setup just mimics the application and the outer
     * Krylov method.
     *
     * We switch to using the level container variables as soon as we enter the
     * recursive part of the algorithm.
     */

    // residual vector
    RCP<Vector> compRes = VectorFactory::Build(dofMap, true);
    {
      A->apply(*X, *compRes, Teuchos::NO_TRANS);
      compRes->update(one, *B, -one);
    }

    // transform composite vectors to regional layout
    std::vector<Teuchos::RCP<Vector> > quasiRegX(maxRegPerProc);
    std::vector<Teuchos::RCP<Vector> > regX(maxRegPerProc);
    compositeToRegional(X, quasiRegX, regX, maxRegPerProc, rowMapPerGrp,
                        revisedRowMapPerGrp, rowImportPerGrp);

    std::vector<RCP<Vector> > quasiRegB(maxRegPerProc);
    std::vector<RCP<Vector> > regB(maxRegPerProc);
    compositeToRegional(B, quasiRegB, regB, maxRegPerProc, rowMapPerGrp,
                        revisedRowMapPerGrp, rowImportPerGrp);

    //    printRegionalObject<Vector>("regB 0", regB, myRank, *fos);

    std::vector<RCP<Vector> > regRes(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) { // step 1
      regRes[j] = VectorFactory::Build(revisedRowMapPerGrp[j], true);
    }

    /////////////////////////////////////////////////////////////////////////
    // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
    /////////////////////////////////////////////////////////////////////////

    // define max iteration counts
    const int maxFineIter = 20;
    const int maxCoarseIter = 100;

    // Prepare output of residual norm to file
    RCP<std::ofstream> log;
    if (myRank == 0)
      {
        std::string s = "residual_norm.txt";
        log = rcp(new std::ofstream(s.c_str()));
        (*log) << "# num procs = " << dofMap->getComm()->getSize() << "\n"
               << "# iteration | res-norm\n"
               << "#\n";
      }

    // Richardson iterations
    for (int cycle = 0; cycle < maxIts; ++cycle) {
      // check for convergence
      {
        ////////////////////////////////////////////////////////////////////////
        // SWITCH BACK TO NON-LEVEL VARIABLES
        ////////////////////////////////////////////////////////////////////////
        {
          computeResidual(regRes, regX, regB, regionGrpMats, dofMap,
              rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

          scaleInterfaceDOFs(regRes, regInterfaceScalings[0], true);
        }

        compRes = VectorFactory::Build(dofMap, true);
        regionalToComposite(regRes, compRes, maxRegPerProc, rowMapPerGrp,
                            rowImportPerGrp, Xpetra::ADD);
        typename Teuchos::ScalarTraits<Scalar>::magnitudeType normRes = compRes->norm2();

        // Output current residual norm to screen (on proc 0 only)
        if (myRank == 0)
          {
            std::cout << cycle << "\t" << normRes << std::endl;
            (*log) << cycle << "\t" << normRes << "\n";
          }

        if (normRes < tol)
          break;
      }

      /////////////////////////////////////////////////////////////////////////
      // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
      /////////////////////////////////////////////////////////////////////////

      //      printRegionalObject<Vector>("regB 2", regB, myRank, *fos);
      vCycle(0, numLevels, smootherIts, maxCoarseIter, smootherDamp, maxRegPerProc,
             regX, regB, regMatrices,
             regProlong, compRowMaps, quasiRegRowMaps, regRowMaps, regRowImporters,
             regInterfaceScalings, coarseCompOp, coarseCompositeDirectSolver);
    }
  }

  comm->barrier();
  tm = Teuchos::null;
  globalTimeMonitor = Teuchos::null;

  RCP<ParameterList> reportParams = rcp(new ParameterList);
  const std::string filter = "";
  std::ios_base::fmtflags ff(out.flags());
  TimeMonitor::report(comm.ptr(), out, filter, reportParams);
  out << std::setiosflags(ff);

  TimeMonitor::clearCounters();

  return EXIT_SUCCESS;
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}
