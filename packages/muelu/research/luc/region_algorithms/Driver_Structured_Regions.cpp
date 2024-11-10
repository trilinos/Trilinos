// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include <MueLu.hpp>

#include <MueLu_BaseClass.hpp>
#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif
#include <MueLu_Level.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosBiCGStabSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>  // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>   // => This header defines Belos::MueLuOp
#include <BelosTpetraAdapter.hpp>  // => This header defines Belos::TpetraOp
#endif

#ifdef HAVE_MUELU_CUDA
#include "cuda_profiler_api.h"
#endif

#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include "Xpetra_TpetraMultiVector.hpp"
#include <KokkosBlas1_abs.hpp>
#include <Tpetra_leftAndOrRightScaleCrsMatrix.hpp>
#include <Tpetra_computeRowAndColumnOneNorms.hpp>

#ifdef HAVE_MUELU_EPETRA
#include "Xpetra_EpetraMultiVector.hpp"
#endif

#include <MueLu_CreateXpetraPreconditioner.hpp>

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  const int myRank                    = comm->getRank();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero(), one = STS::one();
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  GO nx = 10, ny = 10, nz = 10;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

  std::string xmlFileName = "";
  clp.setOption("xml", &xmlFileName, "read parameters from an xml file");
  std::string yamlFileName = "";
  clp.setOption("yaml", &yamlFileName, "read parameters from a yaml file");
  int maxIts = 200;
  clp.setOption("its", &maxIts, "maximum number of solver iterations");
  double tol = 1e-12;
  clp.setOption("tol", &tol, "solver convergence tolerance");
  bool scaleResidualHist = true;
  clp.setOption("scale", "noscale", &scaleResidualHist, "scaled Krylov residual history");
  bool solvePreconditioned = true;
  clp.setOption("solve-preconditioned", "no-solve-preconditioned", &solvePreconditioned, "use MueLu preconditioner in solve");
  std::string equilibrate = "no";
  clp.setOption("equilibrate", &equilibrate, "equilibrate the system (no | diag | 1-norm)");
#ifdef HAVE_MUELU_CUDA
  bool profileSetup = false;
  clp.setOption("cuda-profile-setup", "no-cuda-profile-setup", &profileSetup, "enable CUDA profiling for setup");
  bool profileSolve = false;
  clp.setOption("cuda-profile-solve", "no-cuda-profile-solve", &profileSolve, "enable CUDA profiling for solve");
#endif
  int cacheSize = 0;
  clp.setOption("cachesize", &cacheSize, "cache size (in KB)");

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(xmlFileName != "" && yamlFileName != "", std::runtime_error,
                             "Cannot provide both xml and yaml input files");

  // Instead of checking each time for rank, create a rank 0 stream
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream &out       = *fancy;
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
  if (!comm->getRank() && !node_name.compare("OpenMP/Wrapper"))
    galeriStream << "OpenMP Max Threads = " << omp_get_max_threads() << std::endl;
#endif

  comm->barrier();
  Teuchos::TimeMonitor::setStackedTimer(Teuchos::null);
  RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: S - Global Time")));
  RCP<TimeMonitor> tm                = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1 - Matrix Build")));

  RCP<Matrix> A;
  RCP<const Map> map;
  RCP<RealValuedMultiVector> coordinates;
  typedef typename RealValuedMultiVector::scalar_type Real;
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > nullspace;
  RCP<MultiVector> X, B;

  galeriStream << "========================================================\n"
               << xpetraParameters << galeriParameters;

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
  int numDimensions      = 0;
  int numDofsPerNode     = 1;
  Teuchos::Array<GO> procsPerDim(3);
  Teuchos::Array<GO> gNodesPerDim(3);
  Teuchos::Array<LO> lNodesPerDim(3);

  // Create map and coordinates
  // In the future, we hope to be able to first create a Galeri problem, and then request map and coordinates from it
  // At the moment, however, things are fragile as we hope that the Problem uses same map and coordinates inside
  if (matrixType == "Laplace1D") {
    numDimensions = 1;
    map           = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
    coordinates   = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>("1D", map, galeriList);

  } else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
             matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
    numDimensions = 2;
    map           = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
    coordinates   = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>("2D", map, galeriList);

  } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
    numDimensions = 3;
    map           = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
    coordinates   = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>("3D", map, galeriList);
  }

  // Expand map to do multiple DOF per node for block problems
  if (matrixType == "Elasticity2D")
    map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 2);
  if (matrixType == "Elasticity3D")
    map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 3);

  galeriStream << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
               << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
               << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
               << "========================================================" << std::endl;

  if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
    // Our default test case for elasticity: all boundaries of a square/cube have Neumann b.c. except left which has Dirichlet
    galeriList.set("right boundary", "Neumann");
    galeriList.set("bottom boundary", "Neumann");
    galeriList.set("top boundary", "Neumann");
    galeriList.set("front boundary", "Neumann");
    galeriList.set("back boundary", "Neumann");
  }

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(galeriParameters.GetMatrixType(), map, galeriList);
  A = Pr->BuildMatrix();

  if (matrixType == "Elasticity2D") {
    numDofsPerNode = 2;
  } else if (matrixType == "Elasticity3D") {
    numDofsPerNode = 3;
  }

  A->SetFixedBlockSize(numDofsPerNode);

  X = VectorFactory::Build(map);
  B = VectorFactory::Build(map);

  // we set seed for reproducibility
  Utilities::SetRandomSeed(*comm);
  X->randomize();
  A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

  Teuchos::Array<typename STS::magnitudeType> norms(1);
  B->norm2(norms);
  B->scale(one / norms[0]);
  galeriStream << "Galeri complete.\n========================================================" << std::endl;

  out << galeriStream.str();

  comm->barrier();
  tm = Teuchos::null;

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 2 - Compute region data")));

  // Loading geometric info from galeri
  if (numDimensions == 1) {
    gNodesPerDim[0] = galeriList.get<GO>("nx");
    gNodesPerDim[1] = 1;
    gNodesPerDim[2] = 1;

    lNodesPerDim[0] = galeriList.get<LO>("lnx");
    lNodesPerDim[1] = 1;
    lNodesPerDim[2] = 1;

    procsPerDim[0] = galeriList.get<GO>("mx");
    procsPerDim[1] = 1;
    procsPerDim[2] = 1;
  } else if (numDimensions == 2) {
    gNodesPerDim[0] = galeriList.get<GO>("nx");
    gNodesPerDim[1] = galeriList.get<GO>("ny");
    gNodesPerDim[2] = 1;

    lNodesPerDim[0] = galeriList.get<LO>("lnx");
    lNodesPerDim[1] = galeriList.get<LO>("lny");
    lNodesPerDim[2] = 1;

    procsPerDim[0] = galeriList.get<GO>("mx");
    procsPerDim[1] = galeriList.get<GO>("my");
    procsPerDim[2] = 1;
  } else if (numDimensions == 3) {
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
  const GO startGID = map->getMinGlobalIndex();
  startIndices[2]   = startGID / (gNodesPerDim[1] * gNodesPerDim[0]);
  const GO rem      = startGID % (gNodesPerDim[1] * gNodesPerDim[0]);
  startIndices[1]   = rem / gNodesPerDim[0];
  startIndices[0]   = rem % gNodesPerDim[0];
  endIndices[0]     = startIndices[0] + lNodesPerDim[0] - 1;
  endIndices[1]     = startIndices[1] + lNodesPerDim[1] - 1;
  endIndices[2]     = startIndices[2] + lNodesPerDim[2] - 1;

  int leftBC = 0, rightBC = 0, frontBC = 0, backBC = 0, bottomBC = 0, topBC = 0;
  if (startIndices[0] == 0) {
    leftBC = 1;
  }
  if (startIndices[1] == 0) {
    frontBC = 1;
  }
  if (startIndices[2] == 0) {
    bottomBC = 1;
  }

  if (endIndices[0] == gNodesPerDim[0] - 1) {
    rightBC = 1;
  }
  if (endIndices[1] == gNodesPerDim[1] - 1) {
    backBC = 1;
  }
  if (endIndices[2] == gNodesPerDim[2] - 1) {
    topBC = 1;
  }

  std::cout << "p=" << myRank << " | startGID= " << startGID
            << ", startIndices: " << startIndices
            << ", endIndices: " << endIndices
            << ", gNodesPerDim: " << gNodesPerDim
            << ", BCs={" << leftBC << ", " << rightBC << ", "
            << frontBC << ", " << backBC << ", "
            << bottomBC << ", " << topBC << "}" << std::endl;

  // Rule for boundary duplication
  // For any two ranks that share an interface:
  // the lowest ranks owns the interface and the highest rank gets extra nodes

  // First we count how many nodes the region needs to send and receive
  // and allocate arrays accordingly
  LO numReceive = 0, numSend = 0;
  Teuchos::Array<GO> receiveGIDs;
  Teuchos::Array<int> receivePIDs;
  Teuchos::Array<GO> sendGIDs;
  Teuchos::Array<int> sendPIDs;
  if (numDimensions == 1) {
    if (leftBC == 0) {
      numReceive = 1;
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);

      receiveGIDs[0] = startIndices[0] - 1;
      receivePIDs[0] = myRank - 1;
    }
    if (rightBC == 0) {
      numSend = 1;
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);

      sendGIDs[0] = endIndices[0];
      sendGIDs[0] = myRank + 1;
    }
  } else if (numDimensions == 2) {
    // Received nodes
    if (frontBC == 0 && leftBC == 0) {
      numReceive = lNodesPerDim[0] + lNodesPerDim[1] + 1;
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front-left corner node
      receiveGIDs[countIDs] = startGID - gNodesPerDim[0] - 1;
      receivePIDs[countIDs] = myRank - procsPerDim[0] - 1;
      ++countIDs;
      // Receive front edge nodes
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] + i;
        receivePIDs[countIDs] = myRank - procsPerDim[0];
        ++countIDs;
      }
      // Receive left edge nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        receiveGIDs[countIDs] = startGID - 1 + j * gNodesPerDim[0];
        receivePIDs[countIDs] = myRank - 1;
        ++countIDs;
      }
    } else if (frontBC == 0) {
      numReceive = lNodesPerDim[0];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front edge nodes
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] + i;
        receivePIDs[countIDs] = myRank - procsPerDim[0];
        ++countIDs;
      }
    } else if (leftBC == 0) {
      numReceive = lNodesPerDim[1];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive left edge nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        receiveGIDs[countIDs] = startGID - 1 + j * gNodesPerDim[0];
        receivePIDs[countIDs] = myRank - 1;
        ++countIDs;
      }
    }

    // Sent nodes
    if (rightBC == 0 && backBC == 0) {
      numSend = lNodesPerDim[0] + lNodesPerDim[1] + 1;
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right edge
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        sendGIDs[countIDs] = j * gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + 1;
        ++countIDs;
      }
      // Send nodes of back edge
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        sendGIDs[countIDs] = i + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
        sendPIDs[countIDs] = myRank + procsPerDim[0];
        ++countIDs;
      }
      // Send node of back-right corner
      sendGIDs[countIDs] = startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0] + lNodesPerDim[0] - 1;
      sendPIDs[countIDs] = myRank + procsPerDim[1] + 1;
      ++countIDs;
    } else if (backBC == 0) {
      numSend = lNodesPerDim[0];

      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of back edge
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        sendGIDs[countIDs] = i + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
        sendPIDs[countIDs] = myRank + procsPerDim[0];
        ++countIDs;
      }
    } else if (rightBC == 0) {
      numSend = lNodesPerDim[1];
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right edge
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        sendGIDs[countIDs] = j * gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + 1;
        ++countIDs;
      }
    }
  } else if (numDimensions == 3) {
    // Received nodes
    if ((bottomBC == 0) && (frontBC == 0) && (leftBC == 0)) {
      numReceive = lNodesPerDim[0] * lNodesPerDim[1]                 // bottom face
                   + (lNodesPerDim[0] + 1) * lNodesPerDim[2]         // front face
                   + (lNodesPerDim[1] + 1) * (lNodesPerDim[2] + 1);  // left face
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front-left-bottom corner node
      receiveGIDs[countIDs] = startGID - gNodesPerDim[0] - 1 - gNodesPerDim[1] * gNodesPerDim[0];
      receivePIDs[countIDs] = myRank - procsPerDim[0] - 1 - procsPerDim[1] * procsPerDim[0];
      ++countIDs;
      // Receive front-bottom edge nodes
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] * gNodesPerDim[1] - gNodesPerDim[0] + i;
        receivePIDs[countIDs] = myRank - procsPerDim[0] * procsPerDim[1] - procsPerDim[0];
        ++countIDs;
      }
      // Receive left-bottom edge nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] * gNodesPerDim[1] - 1 + j * gNodesPerDim[0];
        receivePIDs[countIDs] = myRank - procsPerDim[0] * procsPerDim[1] - 1;
        ++countIDs;
      }
      // Receive bottom face nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = startGID - gNodesPerDim[0] * gNodesPerDim[1] + i + j * gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[0] * procsPerDim[1];
          ++countIDs;
        }
      }
      // Receive front-left edge nodes
      for (LO k = 0; k < lNodesPerDim[1]; ++k) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] - 1 + k * gNodesPerDim[0] * gNodesPerDim[1];
        receivePIDs[countIDs] = myRank - procsPerDim[0] - 1;
        ++countIDs;
      }
      // Receive front face nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = startGID - gNodesPerDim[0] + i + k * (gNodesPerDim[1] * gNodesPerDim[0]);
          receivePIDs[countIDs] = myRank - procsPerDim[0];
          ++countIDs;
        }
      }
      // Receive left face nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          receiveGIDs[countIDs] = startGID - 1 + j * gNodesPerDim[0] + k * (gNodesPerDim[1] * gNodesPerDim[0]);
          receivePIDs[countIDs] = myRank - 1;
          ++countIDs;
        }
      }

      // Two faces received
    } else if ((bottomBC == 0) && (frontBC == 0)) {
      numReceive = lNodesPerDim[0] * lNodesPerDim[1]           // bottom face
                   + lNodesPerDim[0] * (lNodesPerDim[2] + 1);  // front face;
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front-bottom edge nodes
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] * gNodesPerDim[1] - gNodesPerDim[0] + i;
        receivePIDs[countIDs] = myRank - procsPerDim[0] * procsPerDim[1] - procsPerDim[0];
        ++countIDs;
      }
      // Receive bottom face nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = startGID - gNodesPerDim[0] * gNodesPerDim[1] + i + j * gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[0] * procsPerDim[1];
          ++countIDs;
        }
      }
      // Receive front face nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = startGID - gNodesPerDim[0] + i + k * (gNodesPerDim[1] * gNodesPerDim[0]);
          receivePIDs[countIDs] = myRank - procsPerDim[0];
          ++countIDs;
        }
      }

    } else if ((bottomBC == 0) && (leftBC == 0)) {
      numReceive = lNodesPerDim[1] * (lNodesPerDim[0] + lNodesPerDim[2] + 1);
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive left-bottom edge nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        receiveGIDs[countIDs] = j * gNodesPerDim[0] + startGID - gNodesPerDim[1] * gNodesPerDim[0] - 1;
        receivePIDs[countIDs] = myRank - procsPerDim[1] * procsPerDim[0] - 1;
        ++countIDs;
      }
      // Receive bottom face nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = j * gNodesPerDim[0] + i + startGID - gNodesPerDim[1] * gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[1] * procsPerDim[0];
          ++countIDs;
        }
      }
      // Receive left face nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          receiveGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + startGID - 1;
          receivePIDs[countIDs] = myRank - 1;
          ++countIDs;
        }
      }

    } else if ((frontBC == 0) && (leftBC == 0)) {
      numReceive = lNodesPerDim[2] * (lNodesPerDim[1] + lNodesPerDim[0] + 1);
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front-left edge nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        receiveGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + startGID - gNodesPerDim[0] - 1;
        receivePIDs[countIDs] = myRank - procsPerDim[0] - 1;
        ++countIDs;
      }
      // Receive front face nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + i + startGID - gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[0];
          ++countIDs;
        }
      }
      // Receive left face nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          receiveGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + startGID - 1;
          receivePIDs[countIDs] = myRank - 1;
          ++countIDs;
        }
      }

      // Single face received
    } else if (bottomBC == 0) {
      numReceive = lNodesPerDim[0] * lNodesPerDim[1];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive bottom face nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = j * gNodesPerDim[0] + i + startGID - gNodesPerDim[1] * gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[1] * procsPerDim[0];
          ++countIDs;
        }
      }

    } else if (frontBC == 0) {
      numReceive = lNodesPerDim[0] * lNodesPerDim[2];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front face nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + i + startGID - gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[0];
          ++countIDs;
        }
      }

    } else if (leftBC == 0) {
      numReceive = lNodesPerDim[1] * lNodesPerDim[2];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);

      LO countIDs = 0;
      // Recive left face nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          receiveGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + startGID - 1;
          receivePIDs[countIDs] = myRank - 1;
          ++countIDs;
        }
      }
    }

    // Sent nodes
    if ((topBC == 0) && (backBC == 0) && (rightBC == 0)) {
      numSend = (lNodesPerDim[0]) * (lNodesPerDim[1]) + (lNodesPerDim[0]) * (lNodesPerDim[2]) + (lNodesPerDim[1]) * (lNodesPerDim[2]) + lNodesPerDim[0] + lNodesPerDim[1] + lNodesPerDim[2] + 1;
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          sendGIDs[countIDs] = k * (gNodesPerDim[1] * gNodesPerDim[0]) + j * gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + 1;
          ++countIDs;
        }
      }
      // Send nodes of back face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = k * (gNodesPerDim[1] * gNodesPerDim[0]) + i + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of right-back edge
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        sendGIDs[countIDs] = k * (gNodesPerDim[1] * gNodesPerDim[0]) + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + procsPerDim[0] + 1;
        ++countIDs;
      }
      // Send nodes of top face
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = j * gNodesPerDim[0] + i + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of top-right edge
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        sendGIDs[countIDs] = j * gNodesPerDim[0] + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0] + 1;
        ++countIDs;
      }
      // Send nodes of top-back edge
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        sendGIDs[countIDs] = i + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0] + (lNodesPerDim[0] - 1) * gNodesPerDim[0];
        sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0] + procsPerDim[1];
        ++countIDs;
      }
      // Send node of top-back-right corner
      sendGIDs[countIDs] = startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0] + (lNodesPerDim[0] - 1) * gNodesPerDim[0] + lNodesPerDim[0] - 1;
      sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0] + procsPerDim[1] + 1;
      ++countIDs;

    } else if ((topBC == 0) && (backBC == 0)) {
      numSend = (lNodesPerDim[0] * lNodesPerDim[2])    // back face
                + (lNodesPerDim[0] * lNodesPerDim[1])  // Top face
                + (lNodesPerDim[0]);                   // top-back edge
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of back face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = k * (gNodesPerDim[1] * gNodesPerDim[0]) + i + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of top face
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = j * gNodesPerDim[0] + i + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of top-back edge
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        sendGIDs[countIDs] = i + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0] + (lNodesPerDim[0] - 1) * gNodesPerDim[0];
        sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0] + procsPerDim[1];
        ++countIDs;
      }

    } else if ((topBC == 0) && (rightBC == 0)) {
      numSend = (lNodesPerDim[1] * lNodesPerDim[2])    // right face
                + (lNodesPerDim[0] * lNodesPerDim[1])  // Top face
                + (lNodesPerDim[1]);                   // top-right edge
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          sendGIDs[countIDs] = k * (gNodesPerDim[1] * gNodesPerDim[0]) + j * gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + 1;
          ++countIDs;
        }
      }
      // Send nodes of top face
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = j * gNodesPerDim[0] + i + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of top-right edge
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        sendGIDs[countIDs] = j * gNodesPerDim[0] + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0] + 1;
        ++countIDs;
      }

    } else if ((backBC == 0) && (rightBC == 0)) {
      numSend = lNodesPerDim[2] * (lNodesPerDim[0] + lNodesPerDim[1] + 1);
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          sendGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + 1;
          ++countIDs;
        }
      }
      // Send nodes of back face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + i + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of back-right edge
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        sendGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + procsPerDim[0] + 1;
        ++countIDs;
      }

    } else if (topBC == 0) {
      numSend = lNodesPerDim[0] * lNodesPerDim[1];
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of top face
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = j * gNodesPerDim[0] + i + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[0] * gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0];
          ++countIDs;
        }
      }

    } else if (backBC == 0) {
      numSend = lNodesPerDim[0] * lNodesPerDim[2];
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of back face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + i + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
          sendPIDs[countIDs] = myRank + procsPerDim[0];
          ++countIDs;
        }
      }

    } else if (rightBC == 0) {
      numSend = lNodesPerDim[1] * lNodesPerDim[2];
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          sendGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + 1;
          ++countIDs;
        }
      }
    }
  }

  // Second we actually fill the send and receive arrays with appropriate data
  // which will allow us to compute the region and composite maps.

  std::cout << "p=" << myRank << " | numReceive=" << numReceive
            << ", numSend=" << numSend << std::endl;
  std::cout << "p=" << myRank << " | receiveGIDs: " << receiveGIDs << std::endl;
  std::cout << "p=" << myRank << " | receivePIDs: " << receivePIDs << std::endl;
  std::cout << "p=" << myRank << " | sendGIDs: " << sendGIDs << std::endl;
  std::cout << "p=" << myRank << " | sendPIDs: " << sendPIDs << std::endl;

  comm->barrier();
  tm = Teuchos::null;

  // #ifdef HAVE_MUELU_CUDA
  //   if(profileSetup) cudaProfilerStart();
  // #endif

  //   tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 2 - MueLu Setup")));
  //   RCP<Hierarchy> H;
  //   RCP<Operator> Prec;
  //   A->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());

  //   const std::string userName = "user data";
  //   Teuchos::ParameterList& userParamList = paramList.sublist(userName);
  //   userParamList.set<int>("int numDimensions", numDimensions);
  //   userParamList.set<Teuchos::Array<LO> >("Array<LO> lNodesPerDim", lNodesPerDim);
  //   userParamList.set<RCP<RealValuedMultiVector> >("Coordinates", coordinates);
  //   H = MueLu::CreateXpetraPreconditioner(A, paramList, paramList);

  //   comm->barrier();
  //   tm = Teuchos::null;

  // #ifdef HAVE_MUELU_CUDA
  //   if(profileSolve) cudaProfilerStop();
  // #endif

  //   tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3 - LHS and RHS initialization")));
  //   X->putScalar(zero);
  //   tm = Teuchos::null;

  // #ifdef HAVE_MUELU_BELOS
  //   tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 5 - Belos Solve")));
  // #ifdef HAVE_MUELU_CUDA
  //   if(profileSolve) cudaProfilerStart();
  // #endif
  //   // Operator and Multivector type that will be used with Belos
  //   typedef MultiVector          MV;
  //   typedef Belos::OperatorT<MV> OP;

  //   // Define Operator and Preconditioner
  //   Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A)); // Turns a Xpetra::Matrix object into a Belos operator
  //   Teuchos::RCP<OP> belosPrec; // Turns a MueLu::Hierarchy object into a Belos operator
  //   H->IsPreconditioner(true);
  //   belosPrec = Teuchos::rcp(new Belos::MueLuOp <SC, LO, GO, NO>(H)); // Turns a MueLu::Hierarchy object into a Belos operator

  //   // Construct a Belos LinearProblem object
  //   RCP<Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
  //   if(solvePreconditioned) belosProblem->setRightPrec(belosPrec);

  //   bool set = belosProblem->setProblem();
  //   if (set == false) {
  //     out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
  //     return EXIT_FAILURE;
  //   }

  //   // Belos parameter list
  //   Teuchos::ParameterList belosList;
  //   belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
  //   belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
  //   belosList.set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  //   belosList.set("Output Frequency",      1);
  //   belosList.set("Output Style",          Belos::Brief);
  //   if (!scaleResidualHist)
  //     belosList.set("Implicit Residual Scaling", "None");

  //   // Create an iterative solver manager
  //   RCP< Belos::SolverManager<SC, MV, OP> > solver;
  //   solver = rcp(new Belos::BlockGmresSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));

  //   // Perform solve
  //   Belos::ReturnType retStatus = Belos::Unconverged;
  //   retStatus = solver->solve();

  //   // Get the number of iterations for this solve.
  //   out << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;
  //   // Check convergence
  //   if (retStatus != Belos::Converged)
  //     out << std::endl << "ERROR:  Belos did not converge! " << std::endl;
  //   else
  //     out << std::endl << "SUCCESS:  Belos converged!" << std::endl;
  // #ifdef HAVE_MUELU_CUDA
  //   if(profileSolve) cudaProfilerStop();
  // #endif
  // #endif //ifdef HAVE_MUELU_BELOS

  comm->barrier();
  tm                = Teuchos::null;
  globalTimeMonitor = Teuchos::null;

  RCP<ParameterList> reportParams = rcp(new ParameterList);
  const std::string filter        = "";
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
  return Automatic_Test_ETI(argc, argv);
}
