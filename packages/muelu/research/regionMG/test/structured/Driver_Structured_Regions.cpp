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
#include <sstream>
#include <unistd.h>

// Teuchos
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
#include <BelosXpetraAdapter.hpp>  // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>   // => This header defines Belos::MueLuOp
#include <BelosTpetraAdapter.hpp>  // => This header defines Belos::TpetraOp
#endif

#ifdef HAVE_MUELU_CUDA
#include "cuda_profiler_api.h"
#endif

// MueLu and Xpetra Tpetra stack
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include "Xpetra_TpetraMultiVector.hpp"
#include <KokkosBlas1_abs.hpp>
#include <Tpetra_leftAndOrRightScaleCrsMatrix.hpp>
#include <Tpetra_computeRowAndColumnOneNorms.hpp>

// Xpetra Epetra stack
#ifdef HAVE_MUELU_EPETRA
#include "Xpetra_EpetraMultiVector.hpp"
#endif

#include <MueLu_CreateXpetraPreconditioner.hpp>

#if defined(HAVE_MUELU_AMESOS2)
#include <Amesos2_config.h>
#include <Amesos2.hpp>
#endif

// Region MG headers
#include "SetupRegionUtilities.hpp"
#include "SetupRegionVector_def.hpp"
#include "SetupRegionMatrix_def.hpp"
#include "SetupRegionHierarchy_def.hpp"
#include "SolveRegionHierarchy_def.hpp"

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
  // const int numRanks = comm->getSize();
  const int myRank = comm->getRank();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  using STS = Teuchos::ScalarTraits<SC>;
  SC zero = STS::zero(), one = STS::one();
  using magnitude_type        = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using real_type             = typename STS::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;

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
  std::string solverType = "region";
  clp.setOption("solverType", &solverType, "iterative solver to be used: (region | Richardson | CG)");
  std::string convergenceLog = "residual_norm.txt";
  clp.setOption("convergence-log", &convergenceLog, "file in which the convergence history of the linear solver is stored");
  int maxIts = 200;
  clp.setOption("its", &maxIts, "maximum number of solver iterations");
  double tol = 1e-12;
  clp.setOption("tol", &tol, "solver convergence tolerance");
  bool scaleResidualHist = true;
  clp.setOption("scale", "noscale", &scaleResidualHist, "scaled Krylov residual history");
  bool serialRandom = false;
  clp.setOption("use-serial-random", "no-use-serial-random", &serialRandom, "generate the random vector serially and then broadcast it");
  std::string cycleType = "V";
  clp.setOption("cycleType", &cycleType, "{Multigrid cycle type. Possible values: V, W.");
  std::string smootherType = "Jacobi";
  clp.setOption("smootherType", &smootherType, "smoother to be used: (None | Jacobi | Gauss | Chebyshev)");
  int smootherIts = 2;
  clp.setOption("smootherIts", &smootherIts, "number of smoother iterations");
  double smootherDamp = 0.67;
  clp.setOption("smootherDamp", &smootherDamp, "damping parameter for the level smoother");
  double smootherChebyEigRatio = 2.0;
  clp.setOption("smootherChebyEigRatio", &smootherChebyEigRatio, "eigenvalue ratio max/min used to approximate the smallest eigenvalue for Chebyshev relaxation");
  double smootherChebyBoostFactor = 1.1;
  clp.setOption("smootherChebyBoostFactor", &smootherChebyBoostFactor, "boost factor for Chebyshev smoother");
  bool keepCoarseCoords = false;
  clp.setOption("keep-coarse-coords", "no-keep-coarse-coords", &keepCoarseCoords, "keep coordinates on coarsest level of region hierarchy");
  bool coarseSolverRebalance = false;
  clp.setOption("rebalance-coarse", "no-rebalance-coarse", &coarseSolverRebalance, "rebalance before AMG coarse grid solve");
  int rebalanceNumPartitions = -1;
  clp.setOption("numPartitions", &rebalanceNumPartitions, "number of partitions for rebalancing the coarse grid AMG solve");
  std::string coarseSolverType = "direct";
  clp.setOption("coarseSolverType", &coarseSolverType, "Type of solver for (composite) coarse level operator (smoother | direct | amg)");
  std::string unstructured = "{}";
  clp.setOption("unstructured", &unstructured, "List of ranks to be treated as unstructured, e.g. {0, 2, 5}");
  std::string coarseAmgXmlFile = "";
  clp.setOption("coarseAmgXml", &coarseAmgXmlFile, "Read parameters for AMG as coarse level solve from this xml file.");
  std::string coarseSmootherXMLFile = "";
  clp.setOption("coarseSmootherXML", &coarseSmootherXMLFile, "File containing the parameters to use with the coarse level smoother.");
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
  bool useStackedTimer = false;
  clp.setOption("stacked-timer", "no-stacked-timer", &useStackedTimer, "use stacked timer");
  bool showTimerSummary = true;
  clp.setOption("show-timer-summary", "no-show-timer-summary", &showTimerSummary, "Switch on/off the timer summary at the end of the run.");

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

  Array<RCP<Teuchos::ParameterList> > smootherParams(1);  // TODO: this is good, resized to numlevel
  smootherParams[0] = rcp(new Teuchos::ParameterList());
  smootherParams[0]->set("smoother: type", smootherType);
  smootherParams[0]->set("smoother: sweeps", smootherIts);
  smootherParams[0]->set("smoother: damping", smootherDamp);
  smootherParams[0]->set("smoother: Chebyshev eigRatio", smootherChebyEigRatio);
  smootherParams[0]->set("smoother: Chebyshev boost factor", smootherChebyBoostFactor);

  bool useUnstructured        = false;
  Array<LO> unstructuredRanks = Teuchos::fromStringToArray<LO>(unstructured);
  for (int idx = 0; idx < unstructuredRanks.size(); ++idx) {
    if (unstructuredRanks[idx] == myRank) {
      useUnstructured = true;
    }
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
  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
  if (useStackedTimer)
    stacked_timer = rcp(new Teuchos::StackedTimer("MueLu_Driver"));
  Teuchos::TimeMonitor::setStackedTimer(stacked_timer);
  RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: S - Global Time")));
  RCP<TimeMonitor> tm                = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1 - Build Composite Matrix")));

  RCP<Matrix> A;
  RCP<Map> nodeMap, dofMap;
  RCP<Vector> X, B;
  RCP<MultiVector> nullspace;
  RCP<RealValuedMultiVector> coordinates;

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
  int numDofsPerNode     = 0;
  Teuchos::Array<GO> procsPerDim(3);
  Teuchos::Array<GO> gNodesPerDim(3);
  Teuchos::Array<LO> lNodesPerDim(3);

  // Create map and coordinates
  // In the future, we hope to be able to first create a Galeri problem, and then request map and coordinates from it
  // At the moment, however, things are fragile as we hope that the Problem uses same map and coordinates inside
  if (matrixType == "Laplace1D") {
    numDimensions = 1;
    nodeMap       = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
    coordinates   = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>("1D", nodeMap, galeriList);

  } else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
             matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
    numDimensions = 2;
    nodeMap       = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
    coordinates   = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>("2D", nodeMap, galeriList);

  } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
    numDimensions = 3;
    nodeMap       = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
    coordinates   = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LO, GO, Map, RealValuedMultiVector>("3D", nodeMap, galeriList);
  }

  // Expand map to do multiple DOF per node for block problems
  if (matrixType == "Elasticity2D") {
    numDofsPerNode = 2;
  } else if (matrixType == "Elasticity3D") {
    numDofsPerNode = 3;
  } else {
    numDofsPerNode = 1;
  }
  dofMap = Xpetra::MapFactory<LO, GO, Node>::Build(nodeMap, numDofsPerNode);

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
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(galeriParameters.GetMatrixType(), dofMap, galeriList);
  A = Pr->BuildMatrix();

  A->SetFixedBlockSize(numDofsPerNode);

  nullspace = Pr->BuildNullspace();

  X = VectorFactory::Build(dofMap);
  B = VectorFactory::Build(dofMap);

  if (serialRandom) {
    // Build the seed on rank zero and broadcast it.
    size_t localNumElements = 0;
    if (comm->getRank() == 0) {
      localNumElements = static_cast<size_t>(dofMap->getGlobalNumElements());
    }
    RCP<Map> serialMap  = MapFactory::Build(dofMap->lib(),
                                            dofMap->getGlobalNumElements(),
                                            localNumElements,
                                            0,
                                            comm);
    RCP<Vector> Xserial = VectorFactory::Build(serialMap);
    Xserial->setSeed(251743369);
    Xserial->randomize(true);  // using xpetra's randomize. Otherwise random vector is only consistent for first 128 entries
    RCP<Import> randomnessImporter = ImportFactory::Build(serialMap, dofMap);
    X->doImport(*Xserial, *randomnessImporter, Xpetra::INSERT);
  } else {
    // we set seed for reproducibility
    Utilities::SetRandomSeed(*comm);
    X->randomize(true);  // using xpetra's randomize. Otherwise random vector is only consistent for first 128 entries
  }
  A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

  Teuchos::Array<typename STS::magnitudeType> norms(1);
  B->norm2(norms);
  B->scale(one / norms[0]);
  galeriStream << "Galeri complete.\n========================================================" << std::endl;

#ifdef MATLAB_COMPARE
  Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("Ax.mm", *B);
  Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("A.mm", *A);
  B->putScalar(zero);
  Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("rhs.mm", *B);
  Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("x.mm", *X);
#endif
  out << galeriStream.str();

  comm->barrier();
  tm = Teuchos::null;

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 2 - Compute region data")));

  // Set aggregation type for each region
  std::string aggregationRegionType;
  RCP<ParameterList> interfaceParams = rcp(new ParameterList());
  if (useUnstructured) {
    aggregationRegionType = "uncoupled";
  } else {
    aggregationRegionType = "structured";
  }

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

  // Rule for boundary duplication
  // For any two ranks that share an interface:
  // the lowest rank owns the interface and the highest rank gets extra nodes

  // 1D example of the relation between Composite, Quasi Region, and Region formats
  //
  // Composite:
  // Rank 0   Rank 1
  // [0 1 2]  [3 4]
  //
  // Quasi Region:
  // Rank 0   Rank 1
  // [0 1 2]  [2 3 4]
  //
  // Region:
  // Rank 0   Rank 1
  // [0 1 2]  [5 3 4]

  // First we count how many nodes the region needs to send and receive
  // and allocate arrays accordingly
  Array<int> boundaryConditions;
  int maxRegPerGID       = 0;
  int numInterfaces      = 0;
  LO numLocalRegionNodes = 0;
  Array<GO> sendGIDs;
  Array<int> sendPIDs;
  Array<LO> rNodesPerDim(3);
  Array<LO> compositeToRegionLIDs(nodeMap->getLocalNumElements() * numDofsPerNode);
  Array<GO> quasiRegionGIDs;
  Array<GO> quasiRegionCoordGIDs;
  Array<GO> interfaceGIDs;
  Array<LO> interfaceLIDsData;

  createRegionData(numDimensions, useUnstructured, numDofsPerNode,
                   gNodesPerDim(), lNodesPerDim(), procsPerDim(), nodeMap, dofMap,
                   maxRegPerGID, numLocalRegionNodes, boundaryConditions,
                   sendGIDs, sendPIDs, numInterfaces, rNodesPerDim,
                   quasiRegionGIDs, quasiRegionCoordGIDs, compositeToRegionLIDs,
                   interfaceGIDs, interfaceLIDsData);

  // std::cout << "p=" << myRank << " | numSend=" << numSend << std::endl;
  // << ", numReceive=" << numReceive << std::endl;
  // std::cout << "p=" << myRank << " | receiveGIDs: " << receiveGIDs << std::endl;
  // std::cout << "p=" << myRank << " | receivePIDs: " << receivePIDs << std::endl;
  // std::cout << "p=" << myRank << " | sendGIDs: " << sendGIDs << std::endl;
  // std::cout << "p=" << myRank << " | sendPIDs: " << sendPIDs << std::endl;

  // Second we actually fill the send and receive arrays with appropriate data
  // which will allow us to compute the region and composite maps.
  // Now we can construct a list of GIDs that corresponds to rowMap
  Array<LO> interfacesDimensions, interfacesLIDs;
  if (useUnstructured) {
    findInterface(numDimensions, rNodesPerDim, boundaryConditions,
                  interfacesDimensions, interfacesLIDs);

    // std::cout << "p=" << myRank << " | numLocalRegionNodes=" << numLocalRegionNodes
    //           << ", rNodesPerDim: " << rNodesPerDim << std::endl;
    // std::cout << "p=" << myRank << " | boundaryConditions: " << boundaryConditions << std::endl
    //           << "p=" << myRank << " | rNodesPerDim: " << rNodesPerDim << std::endl
    //           << "p=" << myRank << " | interfacesDimensions: " << interfacesDimensions << std::endl
    //           << "p=" << myRank << " | interfacesLIDs: " << interfacesLIDs << std::endl;
  }

  interfaceParams->set<Array<LO> >("interfaces: nodes per dimensions", interfacesDimensions);  // nodesPerDimensions);
  interfaceParams->set<Array<LO> >("interfaces: interface nodes", interfacesLIDs);             // interfaceLIDs);

  // std::cout << "p=" << myRank << " | compositeToRegionLIDs: " << compositeToRegionLIDs << std::endl;
  // std::cout << "p=" << myRank << " | quasiRegionGIDs: " << quasiRegionGIDs << std::endl;
  // std::cout << "p=" << myRank << " | interfaceGIDs: " << interfaceGIDs << std::endl;
  // std::cout << "p=" << myRank << " | interfaceLIDsData: " << interfaceLIDsData << std::endl;
  // std::cout << "p=" << myRank << " | interfaceLIDs: " << interfaceLIDs << std::endl;
  // std::cout << "p=" << myRank << " | quasiRegionCoordGIDs: " << quasiRegionCoordGIDs() << std::endl;

  // In our very particular case we know that a node is at most shared by 4 (8) regions in 2D (3D) problems.
  // Other geometries will certainly have different constrains and a parallel reduction using MAX
  // would be appropriate.

  comm->barrier();
  tm = Teuchos::null;

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3 - Build Region Matrix")));

  RCP<TimeMonitor> tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.1 - Build Region Maps")));

  Teuchos::RCP<const Xpetra::Map<LO, GO, NO> > rowMap, colMap;
  Teuchos::RCP<const Xpetra::Map<LO, GO, NO> > revisedRowMap, revisedColMap;
  rowMap        = Xpetra::MapFactory<LO, GO, Node>::Build(dofMap->lib(),
                                                          Teuchos::OrdinalTraits<GO>::invalid(),
                                                          quasiRegionGIDs(),
                                                          dofMap->getIndexBase(),
                                                          dofMap->getComm());
  colMap        = rowMap;
  revisedRowMap = Xpetra::MapFactory<LO, GO, Node>::Build(dofMap->lib(),
                                                          Teuchos::OrdinalTraits<GO>::invalid(),
                                                          numLocalRegionNodes * numDofsPerNode,
                                                          dofMap->getIndexBase(),
                                                          dofMap->getComm());
  revisedColMap = revisedRowMap;

  // Build objects needed to construct the region coordinates
  Teuchos::RCP<Xpetra::Map<LO, GO, NO> > quasiRegCoordMap = Xpetra::MapFactory<LO, GO, Node>::
      Build(nodeMap->lib(),
            Teuchos::OrdinalTraits<GO>::invalid(),
            quasiRegionCoordGIDs(),
            nodeMap->getIndexBase(),
            nodeMap->getComm());
  Teuchos::RCP<Xpetra::Map<LO, GO, NO> > regCoordMap = Xpetra::MapFactory<LO, GO, Node>::
      Build(nodeMap->lib(),
            Teuchos::OrdinalTraits<GO>::invalid(),
            numLocalRegionNodes,
            nodeMap->getIndexBase(),
            nodeMap->getComm());

  comm->barrier();
  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.2 - Build Region Importers")));

  // Setup importers
  RCP<Import> rowImport;
  RCP<Import> colImport;
  rowImport                 = ImportFactory::Build(dofMap, rowMap);
  colImport                 = ImportFactory::Build(dofMap, colMap);
  RCP<Import> coordImporter = ImportFactory::Build(nodeMap, quasiRegCoordMap);

  comm->barrier();
  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.3 - Import ghost GIDs")));

  Array<GO> interfaceCompositeGIDs, interfaceRegionGIDs;
  ExtractListOfInterfaceRegionGIDs(revisedRowMap, interfaceLIDsData, interfaceRegionGIDs);

  RCP<Xpetra::MultiVector<LO, LO, GO, NO> > regionsPerGIDWithGhosts;
  RCP<Xpetra::MultiVector<GO, LO, GO, NO> > interfaceGIDsMV;
  MakeRegionPerGIDWithGhosts(nodeMap, revisedRowMap, rowImport,
                             maxRegPerGID, numDofsPerNode,
                             lNodesPerDim, sendGIDs, sendPIDs, interfaceLIDsData,
                             regionsPerGIDWithGhosts, interfaceGIDsMV);

  Teuchos::ArrayRCP<LO> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter;
  SetupMatVec(interfaceGIDsMV, regionsPerGIDWithGhosts, revisedRowMap, rowImport,
              regionMatVecLIDs, regionInterfaceImporter);

  comm->barrier();
  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.4 - Build QuasiRegion Matrix")));

  std::cout << "About to create quasi region matrix" << std::endl;
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > quasiRegionMats;
  MakeQuasiregionMatrices(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A),
                          regionsPerGIDWithGhosts, rowMap, colMap, rowImport,
                          quasiRegionMats, regionMatVecLIDs);
  std::cout << "Done creating quasi region matrix" << std::endl;

  comm->barrier();
  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.5 - Build Region Matrix")));

  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats;
  MakeRegionMatrices(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A), A->getRowMap(), rowMap,
                     revisedRowMap, revisedColMap,
                     rowImport, quasiRegionMats, regionMats);

  // If we don't need the composite operator on the fine level anymore, free it!
  if (solverType == "region") A = Teuchos::null;

  comm->barrier();
  tmLocal = Teuchos::null;

  comm->barrier();
  tm = Teuchos::null;

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 4 - Build Region Hierarchy")));

  // Setting up parameters before hierarchy construction
  // These need to stay in the driver as they would be provide by an app
  Array<int> regionNodesPerDim;
  RCP<MultiVector> regionNullspace;
  RCP<RealValuedMultiVector> regionCoordinates;

  // Set mesh structure data
  regionNodesPerDim = rNodesPerDim;

  // create nullspace vector
  regionNullspace = MultiVectorFactory::Build(rowMap, nullspace->getNumVectors());
  regionNullspace->doImport(*nullspace, *rowImport, Xpetra::INSERT);
  regionNullspace->replaceMap(revisedRowMap);

  // create region coordinates vector
  regionCoordinates = Xpetra::MultiVectorFactory<real_type, LO, GO, NO>::Build(quasiRegCoordMap,
                                                                               coordinates->getNumVectors());
  regionCoordinates->doImport(*coordinates, *coordImporter, Xpetra::INSERT);
  regionCoordinates->replaceMap(regCoordMap);

  using Tpetra_CrsMatrix   = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Tpetra_MultiVector = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  /* Stuff for multi-level algorithm
   *
   * To allow for multi-level schemes with more than two levels, we need to store
   * maps, matrices, vectors, and stuff like that on each level. Since we call the
   * multi-level scheme recursively, this should be reflected in the design of
   * variables.
   *
   * We use MueLu::Hierarchy and MueLu:Level to store each quantity on each level.
   */
  RCP<ParameterList> coarseSolverData = rcp(new ParameterList());
  coarseSolverData->set<std::string>("coarse solver type", coarseSolverType);
  coarseSolverData->set<bool>("coarse solver rebalance", coarseSolverRebalance);
  coarseSolverData->set<int>("coarse rebalance num partitions", rebalanceNumPartitions);
  coarseSolverData->set<std::string>("amg xml file", coarseAmgXmlFile);
  coarseSolverData->set<std::string>("smoother xml file", coarseSmootherXMLFile);
  RCP<ParameterList> hierarchyData = rcp(new ParameterList());

  // Create MueLu Hierarchy Initially...
  // Read MueLu parameter list form xml file
  RCP<ParameterList> mueluParams = Teuchos::rcp(new ParameterList());
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, mueluParams.ptr(), *dofMap->getComm());

  // Insert region-specific data into parameter list
  const std::string userName            = "user data";
  Teuchos::ParameterList &userParamList = mueluParams->sublist(userName);
  userParamList.set<int>("int numDimensions", numDimensions);
  userParamList.set<Array<LO> >("Array<LO> lNodesPerDim", regionNodesPerDim);
  userParamList.set<std::string>("string aggregationRegionType", aggregationRegionType);
  userParamList.set<Array<LO> >("Array<LO> nodeOnInterface", interfaceParams->get<Array<LO> >("interfaces: interface nodes"));
  userParamList.set<Array<LO> >("Array<LO> interfacesDimensions", interfaceParams->get<Array<LO> >("interfaces: nodes per dimensions"));
  if (Teuchos::nonnull(regionCoordinates)) {
    userParamList.set("Coordinates", regionCoordinates);
  }
  if (Teuchos::nonnull(regionNullspace)) {
    userParamList.set("Nullspace", regionNullspace);
  }

  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("CreateXpetraPreconditioner: Hierarchy")));

  // Create multigrid hierarchy part 1
  RCP<Hierarchy> regHierarchy = MueLu::CreateXpetraPreconditioner(regionMats, *mueluParams);

  {
    RCP<MueLu::Level> level = regHierarchy->GetLevel(0);
    level->Set<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >("rowImport", rowImport);
    level->Set<ArrayView<LocalOrdinal> >("compositeToRegionLIDs", compositeToRegionLIDs());
    level->Set<RCP<Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > >("interfaceGIDs", interfaceGIDsMV);
    level->Set<RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > >("regionsPerGIDWithGhosts", regionsPerGIDWithGhosts);
    level->Set<Teuchos::ArrayRCP<LocalOrdinal> >("regionMatVecLIDs", regionMatVecLIDs);
    level->Set<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >("regionInterfaceImporter", regionInterfaceImporter);
    level->print(std::cout, MueLu::Extreme);
  }

  tmLocal = Teuchos::null;

  // Create multigrid hierarchy part 2
  createRegionHierarchy(numDimensions,
                        regionNodesPerDim,
                        aggregationRegionType,
                        interfaceParams,
                        maxRegPerGID,
                        coarseSolverData,
                        smootherParams,
                        hierarchyData,
                        regHierarchy,
                        keepCoarseCoords);

  hierarchyData->print();

  comm->barrier();
  tm = Teuchos::null;

  // Extract the number of levels from the prolongator data structure
  const int numLevels = regHierarchy->GetNumLevels();

  // Set data for fast MatVec
  for (LO levelIdx = 0; levelIdx < numLevels; ++levelIdx) {
    RCP<MueLu::Level> level                                = regHierarchy->GetLevel(levelIdx);
    RCP<Xpetra::Import<LO, GO, NO> > regionInterfaceImport = level->Get<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >("regionInterfaceImporter");
    Teuchos::ArrayRCP<LO> regionMatVecLIDs1                = level->Get<Teuchos::ArrayRCP<LO> >("regionMatVecLIDs");
    smootherParams[levelIdx]->set("Fast MatVec: interface LIDs",
                                  regionMatVecLIDs1);
    smootherParams[levelIdx]->set("Fast MatVec: interface importer",
                                  regionInterfaceImport);
  }

  // RCP<Teuchos::FancyOStream> fancy2 = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  // Teuchos::FancyOStream& out2 = *fancy2;
  // for(LO levelIdx = 0; levelIdx < numLevels; ++levelIdx) {
  //   out2 << "p=" << myRank << " | regionMatVecLIDs on level " << levelIdx << std::endl;
  //   regionMatVecLIDsPerLevel[levelIdx]->describe(out2, Teuchos::VERB_EXTREME);
  // }

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 5 - Solve with V-cycle")));

#ifdef DUMP_LOCALX_AND_A
  FILE *fp;
  char str[80];
  sprintf(str, "theMatrix.%d", myRank);
  fp = fopen(str, "w");
  fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  LO numNzs = 0;
  for (size_t kkk = 0; kkk < regionMats->getLocalNumRows(); kkk++) {
    ArrayView<const LO> AAcols;
    ArrayView<const SC> AAvals;
    regionMats->getLocalRowView(kkk, AAcols, AAvals);
    const int *Acols = AAcols.getRawPtr();
    const SC *Avals  = AAvals.getRawPtr();
    numNzs += AAvals.size();
  }
  fprintf(fp, "%d %d %d\n", regionMats->getLocalNumRows(), regionMats->getLocalNumRows(), numNzs);

  for (size_t kkk = 0; kkk < regionMats->getLocalNumRows(); kkk++) {
    ArrayView<const LO> AAcols;
    ArrayView<const SC> AAvals;
    regionMats->getLocalRowView(kkk, AAcols, AAvals);
    const int *Acols = AAcols.getRawPtr();
    const SC *Avals  = AAvals.getRawPtr();
    LO RowLeng       = AAvals.size();
    for (LO kk = 0; kk < RowLeng; kk++) {
      fprintf(fp, "%d %d %22.16e\n", kkk + 1, Acols[kk] + 1, Avals[kk]);
    }
  }
  fclose(fp);
  // sprintf(str,"theX.%d",myRank);
  // fp = fopen(str,"w");
  // ArrayRCP<SC> lX= regX->getDataNonConst(0);
  // for (size_t kkk = 0; kkk < regionMats->getLocalNumRows(); kkk++) fprintf(fp, "%22.16e\n",lX[kkk]);
  // fclose(fp);
#endif

  if (solverType == "region") {
    solveRegionProblemRichardson(tol, scaleResidualHist, maxIts,
                                 cycleType, convergenceLog,
                                 coarseSolverData, smootherParams, hierarchyData,
                                 regHierarchy, X, B);
  } else if (solverType == "Richardson") {
    solveCompositeProblemRichardson(tol, scaleResidualHist, maxIts,
                                    cycleType, convergenceLog,
                                    coarseSolverData, smootherParams, hierarchyData,
                                    regHierarchy, A, X, B);
  } else if (solverType == "CG") {
    solveCompositeProblemPCG(tol, scaleResidualHist, maxIts,
                             cycleType, convergenceLog,
                             coarseSolverData, smootherParams, hierarchyData,
                             regHierarchy, A, X, B);
  } else {
    throw std::runtime_error("Unknown solverType: " + solverType);
  }

  comm->barrier();
  tm                = Teuchos::null;
  globalTimeMonitor = Teuchos::null;

  if (showTimerSummary) {
    RCP<ParameterList> reportParams = rcp(new ParameterList);
    const std::string filter        = "";
    if (useStackedTimer) {
      Teuchos::StackedTimer::OutputOptions options;
      options.output_fraction = options.output_histogram = options.output_minmax = true;
      stacked_timer->report(out, comm, options);
    } else {
      std::ios_base::fmtflags ff(out.flags());
      TimeMonitor::report(comm.ptr(), out, filter, reportParams);
      out << std::setiosflags(ff);
    }
  }

  TimeMonitor::clearCounters();

  return EXIT_SUCCESS;
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
