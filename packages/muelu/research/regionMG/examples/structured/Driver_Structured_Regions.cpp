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
#include "SetupRegionUtilities.hpp"
#include "SetupRegionVector_def.hpp"
#include "SetupRegionMatrix_def.hpp"
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
  // const int numRanks = comm->getSize();
  const int myRank   = comm->getRank();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  using STS = Teuchos::ScalarTraits<SC>;
  SC zero = STS::zero(), one = STS::one();
  using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using real_type = typename STS::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type,LO,GO,NO>;

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  GO nx = 10, ny = 10, nz = 10;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName           = "";                  clp.setOption("xml",                   &xmlFileName,           "read parameters from an xml file");
  std::string yamlFileName          = "";                  clp.setOption("yaml",                  &yamlFileName,          "read parameters from a yaml file");
  std::string convergenceLog        = "residual_norm.txt"; clp.setOption("convergence-log",       &convergenceLog,        "file in which the convergence history of the linear solver is stored");
  int         maxIts                = 200;                 clp.setOption("its",                   &maxIts,                "maximum number of solver iterations");
  std::string smootherType          = "Jacobi";            clp.setOption("smootherType",          &smootherType,          "smoother to be used: (None | Jacobi | Gauss | Chebyshev)");
  int         smootherIts           = 2;                   clp.setOption("smootherIts",           &smootherIts,           "number of smoother iterations");
  double      smootherDamp          = 0.67;                clp.setOption("smootherDamp",          &smootherDamp,          "damping parameter for the level smoother");
  double      smootherChebyEigRatio = 2.0;                 clp.setOption("smootherChebyEigRatio", &smootherChebyEigRatio, "eigenvalue ratio max/min used to approximate the smallest eigenvalue for Chebyshev relaxation");
  double      smootherChebyBoostFactor = 1.1;              clp.setOption("smootherChebyBoostFactor", &smootherChebyBoostFactor, "boost factor for Chebyshev smoother");
  double      tol                   = 1e-12;               clp.setOption("tol",                   &tol,                   "solver convergence tolerance");
  bool        scaleResidualHist     = true;                clp.setOption("scale", "noscale",      &scaleResidualHist,     "scaled Krylov residual history");
  bool        serialRandom          = false;               clp.setOption("use-serial-random", "no-use-serial-random", &serialRandom, "generate the random vector serially and then broadcast it");
  bool        keepCoarseCoords      = false;               clp.setOption("keep-coarse-coords", "no-keep-coarse-coords", &keepCoarseCoords, "keep coordinates on coarsest level of region hierarchy");
  std::string coarseSolverType      = "direct";            clp.setOption("coarseSolverType",      &coarseSolverType,      "Type of solver for (composite) coarse level operator (smoother | direct | amg)");
  std::string unstructured          = "{}";                clp.setOption("unstructured",          &unstructured,          "List of ranks to be treated as unstructured, e.g. {0, 2, 5}");
  std::string coarseAmgXmlFile      = "";                  clp.setOption("coarseAmgXml",          &coarseAmgXmlFile,      "Read parameters for AMG as coarse level solve from this xml file.");
  std::string coarseSmootherXMLFile = "";                  clp.setOption("coarseSmootherXML",     &coarseSmootherXMLFile, "File containing the parameters to use with the coarse level smoother.");
#ifdef HAVE_MUELU_TPETRA
  std::string equilibrate = "no" ;                         clp.setOption("equilibrate",           &equilibrate,           "equilibrate the system (no | diag | 1-norm)");
#endif
#ifdef HAVE_MUELU_CUDA
  bool profileSetup = false;                               clp.setOption("cuda-profile-setup", "no-cuda-profile-setup", &profileSetup, "enable CUDA profiling for setup");
  bool profileSolve = false;                               clp.setOption("cuda-profile-solve", "no-cuda-profile-solve", &profileSolve, "enable CUDA profiling for solve");
#endif
  int  cacheSize = 0;                                      clp.setOption("cachesize",               &cacheSize,           "cache size (in KB)");
  bool useStackedTimer   = false;                          clp.setOption("stacked-timer","no-stacked-timer", &useStackedTimer, "use stacked timer");
  bool showTimerSummary = true;                            clp.setOption("show-timer-summary", "no-show-timer-summary", &showTimerSummary, "Switch on/off the timer summary at the end of the run.");
  bool useFastMatVec = true;                               clp.setOption("fastMV", "no-fastMV", &useFastMatVec, "Use the fast MatVec implementation (or not)");

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

  Array<RCP<Teuchos::ParameterList> > smootherParams(1);
  smootherParams[0] = rcp(new Teuchos::ParameterList());
  smootherParams[0]->set("smoother: type",    smootherType);
  smootherParams[0]->set("smoother: sweeps",  smootherIts);
  smootherParams[0]->set("smoother: damping", smootherDamp);
  smootherParams[0]->set("smoother: Chebyshev eigRatio", smootherChebyEigRatio);
  smootherParams[0]->set("smoother: Chebyshev boost factor", smootherChebyBoostFactor);

  bool useUnstructured = false;
  Array<LO> unstructuredRanks = Teuchos::fromStringToArray<LO>(unstructured);
  for(int idx = 0; idx < unstructuredRanks.size(); ++idx) {
    if(unstructuredRanks[idx] == myRank) {useUnstructured = true;}
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
  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
  if(useStackedTimer)
    stacked_timer = rcp(new Teuchos::StackedTimer("MueLu_Driver"));
  Teuchos::TimeMonitor::setStackedTimer(stacked_timer);
  RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: S - Global Time")));
  RCP<TimeMonitor> tm                = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1 - Build Composite Matrix")));


  RCP<Matrix> A;
  RCP<Map>    nodeMap, dofMap;
  RCP<Vector> X, B;
  RCP<MultiVector>           nullspace;
  RCP<RealValuedMultiVector> coordinates;

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
  int numDofsPerNode = 0;
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
    numDofsPerNode = 2;
  } else if (matrixType == "Elasticity3D") {
    numDofsPerNode = 3;
  } else {
    numDofsPerNode = 1;
  }
  dofMap = Xpetra::MapFactory<LO,GO,Node>::Build(nodeMap, numDofsPerNode);

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

  A->SetFixedBlockSize(numDofsPerNode);

  nullspace = Pr->BuildNullspace();

  X = VectorFactory::Build(dofMap);
  B = VectorFactory::Build(dofMap);

  if(serialRandom) {
    //Build the seed on rank zero and broadcast it.
    size_t localNumElements = 0;
    if(comm->getRank() == 0) {
      localNumElements = static_cast<size_t>(dofMap->getGlobalNumElements());
    }
    RCP<Map> serialMap = MapFactory::Build(dofMap->lib(),
                                           dofMap->getGlobalNumElements(),
                                           localNumElements,
                                           0,
                                           comm);
    RCP<Vector> Xserial = VectorFactory::Build(serialMap);
    Xserial->setSeed(251743369);
    Xserial->randomize();
    RCP<Import> randomnessImporter = ImportFactory::Build(serialMap, dofMap);
    X->doImport(*Xserial, *randomnessImporter, Xpetra::INSERT);
  } else {
    // we set seed for reproducibility
    Utilities::SetRandomSeed(*comm);
    X->randomize();
  }
  A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

  Teuchos::Array<typename STS::magnitudeType> norms(1);
  B->norm2(norms);
  B->scale(one/norms[0]);
  galeriStream << "Galeri complete.\n========================================================" << std::endl;

#ifdef MATLAB_COMPARE
  Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("Ax.mm",*B);
  Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("A.mm",*A);
  B->putScalar(zero);
  Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("rhs.mm",*B);
  Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("x.mm",*X);
#endif
  out << galeriStream.str();

  comm->barrier();
  tm = Teuchos::null;

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 2 - Compute region data")));

  // Set aggregation type for each region
  Array<std::string> aggregationRegionType(1);
  RCP<ParameterList> interfaceParams = rcp(new ParameterList());
  if(useUnstructured) {
    aggregationRegionType[0] = "uncoupled";
  } else {
    aggregationRegionType[0] = "structured";
  }

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

  const LO numLocalCompositeNodes = lNodesPerDim[0]*lNodesPerDim[1]*lNodesPerDim[2];

  // Rule for boundary duplication
  // For any two ranks that share an interface:
  // the lowest rank owns the interface and the highest rank gets extra nodes

  // First we count how many nodes the region needs to send and receive
  // and allocate arrays accordingly
  Array<int> boundaryConditions;
  int maxRegPerGID = 0;
  int numInterfaces = 0;
  LO numLocalRegionNodes = 0;
  Array<GO>  sendGIDs;
  Array<int> sendPIDs;
  Array<LO>  rNodesPerDim(3);
  Array<LO>  compositeToRegionLIDs(nodeMap->getNodeNumElements()*numDofsPerNode);
  Array<GO>  quasiRegionGIDs;
  Array<GO>  quasiRegionCoordGIDs;
  Array<GO>  interfaceGIDs;
  Array<LO>  interfaceLIDsData;

  createRegionData(numDimensions, useUnstructured, numDofsPerNode,
                   gNodesPerDim(), lNodesPerDim(), procsPerDim(), nodeMap, dofMap,
                   maxRegPerGID, numLocalRegionNodes, boundaryConditions,
                   sendGIDs, sendPIDs, numInterfaces, rNodesPerDim,
                   quasiRegionGIDs, quasiRegionCoordGIDs, compositeToRegionLIDs,
                   interfaceGIDs, interfaceLIDsData);

  const LO numSend = static_cast<LO>(sendGIDs.size());

  // std::cout << "p=" << myRank << " | numSend=" << numSend << std::endl;
            // << ", numReceive=" << numReceive << std::endl;
  // std::cout << "p=" << myRank << " | receiveGIDs: " << receiveGIDs << std::endl;
  // std::cout << "p=" << myRank << " | receivePIDs: " << receivePIDs << std::endl;
  // std::cout << "p=" << myRank << " | sendGIDs: " << sendGIDs << std::endl;
  // std::cout << "p=" << myRank << " | sendPIDs: " << sendPIDs << std::endl;

  // Second we actually fill the send and receive arrays with appropriate data
  // which will allow us to compute the region and composite maps.
  // Now we can construct a list of GIDs that corresponds to rowMapPerGrp
  Array<LO>  interfacesDimensions, interfacesLIDs;
  if(useUnstructured) {
    findInterface(numDimensions, rNodesPerDim, boundaryConditions,
                  interfacesDimensions, interfacesLIDs);

    // std::cout << "p=" << myRank << " | numLocalRegionNodes=" << numLocalRegionNodes
    //           << ", rNodesPerDim: " << rNodesPerDim << std::endl;
    // std::cout << "p=" << myRank << " | boundaryConditions: " << boundaryConditions << std::endl
    //           << "p=" << myRank << " | rNodesPerDim: " << rNodesPerDim << std::endl
    //           << "p=" << myRank << " | interfacesDimensions: " << interfacesDimensions << std::endl
    //           << "p=" << myRank << " | interfacesLIDs: " << interfacesLIDs << std::endl;
  }

  interfaceParams->set<int>       ("interfaces: number",               numInterfaces);
  interfaceParams->set<Array<LO> >("interfaces: nodes per dimensions", interfacesDimensions); // nodesPerDimensions);
  interfaceParams->set<Array<LO> >("interfaces: interface nodes",      interfacesLIDs); // interfaceLIDs);

  // std::cout << "p=" << myRank << " | compositeToRegionLIDs: " << compositeToRegionLIDs() << std::endl;
  // std::cout << "p=" << myRank << " | quasiRegionGIDs: " << quasiRegionGIDs << std::endl;
  // std::cout << "p=" << myRank << " | interfaceLIDs: " << interfaceLIDs() << std::endl;
  // std::cout << "p=" << myRank << " | quasiRegionCoordGIDs: " << quasiRegionCoordGIDs() << std::endl;

  // In our very particular case we know that a node is at most shared by 4 (8) regions in 2D (3D) problems.
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

  // sleep(1);
  // if(myRank == 0) std::cout << "regionsPerGID:" << std::endl;
  // regionsPerGID->describe(*Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)), Teuchos::VERB_EXTREME);
  // sleep(1);

  comm->barrier();
  tm = Teuchos::null;

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3 - Build Region Matrix")));

  RCP<TimeMonitor> tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.1 - Build Region Maps")));

  const int maxRegPerProc = 1;
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
                                                                 numLocalRegionNodes*numDofsPerNode,
                                                                 dofMap->getIndexBase(),
                                                                 dofMap->getComm());
  revisedColMapPerGrp[0] = revisedRowMapPerGrp[0];

  // Build objects needed to construct the region coordinates
  Teuchos::RCP<Xpetra::Map<LO,GO,NO> > quasiRegCoordMap = Xpetra::MapFactory<LO,GO,Node>::
    Build(nodeMap->lib(),
          Teuchos::OrdinalTraits<GO>::invalid(),
          quasiRegionCoordGIDs(),
          nodeMap->getIndexBase(),
          nodeMap->getComm());
  Teuchos::RCP<Xpetra::Map<LO,GO,NO> > regCoordMap = Xpetra::MapFactory<LO,GO,Node>::
    Build(nodeMap->lib(),
          Teuchos::OrdinalTraits<GO>::invalid(),
          numLocalRegionNodes,
          nodeMap->getIndexBase(),
          nodeMap->getComm());

  comm->barrier();
  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.2 - Build Region Importers")));

  // Setup importers
  std::vector<RCP<Import> > rowImportPerGrp(maxRegPerProc);
  std::vector<RCP<Import> > colImportPerGrp(maxRegPerProc);
  rowImportPerGrp[0] = ImportFactory::Build(dofMap, rowMapPerGrp[0]);
  colImportPerGrp[0] = ImportFactory::Build(dofMap, colMapPerGrp[0]);
  RCP<Import> coordImporter = ImportFactory::Build(nodeMap, quasiRegCoordMap);

  comm->barrier();
  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.3 - Import ghost GIDs")));

  Array<GO>  interfaceCompositeGIDs, interfaceRegionGIDs;
  ExtractListOfInterfaceRegionGIDs(revisedRowMapPerGrp, interfaceLIDsData, interfaceRegionGIDs);

  RCP<Xpetra::MultiVector<LO, LO, GO, NO> > regionsPerGIDWithGhosts;
  RCP<Xpetra::MultiVector<GO, LO, GO, NO> > interfaceGIDsMV;
  MakeRegionPerGIDWithGhosts(nodeMap, revisedRowMapPerGrp[0], rowImportPerGrp[0],
                             maxRegPerGID, numDofsPerNode,
                             lNodesPerDim, sendGIDs, sendPIDs, interfaceLIDsData,
                             regionsPerGIDWithGhosts, interfaceGIDsMV);

  Teuchos::ArrayRCP<LO> regionMatVecLIDs;
  RCP<Import> regionInterfaceImporter;
  SetupMatVec(interfaceGIDsMV, regionsPerGIDWithGhosts, revisedRowMapPerGrp, rowImportPerGrp,
              regionMatVecLIDs, regionInterfaceImporter);

  comm->barrier();
  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.4 - Build QuasiRegion Matrix")));

  std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > quasiRegionGrpMats(1);
  MakeQuasiregionMatrices(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A), maxRegPerProc,
                          regionsPerGIDWithGhosts, rowMapPerGrp, colMapPerGrp, rowImportPerGrp,
                          quasiRegionGrpMats);

  comm->barrier();
  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3.5 - Build Region Matrix")));

  std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionGrpMats(1);
  MakeRegionMatrices(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A), A->getRowMap(), rowMapPerGrp,
                     revisedRowMapPerGrp, revisedColMapPerGrp,
                     rowImportPerGrp, maxRegPerProc, quasiRegionGrpMats, regionGrpMats);

  comm->barrier();
  tmLocal = Teuchos::null;

  comm->barrier();
  tm = Teuchos::null;

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 4 - Build Region Hierarchy")));

  // Setting up parameters before hierarchy construction
  // These need to stay in the driver as they would be provide by an app
  Array<Array<int> > regionNodesPerDim(maxRegPerProc);
  Array<RCP<MultiVector> > regionNullspace(maxRegPerProc);
  Array<RCP<RealValuedMultiVector> > regionCoordinates(maxRegPerProc);

  // Set mesh structure data
  regionNodesPerDim[0] = rNodesPerDim;

  // create nullspace vector
  regionNullspace[0] = MultiVectorFactory::Build(rowMapPerGrp[0], nullspace->getNumVectors());
  regionNullspace[0]->doImport(*nullspace, *rowImportPerGrp[0], Xpetra::INSERT);
  regionNullspace[0]->replaceMap(revisedRowMapPerGrp[0]);

  // create region coordinates vector
  regionCoordinates[0] = Xpetra::MultiVectorFactory<real_type,LO,GO,NO>::Build(quasiRegCoordMap,
                                                                               coordinates->getNumVectors());
  regionCoordinates[0]->doImport(*coordinates, *coordImporter, Xpetra::INSERT);
  regionCoordinates[0]->replaceMap(regCoordMap);

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
  Array<Array<RCP<Vector> > > regInterfaceScalings; // regional interface scaling factors on each level
  Array<RCP<Xpetra::MultiVector<GO, LO, GO, Node> > > interfaceGIDsPerLevel(1);
  Array<RCP<Xpetra::MultiVector<LO, LO, GO, Node> > > regionsPerGIDWithGhostsPerLevel(1);
  interfaceGIDsPerLevel[0] = interfaceGIDsMV;
  regionsPerGIDWithGhostsPerLevel[0] = regionsPerGIDWithGhosts;
  Array<ArrayRCP<LO> > regionMatVecLIDsPerLevel(1);
  Array<RCP<Xpetra::Import<LO, GO, Node> > > regionInterfaceImporterPerLevel(1);
  regionMatVecLIDsPerLevel[0] = regionMatVecLIDs;
  regionInterfaceImporterPerLevel[0] = regionInterfaceImporter;
  RCP<ParameterList> coarseSolverData = rcp(new ParameterList());
  coarseSolverData->set<std::string>("coarse solver type", coarseSolverType);
  coarseSolverData->set<std::string>("amg xml file", coarseAmgXmlFile);
  coarseSolverData->set<std::string>("smoother xml file", coarseSmootherXMLFile);
  RCP<ParameterList> hierarchyData = rcp(new ParameterList());


  // Create multigrid hierarchy
  createRegionHierarchy(maxRegPerProc,
                        numDimensions,
                        regionNodesPerDim,
                        aggregationRegionType,
                        interfaceParams,
                        xmlFileName,
                        regionNullspace,
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
                        interfaceGIDsPerLevel,
                        regionsPerGIDWithGhostsPerLevel,
                        regionMatVecLIDsPerLevel,
                        regionInterfaceImporterPerLevel,
                        maxRegPerGID,
                        compositeToRegionLIDs(),
                        coarseSolverData,
                        smootherParams,
                        hierarchyData,
                        keepCoarseCoords);

  hierarchyData->print();


  comm->barrier();
  tm = Teuchos::null;

  // Extract the number of levels from the prolongator data structure
  const int numLevels = regProlong.size();

  // Set data for fast MatVec
  // for (auto levelSmootherParams : smootherParams) {
  for(LO levelIdx = 0; levelIdx < numLevels; ++levelIdx) {
    smootherParams[levelIdx]->set("Use fast MatVec", useFastMatVec);
    smootherParams[levelIdx]->set("Fast MatVec: interface LIDs",
                                  regionMatVecLIDsPerLevel[levelIdx]);
    smootherParams[levelIdx]->set("Fast MatVec: interface importer",
                                  regionInterfaceImporterPerLevel[levelIdx]);
  }

  // RCP<Teuchos::FancyOStream> fancy2 = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  // Teuchos::FancyOStream& out2 = *fancy2;
  // for(LO levelIdx = 0; levelIdx < numLevels; ++levelIdx) {
  //   out2 << "p=" << myRank << " | regionMatVecLIDs on level " << levelIdx << std::endl;
  //   regionMatVecLIDsPerLevel[levelIdx]->describe(out2, Teuchos::VERB_EXTREME);
  // }

  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 5 - Solve with V-cycle")));

  {
//    std::cout << myRank << " | Running V-cycle ..." << std::endl;

    TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels>0), "We require numLevel > 0. Probably, numLevel has not been set, yet.");

    /* We first use the non-level container variables to setup the fine grid problem.
     * This is ok since the initial setup just mimics the application and the outer
     * Krylov method.
     *
     * We switch to using the level container variables as soon as we enter the
     * recursive part of the algorithm.
     */

    // Composite residual vector
    RCP<Vector> compRes = VectorFactory::Build(dofMap, true);
    {
      A->apply(*X, *compRes, Teuchos::NO_TRANS);
      compRes->update(one, *B, -one);
    }

    // transform composite vectors to regional layout
    Array<Teuchos::RCP<Vector> > quasiRegX(maxRegPerProc);
    Array<Teuchos::RCP<Vector> > regX(maxRegPerProc);
    compositeToRegional(X, quasiRegX, regX,
                        revisedRowMapPerGrp, rowImportPerGrp);

    Array<RCP<Vector> > quasiRegB(maxRegPerProc);
    Array<RCP<Vector> > regB(maxRegPerProc);
    compositeToRegional(B, quasiRegB, regB,
                        revisedRowMapPerGrp, rowImportPerGrp);
#ifdef DUMP_LOCALX_AND_A
    FILE *fp;
    char str[80];
    sprintf(str,"theMatrix.%d",myRank);
    fp = fopen(str,"w");
    fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
    LO numNzs = 0;
    for (size_t kkk = 0; kkk < regionGrpMats[0]->getNodeNumRows(); kkk++) {
      ArrayView<const LO> AAcols;
      ArrayView<const SC> AAvals;
      regionGrpMats[0]->getLocalRowView(kkk, AAcols, AAvals);
      const int *Acols    = AAcols.getRawPtr();
      const SC  *Avals = AAvals.getRawPtr();
      numNzs += AAvals.size();
    }
    fprintf(fp, "%d %d %d\n",regionGrpMats[0]->getNodeNumRows(),regionGrpMats[0]->getNodeNumRows(),numNzs);

    for (size_t kkk = 0; kkk < regionGrpMats[0]->getNodeNumRows(); kkk++) {
      ArrayView<const LO> AAcols;
      ArrayView<const SC> AAvals;
      regionGrpMats[0]->getLocalRowView(kkk, AAcols, AAvals);
      const int *Acols    = AAcols.getRawPtr();
      const SC  *Avals = AAvals.getRawPtr();
      LO RowLeng = AAvals.size();
      for (LO kk = 0; kk < RowLeng; kk++) {
          fprintf(fp, "%d %d %22.16e\n",kkk+1,Acols[kk]+1,Avals[kk]);
      }
    }
    fclose(fp);
    sprintf(str,"theX.%d",myRank);
    fp = fopen(str,"w");
    ArrayRCP<SC> lX= regX[0]->getDataNonConst(0);
    for (size_t kkk = 0; kkk < regionGrpMats[0]->getNodeNumRows(); kkk++) fprintf(fp, "%22.16e\n",lX[kkk]);
    fclose(fp);
#endif

    //    printRegionalObject<Vector>("regB 0", regB, myRank, *fos);

    Array<RCP<Vector> > regRes(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) { // step 1
      regRes[j] = VectorFactory::Build(revisedRowMapPerGrp[j], true);
    }

    /////////////////////////////////////////////////////////////////////////
    // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
    /////////////////////////////////////////////////////////////////////////

    // Prepare output of residual norm to file
    RCP<std::ofstream> log;
    if (myRank == 0)
    {
      log = rcp(new std::ofstream(convergenceLog.c_str()));
      (*log) << "# num procs = " << dofMap->getComm()->getSize() << "\n"
             << "# iteration | res-norm (scaled=" << scaleResidualHist << ")\n"
             << "#\n";
      *log << std::setprecision(16) << std::scientific;
    }

    // Print type of residual norm to the screen
    if (scaleResidualHist)
      out << "Using scaled residual norm." << std::endl;
    else
      out << "Using unscaled residual norm." << std::endl;


    // Richardson iterations
    magnitude_type normResIni = Teuchos::ScalarTraits<magnitude_type>::zero();
    const int old_precision = std::cout.precision();
    std::cout << std::setprecision(8) << std::scientific;
    int cycle = 0;

    Array<Teuchos::RCP<Vector> > regCorrect(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) {
      regCorrect[j] = VectorFactory::Build(revisedRowMapPerGrp[j], true);
    }
    for (cycle = 0; cycle < maxIts; ++cycle)
    {
      const Scalar SC_ZERO = Teuchos::ScalarTraits<SC>::zero();
      regCorrect[0]->putScalar(SC_ZERO);
      // check for convergence
      {
        ////////////////////////////////////////////////////////////////////////
        // SWITCH BACK TO NON-LEVEL VARIABLES
        ////////////////////////////////////////////////////////////////////////
        {
          if (useFastMatVec)
          {
            computeResidual(regRes, regX, regB, regionGrpMats, *smootherParams[0]);
          }
          else
          {
            computeResidual(regRes, regX, regB, regionGrpMats,
                revisedRowMapPerGrp, rowImportPerGrp);
          }

          scaleInterfaceDOFs(regRes, regInterfaceScalings[0], true);
        }

        compRes = VectorFactory::Build(dofMap, true);
        regionalToComposite(regRes, compRes, rowImportPerGrp);

        typename Teuchos::ScalarTraits<Scalar>::magnitudeType normRes = compRes->norm2();
        if(cycle == 0) { normResIni = normRes; }

        if (scaleResidualHist)
          normRes /= normResIni;

        // Output current residual norm to screen (on proc 0 only)
        out << cycle << "\t" << normRes << std::endl;
        if (myRank == 0)
          (*log) << cycle << "\t" << normRes << "\n";

        if (normRes < tol)
          break;
        }

        /////////////////////////////////////////////////////////////////////////
        // SWITCH TO RECURSIVE STYLE --> USE LEVEL CONTAINER VARIABLES
        /////////////////////////////////////////////////////////////////////////

        //      printRegionalObject<Vector>("regB 2", regB, myRank, *fos);

        bool zeroInitGuess = true;
        scaleInterfaceDOFs(regRes, regInterfaceScalings[0], false);
        vCycle(0, numLevels,
               regCorrect, regRes, regMatrices,
               regProlong, compRowMaps, quasiRegRowMaps, regRowMaps, regRowImporters,
               regInterfaceScalings, smootherParams, zeroInitGuess, coarseSolverData, hierarchyData);

        regX[0]->update(one, *regCorrect[0], one);
    }
    out << "Number of iterations performed for this solve: " << cycle << std::endl;

    std::cout << std::setprecision(old_precision);
    std::cout.unsetf(std::ios::fixed | std::ios::scientific);
  }

  comm->barrier();
  tm = Teuchos::null;
  globalTimeMonitor = Teuchos::null;

  if (showTimerSummary)
  {
    RCP<ParameterList> reportParams = rcp(new ParameterList);
    const std::string filter = "";
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
  return Automatic_Test_ETI(argc,argv);
}
