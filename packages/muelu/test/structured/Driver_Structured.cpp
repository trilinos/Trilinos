// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <vector>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_YamlParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Xpetra
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MatrixMatrix.hpp>
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

#ifdef HAVE_MUELU_AMGX
#include <MueLu_AMGXOperator.hpp>
#endif
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include "Xpetra_TpetraMultiVector.hpp"
#include <KokkosBlas1_abs.hpp>
#include <Tpetra_leftAndOrRightScaleCrsMatrix.hpp>
#include <Tpetra_computeRowAndColumnOneNorms.hpp>

#include <MueLu_CreateXpetraPreconditioner.hpp>

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib& lib, int argc, char* argv[]) {
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
  GO nx = 100, ny = 100, nz = 100;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

  std::string xmlFileName = "";
  clp.setOption("xml", &xmlFileName, "read parameters from an xml file");
  std::string yamlFileName = "";
  clp.setOption("yaml", &yamlFileName, "read parameters from a yaml file");
  int rerun = 1;
  clp.setOption("rerun", &rerun, "number of reruns required");
  int maxIts = 200;
  clp.setOption("its", &maxIts, "maximum number of solver iterations");
  double tol = 1e-12;
  clp.setOption("tol", &tol, "solver convergence tolerance");
  bool scaleResidualHist = true;
  clp.setOption("scale", "noscale", &scaleResidualHist, "scaled Krylov residual history");
  bool solvePreconditioned = true;
  clp.setOption("solve-preconditioned", "no-solve-preconditioned", &solvePreconditioned, "use MueLu preconditioner in solve");
  std::string operation = "solve";
  clp.setOption("op", &operation, "operation to perform: (matvec | setup | solver)");
  bool checkRapAgainstSymbolic = false;
  clp.setOption("check-rap-against-symbolic", "no-check-rap-against-symbolic", &checkRapAgainstSymbolic,
                "compare the final prebuilt RAP matrix against symbolic RAP");
  std::string rapGraphComparison = "exact";
  clp.setOption("rap-graph-comparison", &rapGraphComparison,
                "RAP graph comparison policy: exact or prebuilt-superset");
  std::string belosType = "cg";
  clp.setOption("belosType", &belosType, "belos solver type: (Pseudoblock CG | Block CG | Pseudoblock GMRES | Block GMRES | ...) see BelosSolverFactory.hpp for exhaustive list of solvers");
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

  const bool requireExactRapGraph  = rapGraphComparison == "exact";
  const bool allowRapGraphSuperset = rapGraphComparison == "prebuilt-superset";
  TEUCHOS_TEST_FOR_EXCEPTION(
      !requireExactRapGraph && !allowRapGraphSuperset,
      std::runtime_error,
      "Invalid RAP graph comparison policy '" << rapGraphComparison
                                              << "'. Valid values are 'exact' and 'prebuilt-superset'.");

  // Instead of checking each time for rank, create a rank 0 stream
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out       = *fancy;
  out.setOutputToRootOnly(0);

  ParameterList paramList;
  auto inst = xpetraParameters.GetInstantiation();

  if (yamlFileName != "") {
    Teuchos::updateParametersFromYamlFileAndBroadcast(yamlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);
  } else {
    if (inst == Xpetra::COMPLEX_INT_INT)
      xmlFileName = (xmlFileName != "" ? xmlFileName : "structured-complex.xml");
    else
      xmlFileName = (xmlFileName != "" ? xmlFileName : "structured.xml");
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

  Teuchos::ParameterList& rapList =
      paramList.sublist("Factories").sublist("myRAPFact");

  const std::string rapFactory = rapList.get<std::string>("factory", "");
  const bool isStructuredRAPFactory =
      rapFactory == "StructuredRAPFactory";
  if (isStructuredRAPFactory)
    rapList.set("rap: matrix type", matrixType);

  int numDimensions = 0;
  Teuchos::Array<LO> lNodesPerDim(3);

  // Create map and coordinates
  // In the future, we hope to be able to first create a Galeri problem, and then request map and coordinates from it
  // At the moment, however, things are fragile as we hope that the Problem uses same map and coordinates inside
  if (matrixType == "Laplace1D" || matrixType == "Elasticity1D") {
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

  if (numDimensions == 1) {
    lNodesPerDim[0] = galeriList.get<LO>("lnx");
    lNodesPerDim[1] = 1;
    lNodesPerDim[2] = 1;
  } else if (numDimensions == 2) {
    lNodesPerDim[0] = galeriList.get<LO>("lnx");
    lNodesPerDim[1] = galeriList.get<LO>("lny");
    lNodesPerDim[2] = 1;
  } else if (numDimensions == 3) {
    lNodesPerDim[0] = galeriList.get<LO>("lnx");
    lNodesPerDim[1] = galeriList.get<LO>("lny");
    lNodesPerDim[2] = galeriList.get<LO>("lnz");
  }
  if (isStructuredRAPFactory) {
    rapList.set("rap: processor grid x", Teuchos::as<int>(galeriList.get<GO>("mx")));
    rapList.set("rap: processor grid y", Teuchos::as<int>(galeriList.get<GO>("my")));
    rapList.set("rap: processor grid z", Teuchos::as<int>(galeriList.get<GO>("mz")));
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

  if (matrixType == "Elasticity2D" ||
      matrixType == "Elasticity3D") {
    nullspace = Pr->BuildNullspace();
    A->SetFixedBlockSize((galeriParameters.GetMatrixType() == "Elasticity2D") ? 2 : 3);
  }

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

#ifdef HAVE_MUELU_CUDA
  if (profileSetup) cudaProfilerStart();
#endif

  // Extract the diagonal of A (for RAPShiftFactory testing)
  RCP<Vector> Mdiag = Xpetra::VectorFactory<SC, LO, GO, NO>::Build(A->getRowMap(), false);
  A->getLocalDiagCopy(*Mdiag);

  for (int run = 0; run < rerun; ++run) {
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 2 - MueLu Setup")));
    RCP<Hierarchy> H;
    RCP<Operator> Prec;
    A->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());

    const std::string userName            = "user data";
    Teuchos::ParameterList& userParamList = paramList.sublist(userName);
    userParamList.set<int>("int numDimensions", numDimensions);
    userParamList.set<Teuchos::Array<LO> >("Array<LO> lNodesPerDim", lNodesPerDim);
    userParamList.set<RCP<RealValuedMultiVector> >("Coordinates", coordinates);
    userParamList.set<double>("double cfl", 1.0);
    userParamList.set<double>("double deltaT", 1.0);
    userParamList.set("Mdiag", Mdiag);
    userParamList.set<std::string>("string matrixType", matrixType);

    H = MueLu::CreateXpetraPreconditioner(A, paramList);

    if (checkRapAgainstSymbolic) {
      ParameterList referenceParamList(paramList);
      referenceParamList.sublist("Factories")
          .sublist("myRAPFact")
          .set<bool>("rap: prebuild coarse graph", false);

      A->SetMaxEigenvalueEstimate(-one);
      RCP<Hierarchy> referenceH =
          MueLu::CreateXpetraPreconditioner(A, referenceParamList);

      RCP<Matrix> structuredAc;
      RCP<Matrix> referenceAc;
      H->GetLevel(H->GetNumLevels() - 1)->Get("A", structuredAc);
      referenceH->GetLevel(referenceH->GetNumLevels() - 1)->Get("A", referenceAc);

      TEUCHOS_TEST_FOR_EXCEPTION(
          structuredAc.is_null() || referenceAc.is_null(),
          std::runtime_error,
          "Cannot compare RAP matrices because at least one matrix is null.");

      TEUCHOS_TEST_FOR_EXCEPTION(
          !structuredAc->isFillComplete() || !referenceAc->isFillComplete(),
          std::runtime_error,
          "Both RAP matrices must be fill complete before comparison.");

      TEUCHOS_TEST_FOR_EXCEPTION(
          !structuredAc->getRowMap()->isSameAs(*referenceAc->getRowMap()),
          std::runtime_error,
          "RAP matrices have different row maps.");

      TEUCHOS_TEST_FOR_EXCEPTION(
          !structuredAc->getDomainMap()->isSameAs(*referenceAc->getDomainMap()),
          std::runtime_error,
          "RAP matrices have different domain maps.");

      TEUCHOS_TEST_FOR_EXCEPTION(
          !structuredAc->getRangeMap()->isSameAs(*referenceAc->getRangeMap()),
          std::runtime_error,
          "RAP matrices have different range maps.");

      if (requireExactRapGraph) {
        TEUCHOS_TEST_FOR_EXCEPTION(
            structuredAc->getGlobalNumEntries() != referenceAc->getGlobalNumEntries(),
            std::runtime_error,
            "RAP matrices have different global entry counts: "
                << structuredAc->getGlobalNumEntries() << " versus "
                << referenceAc->getGlobalNumEntries());
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(
            structuredAc->getGlobalNumEntries() < referenceAc->getGlobalNumEntries(),
            std::runtime_error,
            "Prebuilt RAP graph has fewer global entries than the symbolic RAP graph: "
                << structuredAc->getGlobalNumEntries() << " versus "
                << referenceAc->getGlobalNumEntries());
      }

      // Local column indices cannot be compared directly because equivalent
      // column maps may assign different local indices to the same global
      // column. Compare the sorted global column IDs in every owned row,
      // requiring equality or symbolic-graph inclusion according to the
      // selected policy.
      const RCP<const Map> structuredRowMap = structuredAc->getRowMap();
      const RCP<const Map> structuredColMap = structuredAc->getColMap();
      const RCP<const Map> referenceColMap  = referenceAc->getColMap();

      int localGraphMismatch = 0;
      std::string localMismatchReason;

      const size_t localNumRows = structuredRowMap->getLocalNumElements();
      for (size_t row = 0; row < localNumRows; ++row) {
        const LO rowLid = Teuchos::as<LO>(row);
        const GO rowGid = structuredRowMap->getGlobalElement(rowLid);

        Teuchos::ArrayView<const LO> structuredIndices;
        Teuchos::ArrayView<const SC> structuredValues;
        Teuchos::ArrayView<const LO> referenceIndices;
        Teuchos::ArrayView<const SC> referenceValues;
        structuredAc->getLocalRowView(rowLid, structuredIndices, structuredValues);
        referenceAc->getLocalRowView(rowLid, referenceIndices, referenceValues);

        std::vector<GO> structuredGids(structuredIndices.size());
        std::vector<GO> referenceGids(referenceIndices.size());
        for (int entry = 0; entry < structuredIndices.size(); ++entry)
          structuredGids[entry] = structuredColMap->getGlobalElement(structuredIndices[entry]);
        for (int entry = 0; entry < referenceIndices.size(); ++entry)
          referenceGids[entry] = referenceColMap->getGlobalElement(referenceIndices[entry]);

        std::sort(structuredGids.begin(), structuredGids.end());
        std::sort(referenceGids.begin(), referenceGids.end());

        const bool structuredHasDuplicates =
            std::adjacent_find(structuredGids.begin(), structuredGids.end()) != structuredGids.end();
        const bool referenceHasDuplicates =
            std::adjacent_find(referenceGids.begin(), referenceGids.end()) != referenceGids.end();

        const bool graphMatches =
            requireExactRapGraph
                ? structuredGids == referenceGids
                : std::includes(structuredGids.begin(), structuredGids.end(),
                                referenceGids.begin(), referenceGids.end());

        if (structuredHasDuplicates || referenceHasDuplicates || !graphMatches) {
          localGraphMismatch = 1;
          std::ostringstream reason;
          reason << "rank " << comm->getRank() << ", row GID " << rowGid
                 << ": prebuilt columns = {";
          for (size_t entry = 0; entry < structuredGids.size(); ++entry)
            reason << (entry == 0 ? "" : ", ") << structuredGids[entry];
          reason << "}, symbolic columns = {";
          for (size_t entry = 0; entry < referenceGids.size(); ++entry)
            reason << (entry == 0 ? "" : ", ") << referenceGids[entry];
          reason << "}";
          if (structuredHasDuplicates)
            reason << "; prebuilt row contains duplicate columns";
          if (referenceHasDuplicates)
            reason << "; symbolic row contains duplicate columns";
          localMismatchReason = reason.str();
          break;
        }
      }

      // Make the decision collectively before any rank throws. Otherwise a
      // rank with a graph mismatch could throw while other ranks enter the
      // collective matrix operations below and hang.
      int globalGraphMismatch = 0;
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1,
                         &localGraphMismatch, &globalGraphMismatch);
      const std::string graphMismatchMessage =
          requireExactRapGraph
              ? "Prebuilt and symbolic RAP graphs differ"
              : "Prebuilt RAP graph does not contain symbolic RAP graph";
      TEUCHOS_TEST_FOR_EXCEPTION(
          globalGraphMismatch != 0,
          std::runtime_error,
          localGraphMismatch
              ? graphMismatchMessage + ": " + localMismatchReason
              : graphMismatchMessage + " on another MPI rank.");

      RCP<Matrix> difference;
      // Compute difference matrix = structuredAc - referenceAc
      Xpetra::MatrixMatrix<SC, LO, GO, NO>::TwoMatrixAdd(
          *structuredAc, false, one,
          *referenceAc, false, -one,
          difference, out);

      if (!difference->isFillComplete())
        difference->fillComplete(structuredAc->getDomainMap(), structuredAc->getRangeMap());

      // Compute Frobenius norm of difference matrix and reference matrix, then compute relative error
      const real_type differenceNorm = difference->getFrobeniusNorm();
      const real_type referenceNorm  = referenceAc->getFrobeniusNorm();
      const real_type scale =
          std::max(referenceNorm, Teuchos::ScalarTraits<real_type>::one());
      const real_type relativeError = differenceNorm / scale;
      const real_type comparisonTolerance =
          1000.0 * Teuchos::ScalarTraits<real_type>::eps();

      TEUCHOS_TEST_FOR_EXCEPTION(
          !(relativeError <= comparisonTolerance),
          std::runtime_error,
          "Final RAP matrices differ: relative Frobenius error = "
              << relativeError << ", tolerance = " << comparisonTolerance);

      out << "RAP matrix comparison passed: global nnz = "
          << structuredAc->getGlobalNumEntries()
          << ", graph comparison = "
          << (requireExactRapGraph ? "exact" : "prebuilt superset");
      if (allowRapGraphSuperset) {
        out << ", symbolic global nnz = " << referenceAc->getGlobalNumEntries()
            << ", extra prebuilt entries = "
            << structuredAc->getGlobalNumEntries() - referenceAc->getGlobalNumEntries();
      }
      out << ", relative Frobenius error = " << relativeError << std::endl;
    }

    comm->barrier();
    tm = Teuchos::null;

#ifdef HAVE_MUELU_CUDA
    if (profileSolve) cudaProfilerStop();
#endif

    if (operation == "setup") {
      continue;
    }

    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 3 - LHS and RHS initialization")));
    X->putScalar(zero);
    tm = Teuchos::null;

#ifdef HAVE_MUELU_BELOS
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 5 - Belos Solve")));
#ifdef HAVE_MUELU_CUDA
    if (profileSolve) cudaProfilerStart();
#endif
    // Operator and Multivector type that will be used with Belos
    typedef MultiVector MV;
    typedef Belos::OperatorT<MV> OP;

    // Define Operator and Preconditioner
    Teuchos::RCP<OP> belosOp = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A));  // Turns a Xpetra::Matrix object into a Belos operator
    Teuchos::RCP<OP> belosPrec;                                                       // Turns a MueLu::Hierarchy object into a Belos operator
    H->IsPreconditioner(true);
    belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO>(H));  // Turns a MueLu::Hierarchy object into a Belos operator

    // Construct a Belos LinearProblem object
    RCP<Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
    if (solvePreconditioned) belosProblem->setRightPrec(belosPrec);

    bool set = belosProblem->setProblem();
    if (set == false) {
      out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }

    // Belos parameter list
    RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();
    belosList->set("Maximum Iterations", maxIts);  // Maximum number of iterations allowed
    belosList->set("Convergence Tolerance", tol);  // Relative convergence tolerance requested
    belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList->set("Output Frequency", 1);
    belosList->set("Output Style", Belos::Brief);
    if (!scaleResidualHist)
      belosList->set("Implicit Residual Scaling", "None");

    // Create an iterative solver manager
    Belos::SolverFactory<SC, MV, OP> solverFactory;
    RCP<Belos::SolverManager<SC, MV, OP> > solver = solverFactory.create(belosType, belosList);
    solver->setProblem(belosProblem);

    // Perform solve
    Belos::ReturnType retStatus = Belos::Unconverged;
    retStatus                   = solver->solve();

    // Get the number of iterations for this solve.
    out << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;
    // Check convergence
    if (retStatus != Belos::Converged)
      out << std::endl
          << "ERROR:  Belos did not converge! " << std::endl;
    else
      out << std::endl
          << "SUCCESS:  Belos converged!" << std::endl;
#ifdef HAVE_MUELU_CUDA
    if (profileSolve) cudaProfilerStop();
#endif
#endif  // ifdef HAVE_MUELU_BELOS

    comm->barrier();
    tm = Teuchos::null;
  }  // Rerun loop

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

int main(int argc, char* argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
