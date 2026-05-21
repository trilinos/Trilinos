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
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_BlockedMap.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include "Galeri_XpetraMaps.hpp"
#include "Galeri_MatrixTraits.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraProblemFactory.hpp"

// MueLu
#include <MueLu_TestHelpers.hpp>
#include <MueLu_MultiPhys.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosXpetraAdapter.hpp>
#endif

#ifdef HAVE_XPETRA_THYRA
#include <Thyra_DefaultBlockedLinearOp.hpp>
#endif

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct MultiPhysTestObjects {
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> monolithicA;
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> blockedA;

  Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>> auxMatrices;
  Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>>> coords;
  Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>> nullspaces;
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
BuildCartesianMap(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                  Xpetra::UnderlyingLib lib,
                  GlobalOrdinal nx,
                  GlobalOrdinal ny) {
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  galeriList.set("mx", 1);
  galeriList.set("my", comm->getSize());

  return Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node>(
      lib, "Cartesian2D", comm, galeriList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
auto BuildLaplace2D(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                    Xpetra::UnderlyingLib lib,
                    GlobalOrdinal nx,
                    GlobalOrdinal ny,
                    const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map) {
  using Map         = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using Matrix      = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using CrsWrap     = Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using MultiVector = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  galeriList.set("mx", 1);
  galeriList.set("my", comm->getSize());

  Teuchos::RCP<Galeri::Xpetra::Problem<Map, CrsWrap, MultiVector>> problem =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsWrap, MultiVector>(
          "Laplace2D", map, galeriList);

  Teuchos::RCP<CrsWrap> A = problem->BuildMatrix();
  auto coords             = problem->BuildCoords();
  auto nullspace          = problem->BuildNullspace();
  return std::make_tuple(Teuchos::rcp_dynamic_cast<Matrix>(A, true), coords, nullspace);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BuildPointCouplingMatrix(
    const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMap,
    const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap,
    Scalar value) {
  using Matrix = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  Teuchos::RCP<Matrix> A =
      Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rangeMap, 1);

  TEUCHOS_TEST_FOR_EXCEPTION(rangeMap->getLocalNumElements() != domainMap->getLocalNumElements(),
                             std::runtime_error,
                             "rangeMap and domainMap local sizes must match");

  Teuchos::ArrayView<const GlobalOrdinal> rowGids = rangeMap->getLocalElementList();
  for (size_t i = 0; i < static_cast<size_t>(rowGids.size()); ++i) {
    GlobalOrdinal row = rowGids[i];
    LocalOrdinal lid  = rangeMap->getLocalElement(row);
    GlobalOrdinal col = domainMap->getGlobalElement(lid);

    Teuchos::Array<GlobalOrdinal> cols(1, col);
    Teuchos::Array<Scalar> vals(1, value);
    A->insertGlobalValues(row, cols(), vals());
  }

  A->fillComplete(domainMap, rangeMap);
  return A;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MultiPhysTestObjects<Scalar, LocalOrdinal, GlobalOrdinal, Node>
BuildProblem(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
             Xpetra::UnderlyingLib lib,
             GlobalOrdinal nx,
             GlobalOrdinal ny,
             Scalar alpha,
             Scalar beta) {
#ifndef HAVE_XPETRA_THYRA
  return MultiPhysTestObjects<Scalar, LocalOrdinal, GlobalOrdinal, Node>{};
#else
  using Matrix           = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using BlockedCrsMatrix = Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using BlockedMap       = Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>;
  using Map              = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using real_type        = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using CoordMV          = Xpetra::MultiVector<real_type, LocalOrdinal, GlobalOrdinal, Node>;
  using MultiVector      = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  MultiPhysTestObjects<Scalar, LocalOrdinal, GlobalOrdinal, Node> objs;

  // Build the base scalar map for field 0
  Teuchos::RCP<const Map> mapU = BuildCartesianMap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(comm, lib, nx, ny);

  // Shift field 1 GIDs for Xpetra-style blocked numbering
  Teuchos::RCP<const Map> mapP = BuildCartesianMap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(comm, lib, nx, ny);

  auto [A00, coords0, nullspace0] = BuildLaplace2D<Scalar, LocalOrdinal, GlobalOrdinal, Node>(comm, lib, nx, ny, mapU);
  auto [A11, coords1, nullspace1] = BuildLaplace2D<Scalar, LocalOrdinal, GlobalOrdinal, Node>(comm, lib, nx, ny, mapP);

  Teuchos::RCP<Matrix> A01 = BuildPointCouplingMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(mapU, mapP, alpha);
  Teuchos::RCP<Matrix> A10 = BuildPointCouplingMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(mapP, mapU, beta);

  auto asThyra = [](const Teuchos::RCP<Matrix>& subblock) {
    auto op        = toTpetra(subblock);
    auto rangeMap  = Thyra::tpetraVectorSpace<Scalar>(op->getRangeMap());
    auto domainMap = Thyra::tpetraVectorSpace<Scalar>(op->getDomainMap());
    return Thyra::tpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
        rangeMap,
        domainMap,
        op);
  };

  auto blockedOperator = Teuchos::make_rcp<Thyra::DefaultBlockedLinearOp<Scalar>>();
  const int nBlock     = 2;

  blockedOperator->beginBlockFill(nBlock, nBlock);
  {
    blockedOperator->setBlock(0, 0, asThyra(A00));
    blockedOperator->setBlock(0, 1, asThyra(A01));
    blockedOperator->setBlock(1, 0, asThyra(A10));
    blockedOperator->setBlock(1, 1, asThyra(A11));
  }
  blockedOperator->endBlockFill();

  Teuchos::RCP<BlockedCrsMatrix> blockedA =
      Teuchos::rcp(new BlockedCrsMatrix(blockedOperator, Teuchos::null));

  objs.blockedA    = blockedA;
  objs.monolithicA = blockedA->Merge();

  // MultiPhys auxiliary data: one auxiliary operator per field
  objs.auxMatrices = Teuchos::ArrayRCP<Teuchos::RCP<Matrix>>(nBlock);
  objs.coords      = Teuchos::ArrayRCP<Teuchos::RCP<CoordMV>>(nBlock);
  objs.nullspaces  = Teuchos::ArrayRCP<Teuchos::RCP<MultiVector>>(nBlock);

  objs.auxMatrices[0] = A00;
  objs.auxMatrices[1] = A11;

  objs.coords[0] = coords0;
  objs.coords[1] = coords1;

  objs.nullspaces[0] = nullspace0;
  objs.nullspaces[1] = nullspace1;

  return objs;
#endif
}

#ifdef HAVE_MUELU_BELOS
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool SolveWithMultiPhys(
    const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
    const Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>& auxMatrices,
    const Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>& nullspaces,
    const Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>>>& coords,
    Teuchos::ParameterList& comboList,
    int expectedMaxIters,
    int& itersOut) {
#include <MueLu_UseShortNames.hpp>

  typedef Belos::OperatorT<MultiVector> OP;

  Teuchos::RCP<Operator> preconditioner =
      Teuchos::rcp(new MueLu::MultiPhys<SC, LO, GO, NO>(A, auxMatrices, nullspaces, coords, int(auxMatrices.size()), comboList, true));

  Teuchos::RCP<Vector> X = VectorFactory::Build(A->getRowMap());
  Teuchos::RCP<Vector> B = VectorFactory::Build(A->getRowMap());

  {
    Utilities::SetRandomSeed(*A->getRowMap()->getComm());
    X->randomize();
    A->apply(*X, *B, Teuchos::NO_TRANS, Teuchos::ScalarTraits<SC>::one(), Teuchos::ScalarTraits<SC>::zero());

    Teuchos::Array<typename Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
    B->norm2(norms);
    B->scale(Teuchos::ScalarTraits<SC>::one() / norms[0]);
    X->putScalar(Teuchos::ScalarTraits<SC>::zero());
  }

  Teuchos::RCP<OP> belosOp =
      Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A));
  Teuchos::RCP<OP> belosPrecOp =
      Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(preconditioner));

  Teuchos::RCP<Belos::LinearProblem<SC, MultiVector, OP>> belosProblem =
      Teuchos::rcp(new Belos::LinearProblem<SC, MultiVector, OP>(belosOp, X, B));

  belosProblem->setRightPrec(belosPrecOp);
  bool set = belosProblem->setProblem();
  if (!set) return false;

  Teuchos::RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();
  belosList->set("Maximum Iterations", 100);
  belosList->set("Convergence Tolerance", 1e-6);
  belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  belosList->set("Output Frequency", 1);
  belosList->set("Output Style", Belos::Brief);
  belosList->set("Implicit Residual Scaling", "None");

  Belos::SolverFactory<SC, MultiVector, OP> solverFactory;
  Teuchos::RCP<Belos::SolverManager<SC, MultiVector, OP>> solver =
      solverFactory.create("gmres", belosList);

  solver->setProblem(belosProblem);
  Belos::ReturnType retStatus = solver->solve();

  itersOut = solver->getNumIters();
  return (retStatus == Belos::Converged && itersOut <= expectedMaxIters);
}
#endif

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>
#ifndef HAVE_XPETRA_THYRA
  std::cout << "MultiphysicsDriver requires Thyra to be enabled. Skipping test.";
  return EXIT_SUCCESS;
#else

  Teuchos::oblackholestream blackhole;
  bool success = true;
  bool verbose = true;

  try {
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

    GlobalOrdinal nx = 100;
    GlobalOrdinal ny = 100;

    double alpha_in = 0.05;
    double beta_in  = -0.05;

    int expectedIters = 10;

    std::string mueluXmlFile;
    std::string pathString = "monolithic";

    clp.setOption("nx", &nx, "Global number of grid points in x");
    clp.setOption("ny", &ny, "Global number of grid points in y");
    clp.setOption("alpha", &alpha_in, "Coupling coefficient for A01");
    clp.setOption("beta", &beta_in, "Coupling coefficient for A10");
    clp.setOption("expected-iters", &expectedIters,
                  "Maximum acceptable iteration count for declaring success");

    clp.setOption("xml", &mueluXmlFile,
                  "Required: XML file containing MueLu parameters");

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseResult = clp.parse(argc, argv);
    if (parseResult != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
      return EXIT_FAILURE;

    TEUCHOS_TEST_FOR_EXCEPTION(
        mueluXmlFile.empty(),
        std::invalid_argument,
        "A MueLu XML file is required. Please provide --xml=<file.xml>.");

    auto params = Teuchos::getParametersFromXmlFile(mueluXmlFile);

    Scalar alpha = Teuchos::as<Scalar>(alpha_in);
    Scalar beta  = Teuchos::as<Scalar>(beta_in);

    auto objs          = BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>(comm, lib, nx, ny, alpha, beta);

#ifdef HAVE_MUELU_BELOS
    int iters          = -1;
    const auto solveOK = SolveWithMultiPhys<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
        objs.monolithicA,
        objs.auxMatrices,
        objs.nullspaces,
        objs.coords,
        *params,
        expectedIters,
        iters);

    if (comm->getRank() == 0) {
      std::cout << "Selected solve path: " << pathString << std::endl;
      std::cout << "Expected maximum iteration count: " << expectedIters << std::endl;

      std::cout << "Using MueLu XML file: " << mueluXmlFile << std::endl;

      if (solveOK)
        std::cout << "SUCCESS! Belos converged in " << iters << " iterations." << std::endl;
      else
        std::cout << "FAILURE! Belos did not converge fast enough." << std::endl;
    }

    success = success && solveOK;
#else
    if (comm->getRank() == 0)
      std::cout << "Belos not enabled; only assembly / merge consistency checked." << std::endl;
#endif
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
#endif
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char* argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
