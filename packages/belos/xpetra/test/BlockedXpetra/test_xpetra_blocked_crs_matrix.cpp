// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver builds a 2x2 Xpetra::BlockedCrsMatrix using Galeri Laplace2D
// diagonal blocks and simple point-wise off-diagonal couplings.  The
// right-hand side corresponds to a randomly generated solution.  The initial
// guess is set to zero.
//
// NOTE: No preconditioner is used in this case.
//

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Tpetra_Core.hpp"

// Xpetra
#include "Xpetra_BlockedCrsMatrix.hpp"
#include "Xpetra_BlockedMap.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_Vector.hpp"
#include "Xpetra_VectorFactory.hpp"

// Galeri
#include "Galeri_MatrixTraits.hpp"
#include "Galeri_XpetraMaps.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraParameters.hpp"
#include "Galeri_XpetraProblemFactory.hpp"

#include "BelosBlockGmresSolMgr.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosXpetraAdapter.hpp"

using std::cout;
using std::endl;
using std::vector;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

// ============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
BuildCartesianMap(const RCP<const Teuchos::Comm<int>> &comm,
                  Xpetra::UnderlyingLib lib, GlobalOrdinal nx,
                  GlobalOrdinal ny) {
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  galeriList.set("mx", 1);
  galeriList.set("my", comm->getSize());

  return Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node>(
      lib, "Cartesian2D", comm, galeriList);
}

// ============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> BuildLaplace2D(
    const RCP<const Teuchos::Comm<int>> &comm, Xpetra::UnderlyingLib lib,
    GlobalOrdinal nx, GlobalOrdinal ny,
    const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &map) {
  using Map = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using Matrix = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using CrsWrap =
      Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using MultiVector =
      Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  galeriList.set("mx", 1);
  galeriList.set("my", comm->getSize());

  RCP<Galeri::Xpetra::Problem<Map, CrsWrap, MultiVector>> problem =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map,
                                   CrsWrap, MultiVector>("Laplace2D", map,
                                                         galeriList);

  RCP<CrsWrap> A = problem->BuildMatrix();
  return Teuchos::rcp_dynamic_cast<Matrix>(A, true);
}

// ============================================================================
template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> BuildShiftedMap(
    const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &baseMap,
    GlobalOrdinal gidOffset) {
  Teuchos::ArrayView<const GlobalOrdinal> gids = baseMap->getLocalElementList();
  Teuchos::Array<GlobalOrdinal> shifted(gids.size());

  for (size_t i = 0; i < static_cast<size_t>(gids.size()); ++i) {
    shifted[i] = gids[i] + gidOffset;
  }

  return Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      baseMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
      shifted(), baseMap->getIndexBase(), baseMap->getComm());
}

// ============================================================================
template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> BuildFullMap(
    const std::vector<RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>>
        &subMaps) {
  TEUCHOS_TEST_FOR_EXCEPTION(subMaps.empty(), std::runtime_error,
                             "subMaps must be nonempty");

  Teuchos::Array<GlobalOrdinal> fullGids;
  for (size_t k = 0; k < subMaps.size(); ++k) {
    Teuchos::ArrayView<const GlobalOrdinal> gids =
        subMaps[k]->getLocalElementList();
    for (size_t i = 0; i < static_cast<size_t>(gids.size()); ++i) {
      fullGids.push_back(gids[i]);
    }
  }

  return Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      subMaps[0]->lib(),
      Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), fullGids(),
      subMaps[0]->getIndexBase(), subMaps[0]->getComm());
}

// ============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BuildPointCouplingMatrix(
    const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &rangeMap,
    const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &domainMap,
    Scalar value) {
  using Matrix = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  RCP<Matrix> A =
      Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
          rangeMap, 1);

  TEUCHOS_TEST_FOR_EXCEPTION(
      rangeMap->getLocalNumElements() != domainMap->getLocalNumElements(),
      std::runtime_error, "rangeMap and domainMap local sizes must match");

  Teuchos::ArrayView<const GlobalOrdinal> rowGids =
      rangeMap->getLocalElementList();
  for (size_t i = 0; i < static_cast<size_t>(rowGids.size()); ++i) {
    const GlobalOrdinal row = rowGids[i];
    const LocalOrdinal lid = rangeMap->getLocalElement(row);
    const GlobalOrdinal col = domainMap->getGlobalElement(lid);

    Teuchos::Array<GlobalOrdinal> cols(1, col);
    Teuchos::Array<Scalar> vals(1, value);
    A->insertGlobalValues(row, cols(), vals());
  }

  A->fillComplete(domainMap, rangeMap);
  return A;
}

// ============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> ShiftMatrixGIDs(
    const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &A,
    GlobalOrdinal rowOffset, GlobalOrdinal colOffset,
    const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
        &newDomainMap,
    const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
        &newRangeMap) {
  using Matrix = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using CrsWrap =
      Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using CrsMatrix =
      Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  RCP<const CrsWrap> Awrap = Teuchos::rcp_dynamic_cast<const CrsWrap>(A, true);
  RCP<const CrsMatrix> Acrs = Awrap->getCrsMatrix();

  RCP<Matrix> Ashift =
      Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
          newRangeMap, Acrs->getLocalMaxNumRowEntries());

  const size_t maxNumEntries = Acrs->getLocalMaxNumRowEntries();
  Teuchos::Array<GlobalOrdinal> inds(maxNumEntries);
  Teuchos::Array<Scalar> vals(maxNumEntries);
  size_t numEntries = 0;

  Teuchos::ArrayView<const GlobalOrdinal> rowGids =
      Acrs->getRowMap()->getLocalElementList();
  for (size_t i = 0; i < static_cast<size_t>(Acrs->getLocalNumRows()); ++i) {
    const GlobalOrdinal oldRow = rowGids[i];
    Acrs->getGlobalRowCopy(oldRow, inds(), vals(), numEntries);

    for (size_t j = 0; j < numEntries; ++j) {
      inds[j] += colOffset;
    }

    const GlobalOrdinal newRow = oldRow + rowOffset;
    Ashift->insertGlobalValues(newRow, inds(0, numEntries),
                               vals(0, numEntries));
  }

  Ashift->fillComplete(newDomainMap, newRangeMap);
  return Ashift;
}

// ============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BuildBlockedMatrix(const RCP<const Teuchos::Comm<int>> &comm,
                   Xpetra::UnderlyingLib lib, GlobalOrdinal nx,
                   GlobalOrdinal ny, Scalar alpha, Scalar beta) {
  using Matrix = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using BlockedCrsMatrix =
      Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using BlockedMap = Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>;
  using Map = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;

  // Build the base scalar map for field 0.
  RCP<const Map> mapU =
      BuildCartesianMap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(comm, lib,
                                                                   nx, ny);
  const GlobalOrdinal N =
      static_cast<GlobalOrdinal>(mapU->getGlobalNumElements());

  // Shift field 1 GIDs for Xpetra-style blocked numbering.
  RCP<const Map> mapPAux =
      BuildCartesianMap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(comm, lib,
                                                                   nx, ny);
  RCP<const Map> mapP =
      BuildShiftedMap<LocalOrdinal, GlobalOrdinal, Node>(mapPAux, N);

  // Diagonal blocks.
  RCP<Matrix> A00 = BuildLaplace2D<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
      comm, lib, nx, ny, mapU);
  RCP<Matrix> A11Aux =
      BuildLaplace2D<Scalar, LocalOrdinal, GlobalOrdinal, Node>(comm, lib, nx,
                                                                ny, mapPAux);
  RCP<Matrix> A11 = ShiftMatrixGIDs(A11Aux, N, N, mapP, mapP);

  // Off-diagonal couplings.
  RCP<Matrix> A01 =
      BuildPointCouplingMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
          mapU, mapP, alpha);
  RCP<Matrix> A10 =
      BuildPointCouplingMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
          mapP, mapU, beta);

  // Build blocked operator.
  std::vector<RCP<const Map>> subMaps(2);
  subMaps[0] = mapU;
  subMaps[1] = mapP;

  RCP<const Map> fullMap =
      BuildFullMap<LocalOrdinal, GlobalOrdinal, Node>(subMaps);
  RCP<const BlockedMap> blockedMap =
      rcp(new BlockedMap(fullMap, subMaps, false));

  RCP<BlockedCrsMatrix> blockedA =
      rcp(new BlockedCrsMatrix(blockedMap, blockedMap, 5));

  blockedA->setMatrix(0, 0, A00);
  blockedA->setMatrix(0, 1, A01);
  blockedA->setMatrix(1, 0, A10);
  blockedA->setMatrix(1, 1, A11);
  blockedA->fillComplete();

  return blockedA;
}

// ============================================================================
int main(int argc, char *argv[]) {
  typedef Tpetra::MultiVector<>::scalar_type ST;
  typedef Tpetra::MultiVector<>::local_ordinal_type LO;
  typedef Tpetra::MultiVector<>::global_ordinal_type GO;
  typedef Tpetra::MultiVector<>::node_type Node;

  typedef Teuchos::ScalarTraits<ST> SCT;
  typedef typename SCT::magnitudeType MT;
  typedef Xpetra::MultiVector<ST, LO, GO, Node> MV;
  typedef Belos::OperatorT<MV> OP;
  typedef Belos::OperatorTraits<ST, MV, OP> OPT;
  typedef Belos::MultiVecTraits<ST, MV> MVT;

  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  bool success = false;
  bool verbose = false;

  try {
    const ST one = SCT::one();
    const ST zero = SCT::zero();

    int MyPID = 0;

    RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

    //
    // Get test parameters from command-line processor
    //
    bool proc_verbose = false;
    bool debug = false;
    int frequency = -1; // how often residuals are printed by solver
    int numrhs = 1;     // total number of right-hand sides to solve for
    int blocksize = 1;  // blocksize used by solver
    int maxiters = 200; // maximum number of iterations for solver to use
    int length = 100;   // max subspace size
    GO nx = 40;         // problem size in x
    GO ny = 40;         // problem size in y
    ST alpha = Teuchos::as<ST>(0.05);
    ST beta = Teuchos::as<ST>(-0.05);
    MT tol = 1.0e-6; // relative residual tolerance

    Teuchos::CommandLineProcessor cmdp(false, true);
    cmdp.setOption("verbose", "quiet", &verbose, "Print messages and results.");
    cmdp.setOption("debug", "nodebug", &debug, "Run debugging checks.");
    cmdp.setOption("frequency", &frequency,
                   "Solver frequency for printing residuals (#iters).");
    cmdp.setOption("tol", &tol,
                   "Relative residual tolerance used by Gmres solver.");
    cmdp.setOption("nx", &nx, "Number of grid points in x-direction.");
    cmdp.setOption("ny", &ny, "Number of grid points in y-direction.");
    cmdp.setOption("alpha", &alpha,
                   "Value for the (0,1) point-coupling block.");
    cmdp.setOption("beta", &beta, "Value for the (1,0) point-coupling block.");
    cmdp.setOption("num-rhs", &numrhs,
                   "Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters", &maxiters,
                   "Maximum number of iterations per linear system.");
    cmdp.setOption(
        "max-subspace", &length,
        "Maximum number of blocks the solver can use for the subspace.");
    cmdp.setOption("block-size", &blocksize,
                   "Block size to be used by the Gmres solver.");
    if (cmdp.parse(argc, argv) !=
        Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }

    if (debug) {
      verbose = true;
    }
    if (!verbose) {
      frequency = -1;
    }

    MyPID = rank(*comm);
    proc_verbose = (verbose && (MyPID == 0));

    if (proc_verbose) {
      std::cout << Belos::Belos_Version() << std::endl << std::endl;
    }

    const Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

    //
    // Build the 2x2 blocked matrix
    //
    RCP<Xpetra::BlockedCrsMatrix<ST, LO, GO, Node>> blockedA =
        BuildBlockedMatrix<ST, LO, GO, Node>(comm, lib, nx, ny, alpha, beta);

    // BelosXpetraAdapter works with Xpetra::Matrix, and BlockedCrsMatrix
    // derives from it.
    RCP<Xpetra::Matrix<ST, LO, GO, Node>> A =
        Teuchos::rcp_dynamic_cast<Xpetra::Matrix<ST, LO, GO, Node>>(blockedA,
                                                                    true);

    RCP<const Xpetra::Map<LO, GO, Node>> map = A->getDomainMap();

    //
    // Create initial vectors
    //
    RCP<MV> B =
        Xpetra::MultiVectorFactory<ST, LO, GO, Node>::Build(map, numrhs);
    RCP<MV> X =
        Xpetra::MultiVectorFactory<ST, LO, GO, Node>::Build(map, numrhs);

    MVT::MvRandom(*X);

    RCP<OP> belosOp = rcp(new Belos::XpetraOp<ST, LO, GO, Node>(A));

    OPT::Apply(*belosOp, *X, *B);
    MVT::MvInit(*X, zero);

    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    const auto NumGlobalElements = B->getGlobalLength();

    ParameterList belosList;
    belosList.set("Block Size", blocksize);
    belosList.set("Maximum Iterations", maxiters);
    belosList.set("Convergence Tolerance", tol);
    belosList.set("Num Blocks", length);

    int verbLevel = Belos::Errors + Belos::Warnings;
    if (debug) {
      verbLevel += Belos::Debug;
    }
    if (verbose) {
      verbLevel +=
          Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails;
    }
    belosList.set("Verbosity", verbLevel);
    if (verbose) {
      if (frequency > 0) {
        belosList.set("Output Frequency", frequency);
      }
    }

    //
    // Construct a linear problem instance.
    //
    Belos::LinearProblem<ST, MV, OP> problem(belosOp, X, B);
    problem.setLabel("Belos Xpetra BlockedCrsMatrix Solve");

    bool set = problem.setProblem();
    if (set == false) {
      if (proc_verbose) {
        std::cout << std::endl
                  << "ERROR: Belos::LinearProblem failed to set up correctly!"
                  << std::endl;
      }
      return -1;
    }

    //
    // *******************************************************************
    // *************Start the block Gmres iteration***********************
    // *******************************************************************
    //
    Belos::BlockGmresSolMgr<ST, MV, OP> solver(rcpFromRef(problem),
                                               rcpFromRef(belosList));

    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Block size used by solver: " << blocksize << std::endl;
      std::cout << "Max number of Gmres iterations: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << "nx: " << nx << ", ny: " << ny << std::endl;
      std::cout << "alpha: " << alpha << ", beta: " << beta << std::endl;
      std::cout << std::endl;
    }

    //
    // Perform solve
    //
    Belos::ReturnType ret = solver.solve();

    //
    // Compute actual residuals.
    //
    bool badRes = false;
    std::vector<MT> actual_resids(numrhs);
    std::vector<MT> rhs_norm(numrhs);
    RCP<MV> resid =
        Xpetra::MultiVectorFactory<ST, LO, GO, Node>::Build(map, numrhs);

    OPT::Apply(*belosOp, *X, *resid);
    MVT::MvAddMv(-one, *resid, one, *B, *resid);
    MVT::MvNorm(*resid, actual_resids);
    MVT::MvNorm(*B, rhs_norm);

    if (proc_verbose) {
      std::cout << "---------- Actual Residuals (normalized) ----------"
                << std::endl
                << std::endl;
    }

    for (int i = 0; i < numrhs; i++) {
      MT actRes = actual_resids[i] / rhs_norm[i];
      if (proc_verbose) {
        std::cout << "Problem " << i << " : \t" << actRes << std::endl;
      }
      if (actRes > tol) {
        badRes = true;
      }
    }

    success = (ret == Belos::Converged && !badRes);

    if (success) {
      if (proc_verbose) {
        std::cout << "\nEnd Result: TEST PASSED" << std::endl;
      }
    } else {
      if (proc_verbose) {
        std::cout << "\nEnd Result: TEST FAILED" << std::endl;
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}
