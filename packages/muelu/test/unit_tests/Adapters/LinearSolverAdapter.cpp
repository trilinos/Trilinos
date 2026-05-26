// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_ParameterList.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_TpetraMultiVector.hpp>

#include <MueLu_Hierarchy.hpp>
#include <MueLu_TpetraOperator.hpp>

#include <MueLu_Details_LinearSolverFactory.hpp>
#if defined(HAVE_MUELU_EXPLICIT_INSTANTIATION)
#include <MueLu_Details_LinearSolverFactory_def.hpp>
#endif

namespace MueLuTests {

// Integration-style coverage for MueLu::Details::LinearSolver (Tpetra): Stratimikos / Trilinos::Details path.
// Exercises setMatrix -> numeric (with/without ParameterList) -> solve, then matrix update + numeric hitting
// changedA_ + ReuseTpetraPreconditioner, and the negative reuse path (Operator that is not CrsMatrix).

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(LinearSolverAdapter, FullCycle_WithoutParameterList, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  typedef Tpetra::CrsMatrix<SC, LO, GO, NO> tpetra_crsmatrix_type;
  typedef Tpetra::Operator<SC, LO, GO, NO> tpetra_operator_type;
  typedef MueLu::Details::LinearSolver<Tpetra::MultiVector<SC, LO, GO, NO>, tpetra_operator_type,
                                       typename Teuchos::ScalarTraits<SC>::magnitudeType>
      linear_solver_type;
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
    RCP<Matrix> Op                     = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(729 * comm->getSize());
    RCP<tpetra_crsmatrix_type> tpA     = Xpetra::toTpetra(Op);
    RCP<const tpetra_operator_type> A  = Teuchos::rcp_implicit_cast<const tpetra_operator_type>(tpA);

    linear_solver_type solver;
    solver.setMatrix(A);
    solver.numeric();

    RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getRowMap(), 1);
    RCP<MultiVector> X   = MultiVectorFactory::Build(Op->getRowMap(), 1);
    RHS->setSeed(271828182);
    RHS->randomize();
    Teuchos::Array<magnitude_type> norms(1);
    RHS->norm2(norms);
    RHS->scale(1 / norms[0]);
    X->putScalar((SC)0.0);

    auto Xtp = Xpetra::toTpetra(X);
    auto Btp = Xpetra::toTpetra(RHS);
    solver.solve(*Xtp, *Btp);

    X->norm2(norms);
    out << "|| X ||_2 after solve = " << norms[0] << std::endl;
    TEST_INEQUALITY(norms[0], magnitude_type(0));
  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(LinearSolverAdapter, FullCycle_WithParameterList, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  typedef Tpetra::CrsMatrix<SC, LO, GO, NO> tpetra_crsmatrix_type;
  typedef Tpetra::Operator<SC, LO, GO, NO> tpetra_operator_type;
  typedef MueLu::Details::LinearSolver<Tpetra::MultiVector<SC, LO, GO, NO>, tpetra_operator_type,
                                       typename Teuchos::ScalarTraits<SC>::magnitudeType>
      linear_solver_type;
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
    RCP<Matrix> Op                     = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(729 * comm->getSize());
    RCP<tpetra_crsmatrix_type> tpA     = Xpetra::toTpetra(Op);
    RCP<const tpetra_operator_type> A  = Teuchos::rcp_implicit_cast<const tpetra_operator_type>(tpA);

    Teuchos::ParameterList mueluList;
    linear_solver_type solver;
    solver.setMatrix(A);
    solver.setParameters(Teuchos::rcp(new Teuchos::ParameterList(mueluList)));
    solver.numeric();

    RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getRowMap(), 1);
    RCP<MultiVector> X   = MultiVectorFactory::Build(Op->getRowMap(), 1);
    RHS->setSeed(314159265);
    RHS->randomize();
    Teuchos::Array<magnitude_type> norms(1);
    RHS->norm2(norms);
    RHS->scale(1 / norms[0]);
    X->putScalar((SC)0.0);

    auto Xtp = Xpetra::toTpetra(X);
    auto Btp = Xpetra::toTpetra(RHS);
    solver.solve(*Xtp, *Btp);

    X->norm2(norms);
    out << "|| X ||_2 after solve = " << norms[0] << std::endl;
    TEST_INEQUALITY(norms[0], magnitude_type(0));
  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}

// Second setMatrix must use a different RCP (to a different CrsMatrix) so that changedA_ is set; same RCP with
// in-place value updates does not trigger reuse in LinearSolver::setMatrix.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(LinearSolverAdapter, MatrixUpdate_ReuseTpetraPreconditioner, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  typedef Tpetra::CrsMatrix<SC, LO, GO, NO> tpetra_crsmatrix_type;
  typedef Tpetra::Operator<SC, LO, GO, NO> tpetra_operator_type;
  typedef MueLu::Details::LinearSolver<Tpetra::MultiVector<SC, LO, GO, NO>, tpetra_operator_type,
                                       typename Teuchos::ScalarTraits<SC>::magnitudeType>
      linear_solver_type;
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType magnitude_type;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
    const GO gblRows                   = static_cast<GO>(243 * comm->getSize());

    RCP<Matrix> Op1 = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(gblRows);
    RCP<Matrix> Op2 = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(gblRows);

    RCP<tpetra_crsmatrix_type> tpA1    = Xpetra::toTpetra(Op1);
    RCP<tpetra_crsmatrix_type> tpA2    = Xpetra::toTpetra(Op2);
    RCP<const tpetra_operator_type> A1 = Teuchos::rcp_implicit_cast<const tpetra_operator_type>(tpA1);
    RCP<const tpetra_operator_type> A2 = Teuchos::rcp_implicit_cast<const tpetra_operator_type>(tpA2);

    Teuchos::ParameterList mueluList;
    linear_solver_type solver;
    solver.setMatrix(A1);
    solver.setParameters(Teuchos::rcp(new Teuchos::ParameterList(mueluList)));
    solver.numeric();

    RCP<MultiVector> RHS = MultiVectorFactory::Build(Op1->getRowMap(), 1);
    RCP<MultiVector> X   = MultiVectorFactory::Build(Op1->getRowMap(), 1);
    RHS->setSeed(846930886);
    RHS->randomize();
    Teuchos::Array<magnitude_type> norms(1);
    RHS->norm2(norms);
    RHS->scale(1 / norms[0]);
    X->putScalar((SC)0.0);

    solver.solve(*Xpetra::toTpetra(X), *Xpetra::toTpetra(RHS));

    // New matrix handle → changedA_ → numeric() calls ReuseTpetraPreconditioner(helperMat, *solver_)
    solver.setMatrix(A2);
    solver.numeric();

    X->putScalar((SC)0.0);
    solver.solve(*Xpetra::toTpetra(X), *Xpetra::toTpetra(RHS));

    X->norm2(norms);
    out << "|| X ||_2 after second solve = " << norms[0] << std::endl;
    TEST_INEQUALITY(norms[0], magnitude_type(0));
  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(LinearSolverAdapter, ReuseWithNonCrsOperatorThrows, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  typedef Tpetra::CrsMatrix<SC, LO, GO, NO> tpetra_crsmatrix_type;
  typedef Tpetra::Operator<SC, LO, GO, NO> tpetra_operator_type;
  typedef MueLu::TpetraOperator<SC, LO, GO, NO> muelu_tpetra_operator_type;
  typedef MueLu::Details::LinearSolver<Tpetra::MultiVector<SC, LO, GO, NO>, tpetra_operator_type,
                                       typename Teuchos::ScalarTraits<SC>::magnitudeType>
      linear_solver_type;

  if (TestHelpers::Parameters::getLib() == Xpetra::UseTpetra) {
    RCP<const Teuchos::Comm<int>> comm   = TestHelpers::Parameters::getDefaultComm();
    RCP<Matrix> Op                       = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(243 * comm->getSize());
    RCP<tpetra_crsmatrix_type> tpA       = Xpetra::toTpetra(Op);
    RCP<const tpetra_operator_type> Acrs = Teuchos::rcp_implicit_cast<const tpetra_operator_type>(tpA);

    linear_solver_type solver;
    solver.setMatrix(Acrs);
    solver.numeric();

    RCP<Hierarchy> H = rcp(new Hierarchy());
    H->setDefaultVerbLevel(Teuchos::VERB_NONE);
    RCP<MueLu::Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_NONE);
    Finest->Set("A", Op);
    H->Setup();

    RCP<muelu_tpetra_operator_type> tOp = rcp(new muelu_tpetra_operator_type(H));
    RCP<const tpetra_operator_type> AnonCrs =
        Teuchos::rcp_implicit_cast<const tpetra_operator_type>(tOp);

    solver.setMatrix(AnonCrs);
    TEST_THROW(solver.numeric(), std::runtime_error);
  } else {
    out << "This test is enabled only for linAlgebra=Tpetra." << std::endl;
  }
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Ifpack2, Amesos2)." << std::endl;
#endif
}

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node)                                                                             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(LinearSolverAdapter, FullCycle_WithoutParameterList, Scalar, LocalOrdinal, GlobalOrdinal, Node)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(LinearSolverAdapter, FullCycle_WithParameterList, Scalar, LocalOrdinal, GlobalOrdinal, Node)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(LinearSolverAdapter, MatrixUpdate_ReuseTpetraPreconditioner, Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(LinearSolverAdapter, ReuseWithNonCrsOperatorThrows, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
