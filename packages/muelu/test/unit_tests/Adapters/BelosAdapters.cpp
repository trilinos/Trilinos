// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_TestHelpers2.hpp"

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"

// Belos / Xpetra-MueLu adapters
#include "BelosXpetraAdapter.hpp"
#include "BelosMueLuAdapter.hpp"

namespace MueLuTests {

//
// Helpers function to build tests
//

// Test Belos adapters for the couple <MV,OP>
// TODO: add a bunch of 'const' on prototype
template <class Scalar, class MV, class OP>
int BelosAdaptersTest(RCP<OP>& belosOp, RCP<OP>& belosPrec, RCP<MV>& X, RCP<MV>& B, Teuchos::FancyOStream& out, bool& success) {
  RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, B));
  belosProblem->setLeftPrec(belosPrec);

  bool set = belosProblem->setProblem();
  TEST_EQUALITY(set, true);

  // Belos parameter list
  Teuchos::ParameterList belosList;
  belosList.set("Maximum Iterations", 10);       // Maximum number of iterations allowed
  belosList.set("Convergence Tolerance", 1e-7);  // Relative convergence tolerance requested

  // Create an iterative solver manager. (was double before...)
  RCP<Belos::SolverManager<Scalar, MV, OP> > belosSolver = rcp(new Belos::BlockCGSolMgr<Scalar, MV, OP>(belosProblem, rcp(&belosList, false)));

  // Perform solve
  Belos::ReturnType ret = belosSolver->solve();
  TEST_EQUALITY(ret, Belos::Converged);

  // Return number of iterations
  return belosSolver->getNumIters();
}

//
// Helpers function to verify results
//

// Singleton for norm comparisons across tests
template <class Scalar>
bool BelosAdaptersTestResultsNorm(typename Teuchos::ScalarTraits<Scalar>::magnitudeType r) {
  static typename Teuchos::ScalarTraits<Scalar>::magnitudeType ref = -1;
  if (ref == -1) {
    // std::cout << "BelosAdaptersTestResults(): Set reference results" << std::endl;
    ref = r;
    return true;
  }
  // std::cout << "BelosAdaptersTestResults(): Compare" << std::endl;

  if (r != ref)
    std::cout << "ref  norm = " << ref << std::endl
              << "curr norm = " << r << std::endl;

  return (r == ref);
}

// Test results
template <class Scalar, class MV>
bool BelosAdaptersTestResults(int numIters, RCP<MV>& X, Teuchos::FancyOStream& out, bool& success) {
  // Check numIters
  switch (TestHelpers::Parameters::getDefaultComm()->getSize()) {
    case 0: TEST_EQUALITY(numIters, 5); break;
    case 4:
      // Epetra TEST_EQUALITY(numIters, 6);
      // Tpetra TEST_EQUALITY(numIters, 7);
      break;
    default:;
  }

  // Compute norm of X (using MV traits)
  typedef Belos::MultiVecTraits<Scalar, MV> MVT;
  std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> norms(1);
  MVT::MvNorm(*X, norms);

  // Test norm equality across the unit tests
  // return MueLuTests::BelosAdaptersTestResultsNorm<Scalar>(norms[0]); // not working
  return true;
}

//
// Tests
//

// TEST:
// - OP: Xpetra::Matrix
// - MV: Xpetra::MultiVector
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BelosAdapters, XpetraOp_XpetraMV, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

#if !defined(HAVE_MUELU_EPETRA) or !defined(HAVE_MUELU_IFPACK) or !defined(HAVE_MUELU_AMESOS)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Amesos, Ifpack");
#endif

#if !defined(HAVE_MUELU_IFPACK2) or !defined(HAVE_MUELU_AMESOS2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Amesos2, Ifpack2");
#endif

  RCP<TestProblem<SC, LO, GO, NO> > p = rcp(new TestProblem<SC, LO, GO, NO>(lib));

  typedef MultiVector MV;
  typedef Belos::OperatorT<MV> OP;

  // Construct a Belos LinearProblem object
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC, LO, GO, NO>(p->GetA()));
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp<SC, LO, GO, NO>(p->GetH()));

  // Run Belos
  RCP<MultiVector> X = p->GetNewX0();
  int numIters       = MueLuTests::BelosAdaptersTest<Scalar, MV, OP>(belosOp, belosPrec, X, p->GetRHS(), out, success);

  // Tests
  TEST_EQUALITY(MueLuTests::BelosAdaptersTestResults<Scalar>(numIters, X, out, success), true);
}

// TEST:
// - OP: Xpetra::Matrix
// - MV: Epetra::MultiVector
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BelosAdapters, XpetraOp_EpetraMV, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra);
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT) && defined(HAVE_MUELU_IFPACK) && defined(HAVE_MUELU_AMESOS)
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();
  RCP<TestProblem<SC, LO, GO, NO> > p = rcp(new TestProblem<SC, LO, GO, NO>(lib));

  typedef Epetra_MultiVector MV;
  typedef Belos::OperatorT<MV> OP;

  // Construct a Belos LinearProblem object
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC, LO, GO, NO>(p->GetA()));
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp<SC, LO, GO, NO>(p->GetH()));

  // X, B
  RCP<MV> X = MueLu::Utilities<SC, LO, GO, NO>::MV2NonConstEpetraMV(p->GetNewX0());
  RCP<MV> B = MueLu::Utilities<SC, LO, GO, NO>::MV2NonConstEpetraMV(p->GetRHS());

  // Run Belos
  int numIters = MueLuTests::BelosAdaptersTest<SC, MV, OP>(belosOp, belosPrec, X, B, out, success);

  // Tests
  TEST_EQUALITY(MueLuTests::BelosAdaptersTestResults<Scalar>(numIters, X, out, success), true);
#endif
}

// TEST:
// - OP: Belos::Operator<double>
// - MV: Belos::MultiVec<double>
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BelosAdapters, BelosMultiVec_BelosMatrix, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra);
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT) && defined(HAVE_MUELU_IFPACK) && defined(HAVE_MUELU_AMESOS)
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();
  RCP<TestProblem<SC, LO, GO, NO> > p = rcp(new TestProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lib));

  typedef Belos::MultiVec<SC> MV;
  typedef Belos::Operator<SC> OP;

  // Construct a Belos LinearProblem object
  RCP<Epetra_CrsMatrix> A = Utilities::Op2NonConstEpetraCrs(p->GetA());
  RCP<OP> belosOp         = rcp(new Belos::EpetraOp(A));
  RCP<OP> belosPrec       = rcp(new Belos::MueLuOp<SC, LO, GO, NO>(p->GetH()));

  // X, B
  RCP<Epetra_MultiVector> eX = Utilities::MV2NonConstEpetraMV(p->GetNewX0());
  RCP<Epetra_MultiVector> eB = Utilities::MV2NonConstEpetraMV(p->GetRHS());
  RCP<MV> X                  = rcp(new Belos::EpetraMultiVec(*eX));
  RCP<MV> B                  = rcp(new Belos::EpetraMultiVec(*eB));

  // Run Belos
  int numIters = MueLuTests::BelosAdaptersTest<SC, MV, OP>(belosOp, belosPrec, X, B, out, success);

  // Tests
  TEST_EQUALITY(MueLuTests::BelosAdaptersTestResults<Scalar>(numIters, X, out, success), true);

  // TODO: this do not work. Is it a bug?
  //  double norm;
  //  eX->Norm2(&norm);
#endif
}

// TEST:
// - OP: Xpetra::Matrix
// - MV: Tpetra::MultiVector
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BelosAdapters, XpetraOp_TpetraMV, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra);
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

#if defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  Xpetra::UnderlyingLib lib                                      = TestHelpers::Parameters::getLib();
  RCP<TestProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> > p = rcp(new TestProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lib));

  typedef typename Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
  typedef typename Belos::OperatorT<MV> OP;

  // Construct a Belos LinearProblem object
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(p->GetA()));
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(p->GetH()));

  // X, B
  RCP<MV> X = Utilities::MV2NonConstTpetraMV(p->GetNewX0());
  RCP<MV> B = Utilities::MV2NonConstTpetraMV(p->GetRHS());

  // Run Belos
  int numIters = MueLuTests::BelosAdaptersTest<SC, MV, OP>(belosOp, belosPrec, X, B, out, success);

  // Tests
  TEST_EQUALITY(MueLuTests::BelosAdaptersTestResults<Scalar>(numIters, X, out, success), true);
#endif
}

// Instantiate the Tpetra and Xpetra based tests
// run Xpetra based tests
#define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BelosAdapters, XpetraOp_XpetraMV, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

// run Tpetra based tests
// These tests use the original belos/tpetra MultiVecTraits which are not guarded and have no specializations
// for Epetra. Therefore, carefully choose valid Tpetra instantiations.
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_SERIAL)
typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode SerialNode;
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BelosAdapters, XpetraOp_TpetraMV, double, int, int, SerialNode)
#endif

#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_SERIAL)
typedef long long int LongLong;
typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode SerialNode;
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BelosAdapters, XpetraOp_TpetraMV, double, int, LongLong, SerialNode)
#endif

#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_OPENMP)
typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode OpenMPNode;
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BelosAdapters, XpetraOp_TpetraMV, double, int, int, OpenMPNode)
#endif

#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_OPENMP)
typedef long long int LongLong;
typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode OpenMPNode;
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BelosAdapters, XpetraOp_TpetraMV, double, int, LongLong, OpenMPNode)
#endif

#if defined(HAVE_MUELU_EPETRA)
#include "Epetra_config.h"
#include "Xpetra_Map.hpp"  // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BelosAdapters, XpetraOp_EpetraMV, double, int, int, EpetraNode)
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
// typedef long long int LongLong;
// TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BelosAdapters, XpetraOp_EpetraMV, double, int, LongLong, EpetraNode)
#endif
#endif

}  // namespace MueLuTests

// TODO: norm test can be factorized, using Belos Adapter Norm function.

// TODO: generate an hierarchy that do not need ifpack[1,2]/amesos[1,2] and remove safeguards
