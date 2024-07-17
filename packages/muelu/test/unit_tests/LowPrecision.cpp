// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <list>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Xpetra_IO.hpp>
#include <MueLu_config.hpp>
#include <MueLu_ConfigDefs.hpp>
#include <MueLu_Version.hpp>
// #include <MueLu_NoFactory.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_LowPrecisionFactory.hpp>

#include <Tpetra_CrsMatrixMultiplyOp.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(LowPrecisionFactory, Basic, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  using Teuchos::TimeMonitor;
  using TST        = Teuchos::ScalarTraits<Scalar>;
  using HalfScalar = typename TST::halfPrecision;
  using HalfTST    = Teuchos::ScalarTraits<HalfScalar>;

  // Modify to read in matrix
  std::string matrix_file = "";

  // for timing
  int numMatvecs = 100;

  RCP<Matrix> A;
  if (matrix_file == "") {
    int nx = 1000;
    A      = TestHelpers::TestFactory<SC, LO, GO, NO>::Build2DPoisson(nx);
  } else {
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    A                                   = Xpetra::IO<SC, LO, GO, NO>::Read(matrix_file, Xpetra::UseTpetra, comm);
  }

  RCP<Operator> lowA;
  const bool isTpetra               = A->getDomainMap()->lib() == Xpetra::UseTpetra;
  std::list<std::string> matrixKeys = {"A", "P", "R"};
  for (auto it = matrixKeys.begin(); it != matrixKeys.end(); ++it) {
    Level aLevel;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(aLevel);
    aLevel.Set(*it, A);

    // Set up LowPrecisionFactory
    RCP<LowPrecisionFactory> LP = rcp(new LowPrecisionFactory());
    Teuchos::ParameterList params;
    params.set("matrix key", *it);
    LP->SetParameterList(params);

    // Build
    aLevel.Request(*it, LP.get());
    LP->Build(aLevel);
    TEST_EQUALITY(aLevel.IsAvailable(*it, LP.get()), true);
    lowA = aLevel.Get<RCP<Operator> >(*it, LP.get());
    aLevel.Release(*it, LP.get());

    // Check that the resulting Op is as expected.
    if (isTpetra && (false
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
                     || std::is_same<Scalar, double>::value
#endif
#if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
                     || std::is_same<Scalar, std::complex<double> >::value
#endif
                     )) {
      auto tpCrsMultOp = rcp_dynamic_cast<Tpetra::CrsMatrixMultiplyOp<Scalar, HalfScalar, LocalOrdinal, GlobalOrdinal, Node> >(rcp_dynamic_cast<TpetraOperator>(lowA)->getOperator());
      TEST_ASSERT(!tpCrsMultOp.is_null());  // Actually converted
    } else {
      TEST_ASSERT(!rcp_dynamic_cast<Matrix>(lowA).is_null());  // Just a regular old matrix
      TEST_EQUALITY(A, lowA);
    }
  }

  // Check apply
  if (isTpetra && std::is_same<Scalar, double>::value) {
    RCP<MultiVector> X = MultiVectorFactory::Build(A->getDomainMap(), 1);
    X->randomize();

    RCP<MultiVector> B  = MultiVectorFactory::Build(A->getRangeMap(), 1);
    RCP<MultiVector> B0 = MultiVectorFactory::Build(lowA->getRangeMap(), 1);

    // warm up
    for (int i = 0; i < 5; i++) {
      A->apply(*X, *B);
    }
    {
      RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MatVec original")));
      for (int i = 0; i < numMatvecs; i++) {
        A->apply(*X, *B);
      }
    }

    // warm up low precision
    for (int i = 0; i < 5; i++) {
      lowA->apply(*X, *B0);
    }
    {
      RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MatVec low precision")));
      for (int i = 0; i < numMatvecs; i++) {
        lowA->apply(*X, *B0);
      }
    }

    TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union, "", true);
    TimeMonitor::zeroOutTimers();

    B->update(-1, *B0, 1);

    Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> norm(1);
    B->norm2(norm);
    TEST_FLOATING_EQUALITY(norm[0], TST::zero(), HalfTST::eps());
  }
}

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(LowPrecisionFactory, Basic, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
