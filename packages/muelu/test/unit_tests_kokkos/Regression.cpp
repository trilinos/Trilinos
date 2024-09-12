// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

#include <MueLu_config.hpp>
#include <MueLu_TestHelpers_kokkos.hpp>
#include <MueLu_Version.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_CoalesceDropFactory_kokkos.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>

#include <Tpetra_Details_KokkosCounter.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Regression, H2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef typename Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers_kokkos::Parameters::getDefaultComm();
  if (comm->getSize() != 1) {
    out << "Skipping test for more than 1 proc" << std::endl;
    return;
  }
  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer = rcp(new Teuchos::StackedTimer("H2D"));
  Teuchos::TimeMonitor::setStackedTimer(stacked_timer);
  int numRows                            = 399;
  RCP<Matrix> A                          = MueLuTests::TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(numRows);
  RCP<const Map> map                     = A->getMap();
  RCP<RealValuedMultiVector> coordinates = Xpetra::MultiVectorFactory<real_type, LO, GO, NO>::Build(map, 1);
  RCP<MultiVector> nullspace             = MultiVectorFactory::Build(map, 1);
  nullspace->putScalar(Teuchos::ScalarTraits<SC>::one());

  Teuchos::ParameterList MueLuList;
  MueLuList.set("verbosity", "high");
  MueLuList.set("coarse: max size", numRows - 1);  // make it so we want two levels
  MueLuList.set("smoother: type", "RELAXATION");
  MueLuList.sublist("smoother: params").set("relaxation: type", "Jacobi");
  MueLuList.set("max levels", 2);
  MueLuList.set("use kokkos refactor", true);

  Teuchos::ParameterList userParamList = MueLuList.sublist("user data");
  userParamList.set<RCP<RealValuedMultiVector> >("Coordinates", coordinates);
  userParamList.set<RCP<MultiVector> >("Nullspace", nullspace);

  Tpetra::Details::DeepCopyCounter::reset();
  Tpetra::Details::DeepCopyCounter::start();

  Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > H =
      MueLu::CreateXpetraPreconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>(A, MueLuList);

  Tpetra::Details::DeepCopyCounter::stop();

  // confirm that we did get a hierarchy with two levels -- a sanity check for this test
  TEST_EQUALITY(2, H->GetGlobalNumLevels());

  // When Kokkos Kernels uses TPLs, some Kokkos::deep_copy in the Kokkos Kernels native implementations are not called.
  int kkNativeDeepCopies = 8;
#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
  if constexpr (Node::is_cpu)
    kkNativeDeepCopies = 0;
#endif
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
  if constexpr (std::is_same_v<typename Node::execution_space, Kokkos::Cuda>)
    kkNativeDeepCopies = 0;
#endif
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)
  if constexpr (std::is_same_v<typename Node::execution_space, Kokkos::HIP>)
    kkNativeDeepCopies = 0;
#endif

  if (Node::is_cpu) {
    TEST_EQUALITY(Tpetra::Details::DeepCopyCounter::get_count_different_space(), 0);
  }
#ifdef KOKKOS_HAS_SHARED_SPACE
  else {
    size_t targetNumDeepCopies = kkNativeDeepCopies + (std::is_same_v<typename Node::memory_space, Kokkos::SharedSpace> ? 20 : 35);
    TEST_EQUALITY(Tpetra::Details::DeepCopyCounter::get_count_different_space(), targetNumDeepCopies);
  }
#else
  else {
    TEST_EQUALITY(Tpetra::Details::DeepCopyCounter::get_count_different_space(), kkNativeDeepCopies + 35);
  }
#endif  // KOKKOS_HAS_SHARED_SPACE

  auto X = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(map, 1);
  auto B = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(map, 1);

  Tpetra::Details::DeepCopyCounter::reset();
  Tpetra::Details::DeepCopyCounter::start();

  H->Iterate(*B, *X, 1, false);

  Tpetra::Details::DeepCopyCounter::stop();

  if (Node::is_cpu) {
    TEST_EQUALITY(Tpetra::Details::DeepCopyCounter::get_count_different_space(), 0);
  } else {
    TEST_EQUALITY(Tpetra::Details::DeepCopyCounter::get_count_different_space(), 2);
  }

  stacked_timer->stopBaseTimer();
  Teuchos::StackedTimer::OutputOptions options;
  options.output_fraction = options.output_histogram = options.output_minmax = true;
  stacked_timer->report(out, comm, options);

}  // H2D

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Regression, Aggregration, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<const Teuchos::Comm<int> > comm               = TestHelpers_kokkos::Parameters::getDefaultComm();
  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer = rcp(new Teuchos::StackedTimer("H2D"));
  Teuchos::TimeMonitor::setStackedTimer(stacked_timer);

  // Make a Matrix with multiple degrees of freedom per node
  GlobalOrdinal nx = 20, ny = 20;

  // Describes the initial layout of matrix rows across processors.
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);

  using test_factory = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>;
  Teuchos::CommandLineProcessor clp(false);
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748);  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                    // manage parameters of xpetra

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build2DPoisson(nx, ny);

  Tpetra::Details::DeepCopyCounter::reset();
  Tpetra::Details::DeepCopyCounter::start();

  Level aLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(aLevel);
  aLevel.Request("A");
  aLevel.Set("A", A);

  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", dropFact);
  aLevel.Request(*aggFact);
  aggFact->Build(aLevel);

  Tpetra::Details::DeepCopyCounter::stop();

  if (Node::is_cpu) {
    TEST_EQUALITY(Tpetra::Details::DeepCopyCounter::get_count_different_space(), 0);
  }
#ifdef KOKKOS_HAS_SHARED_SPACE
  else {
    size_t targetNumDeepCopies = std::is_same_v<typename Node::memory_space, Kokkos::SharedSpace> ? 18 : 24;
    TEST_EQUALITY(Tpetra::Details::DeepCopyCounter::get_count_different_space(), targetNumDeepCopies);
  }
#else
  else {
    TEST_EQUALITY(Tpetra::Details::DeepCopyCounter::get_count_different_space(), 24);
  }
#endif

  stacked_timer->stopBaseTimer();
  Teuchos::StackedTimer::OutputOptions options;
  options.output_fraction = options.output_histogram = options.output_minmax = true;
  stacked_timer->report(out, comm, options);

}  // Aggregration

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Regression, H2D, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Regression, Aggregration, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
