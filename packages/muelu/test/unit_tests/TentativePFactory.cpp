// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TentativePFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<TentativePFactory> tentPFact = rcp(new TentativePFactory);
  TEST_EQUALITY(tentPFact != Teuchos::null, true);

}  // Constructor

// TODO test BuildP

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TentativePFactory, MakeTentative_LapackQR, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real                  = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test QR with user-supplied nullspace" << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  GO nx         = 29;
  RCP<Matrix> A = test_factory::Build1DPoisson(nx);
  A->SetFixedBlockSize(1);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  LocalOrdinal NSdim         = 2;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->randomize();
  fineLevel.Set("Nullspace", nullSpace);
  fineLevel.Set("DofsPerNode", 1);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
  UncoupledAggFact->SetFactory("Graph", dropFact);

  UncoupledAggFact->SetMinNodesPerAggregate(3);
  UncoupledAggFact->SetMaxNeighAlreadySelected(0);
  UncoupledAggFact->SetOrdering("natural");

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", UncoupledAggFact);

  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
  TentativePFact->SetFactory("Aggregates", UncoupledAggFact);
  TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
  TentativePFact->SetFactory("CoarseMap", coarseMapFact);

  coarseLevel.Request("P", TentativePFact.get());          // request Ptent
  coarseLevel.Request("Nullspace", TentativePFact.get());  // request coarse nullspace
  coarseLevel.Request(*TentativePFact);
  TentativePFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> Ptent;
  coarseLevel.Get("P", Ptent, TentativePFact.get());

  RCP<MultiVector> coarseNullSpace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", TentativePFact.get());

  // check interpolation
  RCP<MultiVector> PtN = MultiVectorFactory::Build(Ptent->getRangeMap(), NSdim);
  Ptent->apply(*coarseNullSpace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0);

  RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  diff->putScalar(0.0);

  coarseLevel.Release("P", TentativePFact.get());          // release Ptent
  coarseLevel.Release("Nullspace", TentativePFact.get());  // release coarse nullspace

  // diff = fineNS + (-1.0)*(P*coarseNS) + 0*diff
  diff->update(1.0, *nullSpace, -1.0, *PtN, 0.0);

  Teuchos::Array<typename TST::magnitudeType> norms(NSdim);
  diff->norm2(norms);
  for (LocalOrdinal i = 0; i < NSdim; ++i) {
    out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
    TEST_EQUALITY(norms[i] < 100 * TMT::eps(), true);
  }

  Teuchos::ArrayRCP<const Scalar> col1 = coarseNullSpace->getData(0);
  Teuchos::ArrayRCP<const Scalar> col2 = coarseNullSpace->getData(1);
  TEST_EQUALITY(col1.size() == col2.size(), true);

}  // MakeTentative  Lapack QR

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TentativePFactory, MakeTentative, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real                  = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test QR with user-supplied nullspace" << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  GO nx         = 199;
  RCP<Matrix> A = test_factory::Build1DPoisson(nx);
  fineLevel.Request("A");
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->randomize();
  fineLevel.Set("Nullspace", nullSpace);
  fineLevel.Set("DofsPerNode", 1);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
  UncoupledAggFact->SetFactory("Graph", dropFact);
  UncoupledAggFact->SetMinNodesPerAggregate(3);
  UncoupledAggFact->SetMaxNeighAlreadySelected(0);
  UncoupledAggFact->SetOrdering("natural");

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", UncoupledAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
  TentativePFact->SetFactory("Aggregates", UncoupledAggFact);
  TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
  TentativePFact->SetFactory("CoarseMap", coarseMapFact);

  coarseLevel.Request("P", TentativePFact.get());  // request Ptent
  coarseLevel.Request("Nullspace", TentativePFact.get());
  coarseLevel.Request(*TentativePFact);
  TentativePFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> Ptent;
  coarseLevel.Get("P", Ptent, TentativePFact.get());

  RCP<MultiVector> coarseNullSpace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", TentativePFact.get());

  // check interpolation
  RCP<MultiVector> PtN = MultiVectorFactory::Build(Ptent->getRangeMap(), NSdim);
  Ptent->apply(*coarseNullSpace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0);

  RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  diff->putScalar(0.0);

  coarseLevel.Release("P", TentativePFact.get());  // release Ptent
  coarseLevel.Release("Nullspace", TentativePFact.get());

  // diff = fineNS + (-1.0)*(P*coarseNS) + 0*diff
  diff->update(1.0, *nullSpace, -1.0, *PtN, 0.0);

  Teuchos::Array<typename TST::magnitudeType> norms(NSdim);
  diff->norm2(norms);
  for (LocalOrdinal i = 0; i < NSdim; ++i) {
    out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
    TEST_EQUALITY(norms[i] < 100 * TMT::eps(), true);
  }

  // check normalization and orthogonality of prolongator columns
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Ptent, true, *Ptent, false, out);
  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > diagVec     = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (TST::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_EQUALITY(diagVec->norm1(), diagVec->getGlobalLength());
  TEST_EQUALITY(diagVec->normInf() - 1 < 100 * TMT::eps(), true);
  if (TST::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_EQUALITY(diagVec->meanValue(), TST::one());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

}  // MakeTentative

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TentativePFactory, MakeTentativeUsingDefaultNullSpace, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real                  = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test QR when nullspace isn't supplied by user" << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);

  GO nx         = 199;
  RCP<Matrix> A = test_factory::Build1DPoisson(nx);

  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

  fineLevel.Set("Coordinates", coordinates);

  RCP<TentativePFactory> tentativePFact = rcp(new TentativePFactory());

  coarseLevel.Request("P", tentativePFact.get());          // request Ptent
  coarseLevel.Request("Nullspace", tentativePFact.get());  // request coarse nullspace
  coarseLevel.Request(*tentativePFact);
  tentativePFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> Ptent;
  coarseLevel.Get("P", Ptent, tentativePFact.get());

  RCP<MultiVector> coarseNullSpace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", tentativePFact.get());

  coarseLevel.Release("P", tentativePFact.get());          // release Ptent
  coarseLevel.Release("Nullspace", tentativePFact.get());  // release coarse nullspace

  // grab default fine level nullspace (vector of all ones)
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), 1);
  nullSpace->putScalar(1.0);

  // check interpolation
  LocalOrdinal NSdim   = 1;
  RCP<MultiVector> PtN = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  Ptent->apply(*coarseNullSpace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0);

  RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  diff->putScalar(0.0);

  // diff = fineNS + (-1.0)*(P*coarseNS) + 0*diff
  diff->update(1.0, *nullSpace, -1.0, *PtN, 0.0);

  Teuchos::Array<typename TST::magnitudeType> norms(NSdim);
  diff->norm2(norms);
  for (LocalOrdinal i = 0; i < NSdim; ++i) {
    out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
    TEST_EQUALITY(norms[i] < 100 * TMT::eps(), true);
  }

}  // MakeTentativeUsingDefaultNullSpace

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TentativePFactory, NoQROption, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real                  = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test option that skips local QR factorizations" << std::endl;

  GlobalOrdinal nx = 20, ny = 20;

  // Describes the initial layout of matrix rows across processors.
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<const Map> map                  = Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node>(TestHelpers::Parameters::getLib(), "Cartesian2D", comm, galeriList);

  map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(map, 2);  // expand map for 2 DOFs per node

  galeriList.set("right boundary", "Neumann");
  galeriList.set("bottom boundary", "Neumann");
  galeriList.set("top boundary", "Neumann");
  galeriList.set("front boundary", "Neumann");
  galeriList.set("back boundary", "Neumann");
  galeriList.set("keepBCs", false);

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Elasticity2D", map, galeriList);
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(2);
  RCP<RealValuedMultiVector> coordinates = Pr->BuildCoords();

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", coordinates);

  LocalOrdinal NSdim         = 3;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->randomize();
  fineLevel.Set("Nullspace", nullSpace);
  fineLevel.Set("DofsPerNode", 2);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<UncoupledAggregationFactory> UnCoupledAggFact = rcp(new UncoupledAggregationFactory());
  UnCoupledAggFact->SetFactory("Graph", dropFact);

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", UnCoupledAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
  TentativePFact->SetFactory("Aggregates", UnCoupledAggFact);
  TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
  TentativePFact->SetFactory("CoarseMap", coarseMapFact);

  coarseLevel.Request("P", TentativePFact.get());  // request Ptent
  coarseLevel.Request("Nullspace", TentativePFact.get());
  coarseLevel.Request(*TentativePFact);
  Teuchos::ParameterList paramList;
  paramList.set("tentative: calculate qr", false);
  TentativePFact->SetParameterList(paramList);
  TentativePFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> Ptent;
  coarseLevel.Get("P", Ptent, TentativePFact.get());

  RCP<MultiVector> coarseNullSpace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", TentativePFact.get());

  // check interpolation
  RCP<MultiVector> PtN = MultiVectorFactory::Build(Ptent->getRangeMap(), NSdim);
  Ptent->apply(*coarseNullSpace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0);

  RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  diff->putScalar(0.0);

  coarseLevel.Release("P", TentativePFact.get());  // release Ptent
  coarseLevel.Release("Nullspace", TentativePFact.get());

  // diff = fineNS + (-1.0)*(P*coarseNS) + 0*diff
  diff->update(1.0, *nullSpace, -1.0, *PtN, 0.0);

  Teuchos::Array<typename TST::magnitudeType> norms(NSdim);
  diff->norm2(norms);
  for (LocalOrdinal i = 0; i < NSdim; ++i) {
    out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
    TEST_EQUALITY(norms[i] < 100 * TMT::eps(), true);
  }

}  // NoQR

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TentativePFactory, ConstantColumnSum, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real                  = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test option that skips local QR factorizations & uses constant column sum" << std::endl;

  GlobalOrdinal nx = 20, ny = 20;

  // Describes the initial layout of matrix rows across processors.
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<const Map> map                  = Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node>(TestHelpers::Parameters::getLib(), "Cartesian2D", comm, galeriList);

  map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(map, 1);

  galeriList.set("right boundary", "Neumann");
  galeriList.set("bottom boundary", "Neumann");
  galeriList.set("top boundary", "Neumann");
  galeriList.set("front boundary", "Neumann");
  galeriList.set("back boundary", "Neumann");
  galeriList.set("keepBCs", false);

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Laplace2D", map, galeriList);
  RCP<Matrix> A = Pr->BuildMatrix();

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  fineLevel.Request("A");
  fineLevel.Set("A", A);

  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->randomize();
  fineLevel.Set("Nullspace", nullSpace);
  fineLevel.Set("DofsPerNode", 1);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<UncoupledAggregationFactory> UnCoupledAggFact = rcp(new UncoupledAggregationFactory());
  UnCoupledAggFact->SetFactory("Graph", dropFact);

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", UnCoupledAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
  TentativePFact->SetFactory("Aggregates", UnCoupledAggFact);
  TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
  TentativePFact->SetFactory("CoarseMap", coarseMapFact);

  coarseLevel.Request("P", TentativePFact.get());  // request Ptent
  coarseLevel.Request("Nullspace", TentativePFact.get());
  coarseLevel.Request(*TentativePFact);
  Teuchos::ParameterList paramList;
  paramList.set("tentative: calculate qr", false);
  paramList.set("tentative: constant column sums", true);
  paramList.sublist("matrixmatrix: kernel params").set("compute global constants", true);
  TentativePFact->SetParameterList(paramList);
  TentativePFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> Ptent;
  coarseLevel.Get("P", Ptent, TentativePFact.get());

  RCP<MultiVector> coarseNullSpace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", TentativePFact.get());

  // This is "did it run test" not a "correctness check"

}  // ConstColSum

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TentativePFactory, NonStandardMaps, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  //#warning Unit test PgPFactory NonStandardMaps disabled
  //  return;
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
#if !defined(HAVE_MUELU_IFPACK)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Ifpack");
#endif
#if !defined(HAVE_MUELU_IFPACK2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Ifpack2");
#endif

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  // generate problem
  LocalOrdinal maxLevels = 3;
  // LocalOrdinal its=10;
  GlobalOrdinal nEle       = 63;
  GlobalOrdinal nIndexBase = 10;
  const RCP<const Map> map = MapFactory::Build(lib, nEle, nIndexBase, comm);
  RCP<Matrix> mtx          = Galeri::Xpetra::MatrixTraits<Map, CrsMatrixWrap>::Build(map, 3);

  LocalOrdinal NumMyElements                               = map->getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();
  GlobalOrdinal NumGlobalElements                          = map->getGlobalNumElements();
  assert(NumGlobalElements == nEle);

  GlobalOrdinal NumEntries;
  LocalOrdinal nnz = 2;
  std::vector<Scalar> Values(nnz);
  std::vector<GlobalOrdinal> Indices(nnz);

  Scalar a = 2.0;
  Scalar b = -1.0;
  Scalar c = -1.0;

  for (LocalOrdinal i = 0; i < NumMyElements; ++i) {
    if (MyGlobalElements[i] == nIndexBase) {
      // off-diagonal for first row
      Indices[0] = nIndexBase;
      NumEntries = 1;
      Values[0]  = c;
    } else if (MyGlobalElements[i] == nIndexBase + NumGlobalElements - 1) {
      // off-diagonal for last row
      Indices[0] = nIndexBase + NumGlobalElements - 2;
      NumEntries = 1;
      Values[0]  = b;
    } else {
      // off-diagonal for internal row
      Indices[0] = MyGlobalElements[i] - 1;
      Values[1]  = b;
      Indices[1] = MyGlobalElements[i] + 1;
      Values[0]  = c;
      NumEntries = 2;
    }

    // put the off-diagonal entries
    // Xpetra wants ArrayViews (sigh)
    Teuchos::ArrayView<Scalar> av(&Values[0], NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], NumEntries);
    mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

    // Put in the diagonal entry
    mtx->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                            Teuchos::tuple<Scalar>(a));

  }  // for (LocalOrdinal i = 0; i < NumMyElements; ++i)

  mtx->fillComplete(map, map);

  std::cout << map->getIndexBase() << std::endl;

  RCP<Matrix> Op = Teuchos::rcp_dynamic_cast<Matrix>(mtx);

  // build nullspace
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
  nullSpace->putScalar((Scalar)1.0);

  // fill hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy());
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  RCP<Level> Finest = H->GetLevel();  // first associate level with hierarchy (for defaultFactoryHandler!)

  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A", Op);                 // set fine level matrix
  Finest->Set("Nullspace", nullSpace);  // set null space information for finest level

  // define transfer operators
  RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
  UncoupledAggFact->SetMinNodesPerAggregate(3);
  UncoupledAggFact->SetMaxNeighAlreadySelected(0);
  UncoupledAggFact->SetOrdering("natural");

  RCP<TentativePFactory> Pfact = rcp(new TentativePFactory());
  RCP<Factory> Rfact           = rcp(new TransPFactory());
  RCP<RAPFactory> Acfact       = rcp(new RAPFactory());
  H->SetMaxCoarseSize(1);

  // setup smoothers
  Teuchos::ParameterList smootherParamList;
  smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
  smootherParamList.set("relaxation: sweeps", (LocalOrdinal)1);
  smootherParamList.set("relaxation: damping factor", (Scalar)1.0);
  RCP<SmootherPrototype> smooProto = rcp(new TrilinosSmoother("RELAXATION", smootherParamList));
  RCP<SmootherFactory> SmooFact    = rcp(new SmootherFactory(smooProto));
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  FactoryManager M;
  M.SetKokkosRefactor(false);
  M.SetFactory("P", Pfact);
  M.SetFactory("R", Rfact);
  M.SetFactory("A", Acfact);
  M.SetFactory("Ptent", Pfact);
  M.SetFactory("Aggregates", UncoupledAggFact);
  M.SetFactory("Smoother", SmooFact);
  M.SetFactory("CoarseSolver", SmooFact);

  H->Setup(M, 0, maxLevels);

  RCP<Level> coarseLevel = H->GetLevel(1);
  TEST_EQUALITY(coarseLevel->IsRequested("A", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("P", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("R", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel->IsAvailable("A", MueLu::NoFactory::get()), true);

  TEST_EQUALITY(coarseLevel->IsAvailable("P", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(coarseLevel->IsAvailable("R", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("A", MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("P", MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother", MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("R", MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(coarseLevel->IsRequested("P", Pfact.get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother", SmooFact.get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother", SmooFact.get()), false);
  TEST_EQUALITY(coarseLevel->IsRequested("R", Rfact.get()), false);
  TEST_EQUALITY(coarseLevel->IsAvailable("P", Pfact.get()), false);
  TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother", SmooFact.get()), false);
  TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother", SmooFact.get()), false);
  TEST_EQUALITY(coarseLevel->IsAvailable("R", Rfact.get()), false);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("P", Pfact.get()), 0);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother", SmooFact.get()), 0);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother", SmooFact.get()), 0);
  TEST_EQUALITY(coarseLevel->GetKeepFlag("R", Rfact.get()), 0);
  RCP<Level> coarseLevel2 = H->GetLevel(2);
  TEST_EQUALITY(coarseLevel2->IsRequested("A", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel2->IsRequested("P", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel2->IsRequested("R", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel2->IsRequested("PreSmoother", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel2->IsRequested("PostSmoother", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(coarseLevel2->IsAvailable("A", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(coarseLevel2->IsAvailable("P", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false);  // coarse level only has presmoother = coarse solver
  TEST_EQUALITY(coarseLevel2->IsAvailable("R", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(coarseLevel2->GetKeepFlag("A", MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(coarseLevel2->GetKeepFlag("P", MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(coarseLevel2->GetKeepFlag("PreSmoother", MueLu::NoFactory::get()), MueLu::Final);
  // TEST_EQUALITY(coarseLevel2->GetKeepFlag("PostSmoother",MueLu::NoFactory::get()), MueLu::Final); // coarse level only has presmoother = coarse solver
  TEST_EQUALITY(coarseLevel2->GetKeepFlag("R", MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(coarseLevel2->IsRequested("P", Pfact.get()), false);
  TEST_EQUALITY(coarseLevel2->IsRequested("R", Rfact.get()), false);
  TEST_EQUALITY(coarseLevel2->IsAvailable("P", Pfact.get()), false);
  TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother", SmooFact.get()), false);
  TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother", SmooFact.get()), false);
  TEST_EQUALITY(coarseLevel2->IsAvailable("R", Rfact.get()), false);
  TEST_EQUALITY(coarseLevel2->GetKeepFlag("P", Pfact.get()), 0);
  TEST_EQUALITY(coarseLevel2->GetKeepFlag("PreSmoother", SmooFact.get()), 0);
  TEST_EQUALITY(coarseLevel2->GetKeepFlag("PostSmoother", SmooFact.get()), 0);
  TEST_EQUALITY(coarseLevel2->GetKeepFlag("R", Rfact.get()), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TentativePFactory, PtentEpetraVsTpetra, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_EPETRA_SCOPE_TPETRA_IS_DEFAULT(Scalar, GlobalOrdinal, Node);
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT) && defined(HAVE_MUELU_IFPACK) && defined(HAVE_MUELU_IFPACK2)

  using TST            = Teuchos::ScalarTraits<Scalar>;
  using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;
  using real           = typename TST::coordinateType;
  typedef Xpetra::MultiVector<real, LO, GO, NO> RealValuedMultiVector;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test QR when nullspace isn't supplied by user" << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::Array<magnitude_type> results(2);

  // run test only on 1 proc
  if (comm->getSize() == 1) {
    Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

    // run Epetra and Tpetra test
    for (int run = 0; run < 2; run++) {
      if (run == 0)
        lib = Xpetra::UseEpetra;
      else
        lib = Xpetra::UseTpetra;

      // generate problem
      LocalOrdinal maxLevels   = 3;
      LocalOrdinal its         = 10;
      GlobalOrdinal nEle       = 63;
      const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
      Teuchos::ParameterList matrixParameters;
      matrixParameters.set("nx", nEle);
      RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
          Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Laplace1D", map, matrixParameters);
      RCP<Matrix> Op                         = Pr->BuildMatrix();
      RCP<RealValuedMultiVector> coordinates = Pr->BuildCoords();

      // build nullspace
      RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
      nullSpace->putScalar(TST::one());
      Teuchos::Array<magnitude_type> norms(1);
      nullSpace->norm1(norms);
      if (comm->getRank() == 0)
        out << "||NS|| = " << norms[0] << std::endl;

      // fill hierarchy
      RCP<Hierarchy> H = rcp(new Hierarchy());
      H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
      RCP<Level> Finest = H->GetLevel();  // first associate level with hierarchy (for defaultFactoryHandler!)

      Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
      Finest->Set("A", Op);                 // set fine level matrix
      Finest->Set("Nullspace", nullSpace);  // set null space information for finest level
      // Finest->Set("Coordinates", coordinates);  // set coordinates for finest level

      // define transfer operators
      RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
      UncoupledAggFact->SetMinNodesPerAggregate(3);
      UncoupledAggFact->SetMaxNeighAlreadySelected(0);
      UncoupledAggFact->SetOrdering("natural");

      RCP<TentativePFactory> Pfact = rcp(new TentativePFactory());
      RCP<Factory> Rfact           = rcp(new TransPFactory());
      RCP<RAPFactory> Acfact       = rcp(new RAPFactory());
      H->SetMaxCoarseSize(1);

      // setup smoothers
      Teuchos::ParameterList smootherParamList;
      smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
      smootherParamList.set("relaxation: sweeps", (LocalOrdinal)1);
      smootherParamList.set("relaxation: damping factor", TST::one());
      RCP<SmootherPrototype> smooProto = rcp(new TrilinosSmoother("RELAXATION", smootherParamList));
      RCP<SmootherFactory> SmooFact    = rcp(new SmootherFactory(smooProto));
      Acfact->setVerbLevel(Teuchos::VERB_HIGH);

      RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

      FactoryManager M;
      M.SetKokkosRefactor(false);
      M.SetFactory("P", Pfact);
      M.SetFactory("R", Rfact);
      M.SetFactory("A", Acfact);
      M.SetFactory("Ptent", Pfact);
      M.SetFactory("Aggregates", UncoupledAggFact);
      M.SetFactory("Smoother", SmooFact);
      M.SetFactory("CoarseSolver", coarseSolveFact);
      // M.SetFactory("Coordinates", Pfact);

      H->Setup(M, 0, maxLevels);

      // test some basic multgrid data
      RCP<Level> coarseLevel = H->GetLevel(1);
      RCP<Matrix> P1         = coarseLevel->Get<RCP<Matrix> >("P");
      RCP<Matrix> R1         = coarseLevel->Get<RCP<Matrix> >("R");
      TEST_EQUALITY(P1->getGlobalNumRows(), 63);
      TEST_EQUALITY(P1->getGlobalNumCols(), 21);
      TEST_EQUALITY(R1->getGlobalNumRows(), 21);
      TEST_EQUALITY(R1->getGlobalNumCols(), 63);
      RCP<Level> coarseLevel2 = H->GetLevel(2);
      RCP<Matrix> P2          = coarseLevel2->Get<RCP<Matrix> >("P");
      RCP<Matrix> R2          = coarseLevel2->Get<RCP<Matrix> >("R");
      TEST_EQUALITY(P2->getGlobalNumRows(), 21);
      TEST_EQUALITY(P2->getGlobalNumCols(), 7);
      TEST_EQUALITY(R2->getGlobalNumRows(), 7);
      TEST_EQUALITY(R2->getGlobalNumCols(), 21);

      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*P1, true, *P1, false, out);
      Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > diagVec     = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(PtentTPtent->getRowMap());
      PtentTPtent->getLocalDiagCopy(*diagVec);
      TEST_EQUALITY(diagVec->norm1() - diagVec->getGlobalLength() < 100 * TMT::eps(), true);
      TEST_EQUALITY(diagVec->normInf() - TMT::one() < 100 * TMT::eps(), true);
      TEST_EQUALITY(TST::magnitude(diagVec->meanValue()) - TMT::one() < 100 * TMT::eps(), true);
      TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

      // Define RHS
      RCP<MultiVector> X   = MultiVectorFactory::Build(map, 1);
      RCP<MultiVector> RHS = MultiVectorFactory::Build(map, 1);

      X->putScalar(1.0);
      X->norm2(norms);
      if (comm->getRank() == 0)
        out << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

      Op->apply(*X, *RHS, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);

      // Use AMG directly as an iterative method
      {
        X->putScalar((Scalar)0.0);

        H->Iterate(*RHS, *X, its);

        X->norm2(norms);
        if (comm->getRank() == 0)
          out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
        results[run] = norms[0];
      }
    }

    TEST_FLOATING_EQUALITY(results[0], results[1], 100 * TMT::eps());  // check results of EPETRA vs TPETRA
  }                                                                    // comm->getSize == 1

#else
  out << "Skipping test because some required packages are not enabled (Tpetra, Epetra, EpetraExt, Ifpack, Ifpack2)." << std::endl;
#endif

}  // TentativePFactory_EpetraVsTpetra

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TentativePFactory, Constructor, Scalar, LO, GO, Node)                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TentativePFactory, MakeTentative_LapackQR, Scalar, LO, GO, Node)             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TentativePFactory, MakeTentative, Scalar, LO, GO, Node)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TentativePFactory, MakeTentativeUsingDefaultNullSpace, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TentativePFactory, NoQROption, Scalar, LO, GO, Node)                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TentativePFactory, ConstantColumnSum, Scalar, LO, GO, Node)                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TentativePFactory, NonStandardMaps, Scalar, LO, GO, Node)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TentativePFactory, PtentEpetraVsTpetra, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
