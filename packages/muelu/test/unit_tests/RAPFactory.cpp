// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_FactoryManager.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RAPFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<RAPFactory> rapFactory = rcp(new RAPFactory);
  TEST_EQUALITY(rapFactory != Teuchos::null, true);

  out << *rapFactory << std::endl;
}  // Constructor test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RAPFactory, Correctness, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<Scalar, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

  GO nx          = 27 * comm->getSize();
  RCP<Matrix> Op = TestHelpers::TestFactory<Scalar, LO, GO, NO>::Build1DPoisson(nx);
  fineLevel.Set("A", Op);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", Op->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  TentativePFactory tentpFactory;
  SaPFactory sapFactory;
  sapFactory.SetFactory("P", rcpFromRef(tentpFactory));
  TransPFactory transPFactory;
  transPFactory.SetFactory("P", rcpFromRef(sapFactory));  // todo:rcpFromRef

  coarseLevel.Request("P", &sapFactory);
  coarseLevel.Request("R", &transPFactory);

  coarseLevel.Request(sapFactory);
  coarseLevel.Request(transPFactory);
  sapFactory.Build(fineLevel, coarseLevel);
  transPFactory.Build(fineLevel, coarseLevel);

  RAPFactory rap;
  rap.SetFactory("P", rcpFromRef(sapFactory));
  rap.SetFactory("R", rcpFromRef(transPFactory));

  coarseLevel.Request(rap);

  coarseLevel.Request("A", &rap);
  rap.Build(fineLevel, coarseLevel);

  RCP<Matrix> A = fineLevel.Get<RCP<Matrix> >("A");
  RCP<Matrix> P = coarseLevel.Get<RCP<Matrix> >("P", &sapFactory);
  RCP<Matrix> R = coarseLevel.Get<RCP<Matrix> >("R", &transPFactory);

  RCP<MultiVector> workVec1 = MultiVectorFactory::Build(P->getRangeMap(), 1);
  RCP<MultiVector> workVec2 = MultiVectorFactory::Build(Op->getRangeMap(), 1);
  RCP<MultiVector> result1  = MultiVectorFactory::Build(R->getRangeMap(), 1);
  RCP<MultiVector> X        = MultiVectorFactory::Build(P->getDomainMap(), 1);
  X->randomize();

  // Calculate result1 = R*(A*(P*X))
  P->apply(*X, *workVec1, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  Op->apply(*workVec1, *workVec2, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  R->apply(*workVec2, *result1, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);

  RCP<Matrix> coarseOp = coarseLevel.Get<RCP<Matrix> >("A", &rap);

  // Calculate result2 = (R*A*P)*X
  RCP<MultiVector> result2 = MultiVectorFactory::Build(R->getRangeMap(), 1);
  coarseOp->apply(*X, *result2, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);

  Teuchos::Array<typename TST::magnitudeType> normX(1), normResult1(1), normResult2(1);
  X->norm2(normX);
  out << "This test checks the correctness of the Galerkin triple "
      << "matrix product by comparing (RAP)*X to R(A(P*X))." << std::endl;
  out << "||X||_2 = " << normX << std::endl;
  result1->norm2(normResult1);
  result2->norm2(normResult2);
  TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 100 * TMT::eps());

}  // Correctness test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RAPFactory, ImplicitTranspose, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  if (comm->getSize() > 1 && TestHelpers::Parameters::getLib() == Xpetra::UseEpetra) {
    out << "Skipping ImplicitTranspose test for Epetra and #proc>1" << std::endl;
    return;
  }

  // build test-specific default factory manager
  RCP<FactoryManager> defManager = rcp(new FactoryManager());
  defManager->SetKokkosRefactor(false);
  defManager->SetFactory("A", rcp(MueLu::NoFactory::get(), false));              // dummy factory for A
  defManager->SetFactory("Nullspace", rcp(new NullspaceFactory()));              // real null space factory for Ptent
  defManager->SetFactory("Graph", rcp(new CoalesceDropFactory()));               // real graph factory for Ptent
  defManager->SetFactory("Aggregates", rcp(new UncoupledAggregationFactory()));  // real aggregation factory for Ptent

  Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<Scalar, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

  // overwrite default factory manager
  fineLevel.SetFactoryManager(defManager);
  coarseLevel.SetFactoryManager(defManager);

  GO nx          = 19 * comm->getSize();
  RCP<Matrix> Op = TestHelpers::TestFactory<Scalar, LO, GO, NO>::Build1DPoisson(nx);
  fineLevel.Set("A", Op);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", Op->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  TentativePFactory tentpFactory;
  SaPFactory sapFactory;
  sapFactory.SetFactory("P", rcpFromRef(tentpFactory));
  TransPFactory transPFactory;
  transPFactory.SetFactory("P", rcpFromRef(sapFactory));
  coarseLevel.Request("P", &sapFactory);
  coarseLevel.Request("R", &transPFactory);

  coarseLevel.Request(sapFactory);
  coarseLevel.Request(transPFactory);
  sapFactory.Build(fineLevel, coarseLevel);
  transPFactory.Build(fineLevel, coarseLevel);
  RAPFactory rap;
  Teuchos::ParameterList rapList = *(rap.GetValidParameterList());
  rapList.set("transpose: use implicit", true);
  rap.SetParameterList(rapList);
  rap.SetFactory("P", rcpFromRef(sapFactory));
  rap.SetFactory("R", rcpFromRef(transPFactory));
  coarseLevel.Request("A", &rap);

  coarseLevel.Request(rap);
  rap.Build(fineLevel, coarseLevel);

  RCP<Matrix> A = fineLevel.Get<RCP<Matrix> >("A");
  RCP<Matrix> P = coarseLevel.Get<RCP<Matrix> >("P", &sapFactory);
  RCP<Matrix> R = coarseLevel.Get<RCP<Matrix> >("R", &transPFactory);

  // std::string filename = "A.dat";
  // Utils::Write(filename,Op);
  // filename = "P.dat";
  // Utils::Write(filename,P);

  RCP<MultiVector> workVec1 = MultiVectorFactory::Build(P->getRangeMap(), 1);
  RCP<MultiVector> workVec2 = MultiVectorFactory::Build(Op->getRangeMap(), 1);
  RCP<MultiVector> result1  = MultiVectorFactory::Build(P->getDomainMap(), 1);
  RCP<MultiVector> X        = MultiVectorFactory::Build(P->getDomainMap(), 1);
  X->randomize();
  // out.precision(12);
  // out.setOutputToRootOnly(-1);
  // X->describe(out,Teuchos::VERB_EXTREME);

  // Calculate result1 = P^T*(A*(P*X))
  P->apply(*X, *workVec1, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  Op->apply(*workVec1, *workVec2, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  P->apply(*workVec2, *result1, Teuchos::TRANS, (Scalar)1.0, (Scalar)0.0);

  RCP<Matrix> coarseOp = coarseLevel.Get<RCP<Matrix> >("A", &rap);

  // Calculate result2 = (R*A*P)*X
  RCP<MultiVector> result2 = MultiVectorFactory::Build(P->getDomainMap(), 1);
  coarseOp->apply(*X, *result2, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);

  Teuchos::Array<typename TST::magnitudeType> normX(1), normResult1(1), normResult2(1);
  X->norm2(normX);
  out << "This test checks the correctness of the Galerkin triple "
      << "matrix product by comparing (RAP)*X to R(A(P*X)), where R is the implicit transpose of P." << std::endl;
  out << "||X||_2 = " << normX << std::endl;
  result1->norm2(normResult1);
  result2->norm2(normResult2);
  TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 100 * TMT::eps());

}  // ImplicitTranspose test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RAPFactory, FixZeroDiagonals, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  out << "This test checks that the option \"rap: fix zero diagonals\" works properly." << std::endl;

  typedef typename Teuchos::ScalarTraits<Scalar> TST;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  // build test-specific default factory manager
  RCP<FactoryManager> defManager = rcp(new FactoryManager());
  defManager->SetKokkosRefactor(false);
  defManager->SetFactory("A", rcp(MueLu::NoFactory::get(), false));              // dummy factory for A
  defManager->SetFactory("Nullspace", rcp(new NullspaceFactory()));              // real null space factory for Ptent
  defManager->SetFactory("Graph", rcp(new CoalesceDropFactory()));               // real graph factory for Ptent
  defManager->SetFactory("Aggregates", rcp(new UncoupledAggregationFactory()));  // real aggregation factory for Ptent
  defManager->SetFactory("P", rcp(MueLu::NoFactory::get(), false));              // dummy factory for P
  defManager->SetFactory("R", rcp(MueLu::NoFactory::get(), false));              // dummy factory for R

  Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<Scalar, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

  // overwrite default factory manager
  fineLevel.SetFactoryManager(defManager);
  coarseLevel.SetFactoryManager(defManager);

  GO nx          = 19 * comm->getSize();
  RCP<Matrix> Op = TestHelpers::TestFactory<Scalar, LO, GO, NO>::Build1DPoisson(nx);
  fineLevel.Set("A", Op);

  // Create diagonal matrix, same dimensions as A, but with some zeros on the diagonal
  RCP<const Map> rowMap = Op->getRowMap();
  Teuchos::Array<GO> indout(1);
  Teuchos::Array<SC> valout(1);
  RCP<Matrix> diagMat  = MatrixFactory::Build(rowMap, 1);
  RCP<Matrix> zdiagMat = MatrixFactory::Build(rowMap, 1);
  for (size_t r = 0; r < rowMap->getLocalNumElements(); r++) {
    GO grid   = rowMap->getGlobalElement(r);
    indout[0] = grid;
    (r == 0) ? valout[0] = TST::zero() : valout[0] = TST::one();
    zdiagMat->insertGlobalValues(grid, indout(), valout());
    valout[0] = TST::one();
    diagMat->insertGlobalValues(grid, indout(), valout());
  }
  zdiagMat->fillComplete();
  diagMat->fillComplete();

  // this effectively zero's just the row of R*A*P
  coarseLevel.Set("P", diagMat, MueLu::NoFactory::get());
  coarseLevel.Set("R", zdiagMat, MueLu::NoFactory::get());

  RAPFactory rap;
  Teuchos::ParameterList rapList = *(rap.GetValidParameterList());
  rapList.set("rap: fix zero diagonals", true);
  rap.SetParameterList(rapList);
  coarseLevel.Request("A", &rap);

  coarseLevel.Request(rap);
  rap.Build(fineLevel, coarseLevel);

  RCP<Matrix> A = fineLevel.Get<RCP<Matrix> >("A");
  RCP<Matrix> P = coarseLevel.Get<RCP<Matrix> >("P");
  RCP<Matrix> R = coarseLevel.Get<RCP<Matrix> >("R");

  RCP<Matrix> coarseOp = coarseLevel.Get<RCP<Matrix> >("A", &rap);
  // filename = "Ac.dat";
  // Xpetra::IO<SC,LO,GO,NO>::Write(filename,*coarseOp);

  // The diagonal repair should result in a single nonzero entry with scalar value one in local row 0,
  // There should be no repeated column indices.
  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const Scalar> vals;
  coarseOp->getLocalRowView(0, indices, vals);
  for (int j = 0; j < indices.size(); ++j) {
    for (int i = j + 1; i < indices.size(); ++i) {
      out << "checking indices[" << j << "] and indices[" << i << "]" << std::endl;
      TEST_INEQUALITY(indices[i], indices[j]);
    }
  }
  Scalar sum = Teuchos::ScalarTraits<SC>::zero();
  for (int j = 0; j < vals.size(); ++j)
    sum += vals[j];
  TEST_EQUALITY(sum, Teuchos::ScalarTraits<SC>::one());
  const RCP<const Map> colmap = coarseOp->getColMap();
  const RCP<const Map> rowmap = coarseOp->getRowMap();
  for (int j = 0; j < vals.size(); ++j) {
    if (vals[j] == Teuchos::ScalarTraits<SC>::one()) {
      out << "checking that nonzero entry is on diagonal" << std::endl;
      TEST_EQUALITY(colmap->getGlobalElement(indices[j]), rowmap->getGlobalElement(0));
    }
  }

}  // FixZeroDiagonals test

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RAPFactory, RelativeBoost, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  out << "This test checks that the option \"rap: relative diagonal floor\" works properly." << std::endl;

  typedef typename Teuchos::ScalarTraits<Scalar> TST;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  // build test-specific default factory manager
  RCP<FactoryManager> defManager = rcp(new FactoryManager());
  defManager->SetKokkosRefactor(false);
  defManager->SetFactory("A", rcp(MueLu::NoFactory::get(), false));              // dummy factory for A
  defManager->SetFactory("Nullspace", rcp(new NullspaceFactory()));              // real null space factory for Ptent
  defManager->SetFactory("Graph", rcp(new CoalesceDropFactory()));               // real graph factory for Ptent
  defManager->SetFactory("Aggregates", rcp(new UncoupledAggregationFactory()));  // real aggregation factory for Ptent
  defManager->SetFactory("P", rcp(MueLu::NoFactory::get(), false));              // dummy factory for P
  defManager->SetFactory("R", rcp(MueLu::NoFactory::get(), false));              // dummy factory for R

  Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<Scalar, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

  // overwrite default factory manager
  fineLevel.SetFactoryManager(defManager);
  coarseLevel.SetFactoryManager(defManager);

  GO nx          = 19 * comm->getSize();
  RCP<Matrix> Op = TestHelpers::TestFactory<Scalar, LO, GO, NO>::Build1DPoisson(nx);
  fineLevel.Set("A", Op);

  // Create diagonal matrix, same dimensions as A, but with some zeros on the diagonal
  RCP<const Map> rowMap = Op->getRowMap();
  Teuchos::Array<GO> indout(1);
  Teuchos::Array<SC> valout(1);
  RCP<Matrix> diagMat  = MatrixFactory::Build(rowMap, 1);
  RCP<Matrix> zdiagMat = MatrixFactory::Build(rowMap, 1);
  for (size_t r = 0; r < rowMap->getLocalNumElements(); r++) {
    GO grid   = rowMap->getGlobalElement(r);
    indout[0] = grid;
    (r == 0) ? valout[0] = TST::zero() : valout[0] = TST::one();
    zdiagMat->insertGlobalValues(grid, indout(), valout());
    valout[0] = TST::one();
    diagMat->insertGlobalValues(grid, indout(), valout());
  }
  zdiagMat->fillComplete();
  diagMat->fillComplete();

  // this effectively zero's just the row of R*A*P
  coarseLevel.Set("P", diagMat, MueLu::NoFactory::get());
  coarseLevel.Set("R", zdiagMat, MueLu::NoFactory::get());

  RAPFactory rap;
  Teuchos::ParameterList rapList = *(rap.GetValidParameterList());
  Teuchos::Array<double> relativeFloor(1);
  relativeFloor[0] = 0.5;
  rapList.set("rap: relative diagonal floor", relativeFloor);
  rap.SetParameterList(rapList);
  coarseLevel.Request("A", &rap);

  coarseLevel.Request(rap);
  rap.Build(fineLevel, coarseLevel);

  RCP<Matrix> A = fineLevel.Get<RCP<Matrix> >("A");
  RCP<Matrix> P = coarseLevel.Get<RCP<Matrix> >("P");
  RCP<Matrix> R = coarseLevel.Get<RCP<Matrix> >("R");

  RCP<Matrix> coarseOp = coarseLevel.Get<RCP<Matrix> >("A", &rap);
  std::string filename = "Ac.dat";
  Xpetra::IO<SC, LO, GO, NO>::Write(filename, *coarseOp);

  // The diagonal repair should result in a single nonzero entry with scalar value one in local row 0,
  // There should be no repeated column indices.
  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const Scalar> vals;
  coarseOp->getLocalRowView(0, indices, vals);
  for (int j = 0; j < indices.size(); ++j) {
    for (int i = j + 1; i < indices.size(); ++i) {
      out << "checking indices[" << j << "] and indices[" << i << "]" << std::endl;
      TEST_INEQUALITY(indices[i], indices[j]);
    }
  }
  Scalar sum = Teuchos::ScalarTraits<SC>::zero();
  for (int j = 0; j < vals.size(); ++j)
    sum += vals[j];
  TEST_EQUALITY(sum, Teuchos::ScalarTraits<SC>::one());
  const RCP<const Map> colmap = coarseOp->getColMap();
  const RCP<const Map> rowmap = coarseOp->getRowMap();
  for (int j = 0; j < vals.size(); ++j) {
    if (vals[j] == Teuchos::ScalarTraits<SC>::one()) {
      out << "checking that nonzero entry is on diagonal" << std::endl;
      TEST_EQUALITY(colmap->getGlobalElement(indices[j]), rowmap->getGlobalElement(0));
    }
  }

}  // Relative Boost test

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RAPFactory, Constructor, Scalar, LO, GO, Node)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RAPFactory, Correctness, Scalar, LO, GO, Node)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RAPFactory, ImplicitTranspose, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RAPFactory, FixZeroDiagonals, Scalar, LO, GO, Node)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RAPFactory, RelativeBoost, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
