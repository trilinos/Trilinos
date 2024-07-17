// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>

#include <MueLu_ConfigDefs.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_Utilities.hpp>

#include <MueLu_NoFactory.hpp>
#include <MueLu_Factory.hpp>

#include <MueLu_TestHelpers.hpp>

#include <MueLu_Level.hpp>
#include <MueLu_NullspaceFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>

#include <MueLu_SingleLevelFactoryBase.hpp>
#include <MueLu_Factory.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Level, SetCoreData, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  Level aLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(aLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(2);  // can be an empty operator

  aLevel.Set("Hitchhiker's Guide", 42);
  int fff = aLevel.Get<int>("Hitchhiker's Guide");
  TEST_EQUALITY(fff, 42);

  aLevel.Set("PI", 3.14159265);
  double ggg = aLevel.Get<double>("PI");
  TEST_EQUALITY(ggg, 3.14159265);
  TEST_EQUALITY(aLevel.IsAvailable("PI"), true);

  aLevel.Delete("PI", MueLu::NoFactory::get());
  TEST_EQUALITY(aLevel.IsAvailable("PI"), false);

  aLevel.Set("Hello MueLu", std::string("Greetings to MueMat"));
  std::string hhh = aLevel.Get<std::string>("Hello MueLu");
  TEST_EQUALITY(hhh, "Greetings to MueMat");

  aLevel.Set("A", A);
  RCP<Matrix> newA = aLevel.Get<RCP<Matrix> >("A");
  TEST_EQUALITY(newA, A);

  aLevel.Set("R", A);
  RCP<Matrix> newR = aLevel.Get<RCP<Matrix> >("R");
  TEST_EQUALITY(newR, A);  // TODO from JG: must be tested using another matrix !

  aLevel.Set("P", A);
  RCP<Matrix> newP = aLevel.Get<RCP<Matrix> >("P");
  TEST_EQUALITY(newP, A);

  aLevel.SetLevelID(42);
  TEST_EQUALITY(aLevel.GetLevelID(), 42);  // TODO: test default value of LevelID
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Level, NumRequests, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  Level aLevel;
  aLevel.SetLevelID(0);
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(2);
  aLevel.Set("A", A);

  RCP<FactoryManager> facManager = rcp(new FactoryManager());
  aLevel.SetFactoryManager(facManager);
  RCP<FactoryBase> factory = rcp(new CoalesceDropFactory());

  aLevel.Request("Graph", factory.get());
  aLevel.Request("Graph", factory.get());

  aLevel.Release("Graph", factory.get());
  TEST_EQUALITY(aLevel.IsRequested("Graph", factory.get()), true);

  aLevel.Release("Graph", factory.get());
  TEST_EQUALITY(aLevel.IsRequested("Graph", factory.get()), false);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Level, RequestRelease, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  Level l;
  l.SetLevelID(0);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(2);
  l.Set("A", A);

  RCP<FactoryManager> facManager = rcp(new FactoryManager());
  l.SetFactoryManager(facManager);

  RCP<FactoryBase> factory = rcp(new CoalesceDropFactory());

  l.Request("Graph", factory.get());
  TEST_EQUALITY(l.IsRequested("Graph", factory.get()), true);
  TEST_EQUALITY(l.IsAvailable("Graph", factory.get()), false);
  l.Release("Graph", factory.get());
  TEST_EQUALITY(l.IsRequested("Graph", factory.get()), false);
  TEST_EQUALITY(l.IsAvailable("Graph", factory.get()), false);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Level, RequestReleaseFactory, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  Level l;
  l.SetLevelID(0);

  RCP<FactoryManager> facManager = rcp(new FactoryManager());
  l.SetFactoryManager(facManager);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(2);
  l.Set("A", A);

  RCP<FactoryBase> graphFact = rcp(new CoalesceDropFactory());
  RCP<Factory> aggFact       = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", graphFact);

  l.Request("Aggregates", aggFact.get());
  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), true);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), false);

  TEST_EQUALITY(l.IsRequested("Graph", graphFact.get()), true);
  TEST_EQUALITY(l.IsAvailable("Graph", graphFact.get()), false);

  l.Release("Aggregates", aggFact.get());
  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), false);

  TEST_EQUALITY(l.IsRequested("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Graph", graphFact.get()), false);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Level, KeepFactory, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  Level l;
  l.SetLevelID(0);

  RCP<FactoryManager> facManager = rcp(new FactoryManager());
  l.SetFactoryManager(facManager);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(2);
  l.Set("A", A);

  RCP<Factory> graphFact = rcp(new CoalesceDropFactory());
  RCP<Factory> aggFact   = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", graphFact);

  l.Keep("Aggregates", aggFact.get());  // set keep flag
  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()), MueLu::Keep);
  l.Request("Aggregates", aggFact.get());
  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), true);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()), MueLu::Keep);

  TEST_EQUALITY(l.IsRequested("Graph", graphFact.get()), true);
  TEST_EQUALITY(l.IsAvailable("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.GetKeepFlag("Graph", graphFact.get()), 0);

  l.Release("Aggregates", aggFact.get());
  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()), MueLu::Keep);

  TEST_EQUALITY(l.IsRequested("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Graph", graphFact.get()), false);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Level, KeepAndBuildFactory, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  Level l;
  l.SetLevelID(0);  // level 0 necessary because of Nullspace factory

  RCP<FactoryManager> facManager = rcp(new FactoryManager());
  l.SetFactoryManager(facManager);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(144);
  l.Set("A", A);

  RCP<CoalesceDropFactory> graphFact       = rcp(new CoalesceDropFactory());
  RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", graphFact);

  l.Keep("Aggregates", aggFact.get());  // set keep flag
  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()), MueLu::Keep);
  l.Request("Aggregates", aggFact.get());
  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), true);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()), MueLu::Keep);

  aggFact->Build(l);

  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), true);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), true);
  TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()), MueLu::Keep);

  TEST_EQUALITY(l.IsRequested("Graph", graphFact.get()), true);
  TEST_EQUALITY(l.IsAvailable("Graph", graphFact.get()), true);
  TEST_EQUALITY(l.GetKeepFlag("Graph", graphFact.get()), 0);

  l.Release(*aggFact);  // release dependencies only

  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), true);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), true);
  TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()), MueLu::Keep);

  TEST_EQUALITY(l.IsRequested("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.GetKeepFlag("Graph", graphFact.get()), 0);

  l.Release("Aggregates", aggFact.get());

  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), true);
  TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()), MueLu::Keep);

  TEST_EQUALITY(l.IsRequested("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.GetKeepFlag("Graph", graphFact.get()), 0);

  l.RemoveKeepFlag("Aggregates", aggFact.get(), MueLu::Keep);

  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), false);

  TEST_EQUALITY(l.IsRequested("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.GetKeepFlag("Graph", graphFact.get()), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Level, KeepAndBuildFactory2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  Level l;
  l.SetLevelID(0);  // level 0 necessary because of Nullspace factory

  RCP<FactoryManager> facManager = rcp(new FactoryManager());
  l.SetFactoryManager(facManager);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(144);
  l.Set("A", A);

  RCP<CoalesceDropFactory> graphFact       = rcp(new CoalesceDropFactory());
  RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", graphFact);

  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), false);

  l.Request("Aggregates", aggFact.get());
  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), true);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), false);

  aggFact->Build(l);

  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), true);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), true);

  TEST_EQUALITY(l.IsRequested("Graph", graphFact.get()), true);
  TEST_EQUALITY(l.IsAvailable("Graph", graphFact.get()), true);
  TEST_EQUALITY(l.GetKeepFlag("Graph", graphFact.get()), 0);

  l.Release(*aggFact);

  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), true);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), true);

  TEST_EQUALITY(l.IsRequested("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.GetKeepFlag("Graph", graphFact.get()), 0);

  l.Release("Aggregates", aggFact.get());

  TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()), false);

  TEST_EQUALITY(l.IsRequested("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.IsAvailable("Graph", graphFact.get()), false);
  TEST_EQUALITY(l.GetKeepFlag("Graph", graphFact.get()), 0);

  /*l.RemoveKeepFlag("Aggregates", aggFact.get(), MueLu::Keep);

    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.GetKeepFlag("Graph",      graphFact.get()), 0);*/
}

// Helper class for unit test 'Level/CircularDependency'
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class CircularFactory : public MueLu::SingleLevelFactoryBase {
#include <MueLu_UseShortNames.hpp>

 public:
  CircularFactory(int value)
    : value_(value) {}

  virtual ~CircularFactory() {}

  void SetCircularFactory(RCP<FactoryBase> circular) { circular_ = circular; }

  void DeclareInput(Level& level) const {
    level.DeclareInput("data", circular_.get(), this);
  }

  void Build(Level& level) const {
    level.Set("data", value_, this);
    int value = level.Get<int>("data", circular_.get());
    level.Set("data", value + value_, this);
  }

 private:
  int value_;
  RCP<FactoryBase> circular_;

};  // class CircularFactory

//! Even though it's very special, a factory can generate data, that it requests itself.
//  Level must avoid self-recursive calls of Request
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Level, CircularDependencyWith1Factory, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef CircularFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> circular_factory_type;
  circular_factory_type A(2);

  A.SetCircularFactory(rcpFromRef(A));

  Level level;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);

  level.Request("data", &A);

  TEST_EQUALITY(level.Get<int>("data", &A), (2 + 2));

  level.Release("data", &A);
}

//! Test if circular dependencies between factories are allowed
//  This test corresponds to a use-case found developping repartitionning capability
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Level, CircularDependencyWithTwoFactories, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  typedef CircularFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> circular_factory_type;
  circular_factory_type A(2);
  circular_factory_type B(3);

  A.SetCircularFactory(rcpFromRef(B));
  B.SetCircularFactory(rcpFromRef(A));

  Level level;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);

  level.Request("data", &A);

  A.Build(level);

  TEST_EQUALITY(level.Get<int>("data", &A), (2 + 3) + 2);
  TEST_EQUALITY(level.Get<int>("data", &B), (2 + 3));

  level.Release(A);  // needed because A.Build(level) have been called manually
  level.Release("data", &A);
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Level, SetCoreData, SC, LO, GO, NO)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Level, NumRequests, SC, LO, GO, NO)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Level, RequestRelease, SC, LO, GO, NO)                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Level, RequestReleaseFactory, SC, LO, GO, NO)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Level, KeepFactory, SC, LO, GO, NO)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Level, KeepAndBuildFactory, SC, LO, GO, NO)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Level, KeepAndBuildFactory2, SC, LO, GO, NO)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Level, CircularDependencyWith1Factory, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Level, CircularDependencyWithTwoFactories, SC, LO, GO, NO)

}  // namespace MueLuTests
