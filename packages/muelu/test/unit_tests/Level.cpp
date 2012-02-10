#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_config.hpp"
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Utilities.hpp"

#include "MueLu_NoFactory.hpp"

#include "MueLu_TestHelpers.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(Level, SetCoreData)
  {
    out << "version: " << MueLu::Version() << std::endl;

    Level aLevel;
    TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSingleLevelHierarchy(aLevel);

    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(2); //can be an empty operator

    aLevel.Set("Hitchhiker's Guide", 42);
    int fff = aLevel.Get<int>("Hitchhiker's Guide");
    TEST_EQUALITY(fff, 42);

    aLevel.Set("PI",3.14159265);
    double ggg = aLevel.Get<double>("PI");
    TEST_EQUALITY(ggg, 3.14159265);
    TEST_EQUALITY(aLevel.IsAvailable("PI"), true);

    aLevel.Delete("PI");
    TEST_EQUALITY(aLevel.IsAvailable("PI"), false);

    aLevel.Set("Hello MueLu", std::string("Greetings to MueMat"));
    std::string hhh = aLevel.Get<std::string>("Hello MueLu");
    TEST_EQUALITY(hhh, "Greetings to MueMat");

    aLevel.Set("A",A);
    RCP<Operator> newA = aLevel.Get< RCP<Operator> >("A");
    TEST_EQUALITY(newA, A);

    aLevel.Set("R", A);
    RCP<Operator> newR = aLevel.Get< RCP<Operator> >("R");
    TEST_EQUALITY(newR, A); //TODO from JG: must be tested using another matrix !

    aLevel.Set("P", A);
    RCP<Operator> newP = aLevel.Get< RCP<Operator> >("P");
    TEST_EQUALITY(newP, A);

    aLevel.SetLevelID(42);
    TEST_EQUALITY(aLevel.GetLevelID(), 42); //TODO: test default value of LevelID
  }

  TEUCHOS_UNIT_TEST(Level, RequestRelease)
  {
    Level l;
    
    RCP<FactoryManager> facManager = rcp(new FactoryManager());
    l.SetFactoryManager(facManager);

    RCP<FactoryBase> factory = rcp(new CoalesceDropFactory(rcpFromRef(*MueLu::NoFactory::get())));
    
    l.Request("Graph", factory.get());
    TEST_EQUALITY(l.IsRequested("Graph", factory.get()), true);
    TEST_EQUALITY(l.IsAvailable("Graph", factory.get()), false);
    l.Release("Graph", factory.get());
    TEST_EQUALITY(l.IsRequested("Graph", factory.get()), false);
    TEST_EQUALITY(l.IsAvailable("Graph", factory.get()), false);
  }

  TEUCHOS_UNIT_TEST(Level, RequestReleaseFactory)
  {
    Level l;

    RCP<FactoryManager> facManager = rcp(new FactoryManager());
    l.SetFactoryManager(facManager);

    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(2);
    l.Set("A", A);

    RCP<FactoryBase> graphFact = rcp(new CoalesceDropFactory(rcpFromRef(*MueLu::NoFactory::get())));
    RCP<FactoryBase> aggFact   = rcp(new UCAggregationFactory(graphFact));
    
    l.Request("Aggregates", aggFact.get());
    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   true);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), true);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), false);

    l.Release("Aggregates", aggFact.get());
    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), false);
  }

  TEUCHOS_UNIT_TEST(Level, KeepFactory)
  {
    Level l;

    RCP<FactoryManager> facManager = rcp(new FactoryManager());
    l.SetFactoryManager(facManager);

    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(2);
    l.Set("A", A);

    RCP<FactoryBase> graphFact = rcp(new CoalesceDropFactory(rcpFromRef(*MueLu::NoFactory::get())));
    RCP<FactoryBase> aggFact   = rcp(new UCAggregationFactory(graphFact));

    l.Keep("Aggregates", aggFact.get());      // set keep flag
    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()),   MueLu::Keep);
    l.Request("Aggregates", aggFact.get());
    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   true);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()),   MueLu::Keep);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), true);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.GetKeepFlag("Graph",      graphFact.get()), 0);

    l.Release("Aggregates", aggFact.get());
    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()),   MueLu::Keep);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), false);
  }

  TEUCHOS_UNIT_TEST(Level, KeepAndBuildFactory)
  {
    Level l;
    l.SetLevelID(0); // level 0 necessary because of Nullspace factory

    RCP<FactoryManager> facManager = rcp(new FactoryManager());
    l.SetFactoryManager(facManager);

    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(144);
    l.Set("A", A);

    RCP<CoalesceDropFactory>  graphFact = rcp(new CoalesceDropFactory(rcpFromRef(*MueLu::NoFactory::get())));
    RCP<UCAggregationFactory> aggFact   = rcp(new UCAggregationFactory(graphFact));

    l.Keep("Aggregates", aggFact.get());      // set keep flag
    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()),   MueLu::Keep);
    l.Request("Aggregates", aggFact.get());
    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   true);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()),   MueLu::Keep);

    aggFact->Build(l);

    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   true);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   true);
    TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()),   MueLu::Keep);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), true);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), true);
    TEST_EQUALITY(l.GetKeepFlag("Graph",      graphFact.get()), 0);

    l.Release(*aggFact); // release dependencies only

    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   true);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   true);
    TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()),   MueLu::Keep);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.GetKeepFlag("Graph",      graphFact.get()), 0);

    l.Release("Aggregates", aggFact.get());

    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   true);
    TEST_EQUALITY(l.GetKeepFlag("Aggregates", aggFact.get()),   MueLu::Keep);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.GetKeepFlag("Graph",      graphFact.get()), 0);

    l.RemoveKeepFlag("Aggregates", aggFact.get(), MueLu::Keep);

    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.GetKeepFlag("Graph",      graphFact.get()), 0);

  }

  TEUCHOS_UNIT_TEST(Level, KeepAndBuildFactory2)
  {
    Level l;
    l.SetLevelID(0); // level 0 necessary because of Nullspace factory

    RCP<FactoryManager> facManager = rcp(new FactoryManager());
    l.SetFactoryManager(facManager);

    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(144);
    l.Set("A", A);

    RCP<CoalesceDropFactory>  graphFact = rcp(new CoalesceDropFactory(rcpFromRef(*MueLu::NoFactory::get())));
    RCP<UCAggregationFactory> aggFact   = rcp(new UCAggregationFactory(graphFact));

    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);

    l.Request("Aggregates", aggFact.get());
    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   true);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);

    aggFact->Build(l);

    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   true);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   true);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), true);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), true);
    TEST_EQUALITY(l.GetKeepFlag("Graph",      graphFact.get()), 0);

    l.Release(*aggFact);

    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   true);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   true);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.GetKeepFlag("Graph",      graphFact.get()), 0);

    l.Release("Aggregates", aggFact.get());

    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.GetKeepFlag("Graph",      graphFact.get()), 0);

    /*l.RemoveKeepFlag("Aggregates", aggFact.get(), MueLu::Keep);

    TEST_EQUALITY(l.IsRequested("Aggregates", aggFact.get()),   false);
    TEST_EQUALITY(l.IsAvailable("Aggregates", aggFact.get()),   false);

    TEST_EQUALITY(l.IsRequested("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.IsAvailable("Graph",      graphFact.get()), false);
    TEST_EQUALITY(l.GetKeepFlag("Graph",      graphFact.get()), 0);*/

  }


  // Helper class for unit test 'Level/CircularDependency'
  class CircularFactory : public MueLu::SingleLevelFactoryBase {
    
  public:
    
    CircularFactory(int value) : value_(value) { }
    
    virtual ~CircularFactory() { }
 
    void SetCircularFactory(RCP<FactoryBase> circular) { circular_ = circular; }

    void DeclareInput(Level &level) const {
      level.DeclareInput("data", circular_.get(), this);
    }

    void Build(Level& level) const {
      level.Set("data", value_, this);
      int value = level.Get<int>("data", circular_.get());
      level.Set("data", value+value_, this);
    }

  private:

    int value_;
    RCP<FactoryBase> circular_;

  }; //class CircularFactory

  //! Even though it's very special, a factory can generate data, that it requests itself.
  //  Level must avoid self-recursive calls of Request
  TEUCHOS_UNIT_TEST(Level, CircularDependencyWith1Factory) {
    CircularFactory A(2);
    
    A.SetCircularFactory(rcpFromRef(A));

    Level level; TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSingleLevelHierarchy(level);

    level.Request("data", &A);

    TEST_EQUALITY(level.Get<int>("data", &A), (2 + 2));

    level.Release("data", &A);
  }

  //! Test if circular dependencies between factories are allowed
  //  This test Corresponds to a use-case found developping repartitionning capability
  TEUCHOS_UNIT_TEST(Level, CircularDependencyWithTwoFactories) {
    CircularFactory A(2);
    CircularFactory B(3);
    
    A.SetCircularFactory(rcpFromRef(B));
    B.SetCircularFactory(rcpFromRef(A));

    Level level; TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSingleLevelHierarchy(level);

    level.Request("data", &A);

    A.Build(level);

    TEST_EQUALITY(level.Get<int>("data", &A), (2 + 3) + 2);
    TEST_EQUALITY(level.Get<int>("data", &B), (2 + 3));

    level.Release(A); // needed because A.Build(level) have been called manually
    level.Release("data", &A);
  }

} // namespace MueLuTests

