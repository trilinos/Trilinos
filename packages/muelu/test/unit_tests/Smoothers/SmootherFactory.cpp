#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_TestHelpers.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_FakeSmootherPrototype.hpp"

#include "MueLu_UseDefaultTypes.hpp"

namespace MueLuTests {

  namespace SmootherFactory {
#include "MueLu_UseShortNames.hpp"

    using Teuchos::RCP;
    using Teuchos::rcp;
    using namespace MueLu::Exceptions;
    
    ////////////////////////////////////////////////////////////////////////////////
    // Helper functions: class invariant
    ////////////////////////////////////////////////////////////////////////////////

    // Check the consistency of the internal state of a SmootherFactory object
    void check(RCP<const SmootherFactory> smootherFactory, RCP<const SmootherPrototype> preSmooProto, RCP<const SmootherPrototype> postSmooProto, Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherPrototype> getPreSmooProto, getPostSmooProto;
      smootherFactory->GetSmootherPrototypes(getPreSmooProto, getPostSmooProto);

      TEST_EQUALITY(preSmooProto,  getPreSmooProto);
      TEST_EQUALITY(postSmooProto, getPostSmooProto);

      // Check if internal prototypes are still prototypes (IsSetup() == false)
      if (preSmooProto  != Teuchos::null) TEST_EQUALITY(preSmooProto->IsSetup(),  false);
      if (postSmooProto != Teuchos::null) TEST_EQUALITY(postSmooProto->IsSetup(), false);
    }

    // Similar to check(), but with INEQUALITY.
    // Do not check if internal prototypes verifies IsSetup() == false
    void ineqcheck(RCP<const SmootherFactory> smootherFactory, RCP<const SmootherBase> preSmooProto, RCP<const SmootherBase> postSmooProto, Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherPrototype> getPreSmooProto, getPostSmooProto;
      smootherFactory->GetSmootherPrototypes(getPreSmooProto, getPostSmooProto);
      TEST_INEQUALITY(preSmooProto,  getPreSmooProto);
      TEST_INEQUALITY(postSmooProto, getPostSmooProto);
    }
    
    ////////////////////////////////////////////////////////////////////////////////
    // Helper functions: apply a test to a collection of test cases
    ////////////////////////////////////////////////////////////////////////////////

    // Apply a test to a list of cases (one argument)
    void testCollection(void (*func)(RCP<SmootherPrototype>, Teuchos::FancyOStream &, bool &), Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherPrototype> smooProto = rcp( new FakeSmootherPrototype() );
      
      // tests
      func(smooProto, out, success);
      func(Teuchos::null, out, success);
    }

    // Apply a test to a list of invalid cases (one argument)
    void testInvalidCollection(void (*func)(RCP<SmootherPrototype>, Teuchos::FancyOStream &, bool &), Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherPrototype> invalidSmooProto = rcp( new FakeSmootherPrototype() ); Level level; invalidSmooProto->Setup(level);

      // tests
      func(invalidSmooProto, out, success);
    }

    // Apply a test to a list of cases (two argument)
    void testCollection(void (*func)(RCP<SmootherPrototype>, RCP<SmootherPrototype>, Teuchos::FancyOStream &, bool &), Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherPrototype> smooProtoA = rcp( new FakeSmootherPrototype() );
      RCP<SmootherPrototype> smooProtoB = rcp( new FakeSmootherPrototype() );

      // tests
      func(Teuchos::null, Teuchos::null, out, success);
      func(smooProtoA, smooProtoB, out, success);
      func(smooProtoA, Teuchos::null, out, success);
      func(Teuchos::null, smooProtoB, out, success);
      func(smooProtoA, smooProtoA, out, success);
    }

    // Apply a test to a list of invalid cases (two argument)
    void testInvalidCollection(void (*func)(RCP<SmootherPrototype>, RCP<SmootherPrototype>, Teuchos::FancyOStream &, bool &), Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherPrototype> smooProtoA = rcp( new FakeSmootherPrototype() );
      RCP<SmootherPrototype> smooProtoB = rcp( new FakeSmootherPrototype() );
      RCP<SmootherPrototype> invalidSmooProto = rcp( new FakeSmootherPrototype() ); Level level; invalidSmooProto->Setup(level);

      // tests
      func(invalidSmooProto, smooProtoB, out, success);
      func(smooProtoA, invalidSmooProto, out, success);
      func(invalidSmooProto, Teuchos::null, out, success);
      func(Teuchos::null, invalidSmooProto, out, success);
      func(invalidSmooProto, invalidSmooProto, out, success);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Test: Constructor_OneArg
    ////////////////////////////////////////////////////////////////////////////////

    // Build an object and test internal state (one argument)
    void testConstructor1(RCP<SmootherPrototype> smooProto, Teuchos::FancyOStream & out, bool & success) {
      check(rcp( new SmootherFactory(smooProto) ), smooProto, smooProto, out, success); // One argument == same pre and post smoother prototype
    }

    // Test with invalid input argument
    void testInvalidConstructor1(RCP<SmootherPrototype> smooProto, Teuchos::FancyOStream & out, bool & success) {
      TEST_THROW(rcp( new SmootherFactory(smooProto) ), RuntimeError);
    }
  
    TEUCHOS_UNIT_TEST(SmootherFactory, Constructor_OneArg) 
    {
      testCollection(&testConstructor1, out, success);               // TEST: Valid input parameter
      testInvalidCollection(&testInvalidConstructor1, out, success); // TEST: Invalid input parameter
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Test: Constructor_TwoArg
    ////////////////////////////////////////////////////////////////////////////////

    // Build an object and test internal state (two arguments)
    void testConstructor2(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream & out, bool & success) {
      check(rcp( new SmootherFactory(smooProtoA, smooProtoB) ), smooProtoA, smooProtoB, out, success);  
    }
  
    // Test with invalid input arguments
    void testInvalidConstructor2(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream & out, bool & success) {
      TEST_THROW(rcp( new SmootherFactory(smooProtoA, smooProtoB) ), RuntimeError);
    }

    TEUCHOS_UNIT_TEST(SmootherFactory, Constructor_TwoArg) 
    {
      testCollection(&testConstructor2, out, success);               // TEST: Valid input parameter
      testInvalidCollection(&testInvalidConstructor2, out, success); // TEST: Valid input parameter
    }
  
    ////////////////////////////////////////////////////////////////////////////////
    // Test: SetSmootherPrototypes_OneArg
    ////////////////////////////////////////////////////////////////////////////////

    void testSetSmootherPrototypes1(RCP<SmootherPrototype> smooProto, Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherPrototype> smooProto0 = rcp( new FakeSmootherPrototype() );
      RCP<SmootherFactory>   smooFact   = rcp( new SmootherFactory(smooProto0) );
    
      ineqcheck(smooFact, smooProto, smooProto, out, success);
      smooFact->SetSmootherPrototypes(smooProto);
      check(smooFact, smooProto, smooProto, out, success); // One argument == same pre and post smoother prototype
    }

    void testInvalidSetSmootherPrototypes1(RCP<SmootherPrototype> smooProto, Teuchos::FancyOStream & out, bool & success) {
      SmootherFactory smooFact(Teuchos::null);
      TEST_THROW(smooFact.SetSmootherPrototypes(smooProto), RuntimeError);
    }

    TEUCHOS_UNIT_TEST(SmootherFactory, SetSmootherPrototypes_OneArg)
    {    
      testCollection(&testSetSmootherPrototypes1, out, success);               // TEST: Valid input parameter
      testInvalidCollection(&testInvalidSetSmootherPrototypes1, out, success); // TEST: Invalid input parameter
    }
  
    ////////////////////////////////////////////////////////////////////////////////
    // Test: SetSmootherPrototypes_TwoArg
    ////////////////////////////////////////////////////////////////////////////////

    void testSetSmootherPrototypes2(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherPrototype> smooProto0 = rcp( new FakeSmootherPrototype() );
      RCP<SmootherFactory>   smooFact   = rcp( new SmootherFactory(smooProto0) );
    
      ineqcheck(smooFact, smooProtoA, smooProtoB, out, success);
      smooFact->SetSmootherPrototypes(smooProtoA, smooProtoB);
      check(smooFact, smooProtoA, smooProtoB, out, success);  
    }

    void testInvalidSetSmootherPrototypes2(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream & out, bool & success) {
      SmootherFactory smooFact(Teuchos::null);
      TEST_THROW(smooFact.SetSmootherPrototypes(smooProtoA, smooProtoB), RuntimeError);
    }

    TEUCHOS_UNIT_TEST(SmootherFactory, SetSmootherPrototypes_TwoArg)
    {
      testCollection(&testSetSmootherPrototypes2, out, success);               // TEST: Valid input parameter
      testInvalidCollection(&testInvalidSetSmootherPrototypes2, out, success); // TEST: Invalid input parameter
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Test: GetSmootherPrototypes
    ////////////////////////////////////////////////////////////////////////////////

    TEUCHOS_UNIT_TEST(SmootherFactory, GetSmootherPrototypes)
    {
      // Already heavily tested by previous tests.
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Test: DeclareInput
    ////////////////////////////////////////////////////////////////////////////////

    TEUCHOS_UNIT_TEST(SmootherFactory, DeclareInput)
    {
      // TODO:
      // Create a level, call DeclareInput() and Build() and then check that every counter == 0 at the end
      // Same test have to be done for each factory...
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Test: Build
    ////////////////////////////////////////////////////////////////////////////////

    // test if a smoother created by Build() is correct (check if it corresponds to the prototype)
    void testBuildCheckOutput(RCP<SmootherFactory> smooFact, Level& level, RCP<SmootherPrototype> smooProto, const std::string& tag, Teuchos::FancyOStream & out, bool & success) {
      if (smooProto == Teuchos::null) {
        TEST_EQUALITY(level.IsAvailable(tag, smooFact), false);
      } else {
        RCP<SmootherBase> smoother;
        TEST_NOTHROW(smoother = level.Get< RCP<SmootherBase> >(tag, smooFact));

        TEST_INEQUALITY(smoother, Teuchos::null);
        TEST_INEQUALITY(smoother, smooProto);

        if (smooProto != Teuchos::null) {
          RCP<FakeSmootherPrototype> smootherF;

          // ouput test: smoothers same derived class as prototypes
          TEST_NOTHROW(smootherF = rcp_dynamic_cast<FakeSmootherPrototype>(smoother,true));
          
          if (smootherF != Teuchos::null) {
            // output test: smoother parameters == prototype parameters
            RCP<FakeSmootherPrototype> smooProtoF = rcp_dynamic_cast<FakeSmootherPrototype>(smooProto,true);
            TEST_EQUALITY(smootherF->GetParam(), smooProtoF->GetParam());
            
            // output test: smoothers are ready to be apply
            TEST_EQUALITY(smootherF->IsSetup(), true);
            
            // setup done only once.
            TEST_EQUALITY(smootherF->GetNumOfSetupCall(), 1);
            TEST_EQUALITY(smootherF->GetNumOfSetup(), 1);
          }
        }
      }
    }

    // test outpout of Build() and BuildSmoother()
    void testBuildCheck(RCP<SmootherFactory> smooFact, Level& level, RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, MueLu::PreOrPost preOrPost, Teuchos::FancyOStream & out, bool & success) {
      // invariant: smoother prototypes kept unchanged
      check(smooFact, smooProtoA, smooProtoB, out, success);
      
      // output test
      if (preOrPost == MueLu::BOTH) {

        testBuildCheckOutput(smooFact, level, smooProtoA, "PreSmoother", out, success);
        testBuildCheckOutput(smooFact, level, smooProtoB, "PostSmoother", out, success);
        
        // ReUse: if pre and post prototype are the same, then pre smoother == post smoother
        // otherwise, they are different (have not been tested by previous tests)
        RCP<SmootherBase> smooA, smooB;
        if (smooProtoA != Teuchos::null) { smooA = level.Get< RCP<SmootherBase> >("PreSmoother",  smooFact); }
        if (smooProtoB != Teuchos::null) { smooB = level.Get< RCP<SmootherBase> >("PostSmoother", smooFact); }
        if (smooProtoA == smooProtoB) { TEST_EQUALITY(smooA, smooB); } else { TEST_INEQUALITY(smooA, smooB); }

      } else if (preOrPost == MueLu::PRE) {

        testBuildCheckOutput(smooFact, level, smooProtoA, "PreSmoother", out, success);
        TEST_EQUALITY(level.IsAvailable("PostSmoother", smooFact), false);

      } else if (preOrPost == MueLu::POST) {

        TEST_EQUALITY(level.IsAvailable("PreSmoother", smooFact), false);
        testBuildCheckOutput(smooFact, level, smooProtoB, "PostSmoother", out, success);

      } else { TEST_EQUALITY(true, false); }

    }

    // 
    void testBuild(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherFactory> smooFact = rcp( new SmootherFactory(smooProtoA, smooProtoB) );

      Level level; //level.SetupPhase(true);
      TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSingleLevelHierarchy(level);

      level.Request("PreSmoother",smooFact.get());
      level.Request("PostSmoother", smooFact.get());

      smooFact->Build(level);

      testBuildCheck(smooFact, level, smooProtoA, smooProtoB, MueLu::BOTH, out, success);
    }

    TEUCHOS_UNIT_TEST(SmootherFactory, Build)
    {
      testCollection(&testBuild, out, success);
    }
  
    ////////////////////////////////////////////////////////////////////////////////
    // Test: BuildSmoother
    ////////////////////////////////////////////////////////////////////////////////

    void testBuildSmootherPreOrPost(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, MueLu::PreOrPost preOrPost, Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherFactory> smooFact = rcp( new SmootherFactory(smooProtoA, smooProtoB) );
      
      Level level; //level.SetupPhase(true);
      TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSingleLevelHierarchy(level);

      level.Request("PreSmoother",smooFact.get());
      level.Request("PostSmoother", smooFact.get());

      smooFact->BuildSmoother(level, preOrPost);

      testBuildCheck(smooFact, level, smooProtoA, smooProtoB, preOrPost, out, success);
    }

    void testBuildSmoother(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream & out, bool & success) {
      testBuildSmootherPreOrPost(smooProtoA, smooProtoB, MueLu::PRE,  out, success);
      testBuildSmootherPreOrPost(smooProtoA, smooProtoB, MueLu::POST, out, success);
      testBuildSmootherPreOrPost(smooProtoA, smooProtoB, MueLu::BOTH, out, success);
    }
    
    void testBuildSmootherDefaultArg(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherFactory> smooFact = rcp( new SmootherFactory(smooProtoA, smooProtoB) );
      
      Level level; //level.SetupPhase(true);
      TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSingleLevelHierarchy(level);

      level.Request("PreSmoother",smooFact.get());
      level.Request("PostSmoother", smooFact.get());

      smooFact->BuildSmoother(level);

      testBuildCheck(smooFact, level, smooProtoA, smooProtoB, MueLu::BOTH, out, success);
    }

    TEUCHOS_UNIT_TEST(SmootherFactory, BuildSmoother)
    {
      testCollection(&testBuildSmootherDefaultArg, out, success); // TEST: default arg
      testCollection(&testBuildSmoother, out, success);           // TEST: PRE, POST and BOTH
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Test: Compilation
    ////////////////////////////////////////////////////////////////////////////////

    // Make sure that Constructors and SetSmootherPrototypes() are not using references as input arguments
    // (to allow Teuchos::null as input arguments)
    TEUCHOS_UNIT_TEST(SmootherFactory, Compilation)
    {
      SmootherFactory smooFact1(Teuchos::null);
      SmootherFactory smooFact2(Teuchos::null, Teuchos::null);
      smooFact1.SetSmootherPrototypes(Teuchos::null);
      smooFact1.SetSmootherPrototypes(Teuchos::null, Teuchos::null);
    }

  } // namespace SmootherFactory

} // namespace MueLuTests

//TODO: advance reuse test
//TODO: remove some useless RCP in helper function prototypes?
