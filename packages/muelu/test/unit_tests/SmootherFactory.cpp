#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_TestHelpers.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_FakeSmootherPrototype.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  namespace SmootherFactoryHelpers {
    // This helper function checks the consistency of the internal state of a SmootherFactory object
    void check(RCP<const SmootherFactory> smootherFactory, RCP<const SmootherPrototype> preSmooProto, RCP<const SmootherPrototype> postSmooProto, Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherPrototype> getPreSmooProto, getPostSmooProto;
      smootherFactory->GetSmootherPrototypes(getPreSmooProto, getPostSmooProto);
      TEST_EQUALITY(preSmooProto, getPreSmooProto);
      TEST_EQUALITY(postSmooProto, getPostSmooProto);
    }

    void ineqcheck(RCP<const SmootherFactory> smootherFactory, RCP<const SmootherPrototype> preSmooProto, RCP<const SmootherPrototype> postSmooProto, Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherPrototype> getPreSmooProto, getPostSmooProto;
      smootherFactory->GetSmootherPrototypes(getPreSmooProto, getPostSmooProto);
      TEST_INEQUALITY(preSmooProto, getPreSmooProto);
      TEST_INEQUALITY(postSmooProto, getPostSmooProto);
    }
    
    void testConstructor1(RCP<SmootherPrototype> smooProto, Teuchos::FancyOStream & out, bool & success) {
      check(rcp( new SmootherFactory(smooProto) ), smooProto, smooProto, out, success);  
    }
    
    void testConstructor2(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream & out, bool & success) {
      check(rcp( new SmootherFactory(smooProtoA, smooProtoB) ), smooProtoA, smooProtoB, out, success);  
    }
    
   void testSet1(RCP<SmootherPrototype> smooProto, Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherPrototype> smooProto0 = rcp( new FakeSmootherPrototype() );
      RCP<SmootherFactory>   smooFact  = rcp( new SmootherFactory(smooProto0) );

      ineqcheck(smooFact, smooProto, smooProto, out, success);
      smooFact->SetSmootherPrototypes(smooProto);
      check(smooFact, smooProto, smooProto, out, success);  
    }

    void testSet2(RCP<SmootherPrototype> smooProtoA, RCP<SmootherPrototype> smooProtoB, Teuchos::FancyOStream & out, bool & success) {
      RCP<SmootherPrototype> smooProto0 = rcp( new FakeSmootherPrototype() );
      RCP<SmootherFactory>   smooFact  = rcp( new SmootherFactory(smooProto0) );

      ineqcheck(smooFact, smooProtoA, smooProtoB, out, success);
      smooFact->SetSmootherPrototypes(smooProtoA, smooProtoB);
      check(smooFact, smooProtoA, smooProtoB, out, success);  
    }
  }

#define testConstructor1(proto)          SmootherFactoryHelpers::testConstructor1(proto, out, success)
#define testConstructor2(protoA, protoB) SmootherFactoryHelpers::testConstructor2(protoA, protoB, out, success)
#define testSet1(proto)                  SmootherFactoryHelpers::testSet1(proto, out, success)
#define testSet2(protoA, protoB)         SmootherFactoryHelpers::testSet2(protoA, protoB, out, success)

  TEUCHOS_UNIT_TEST(SmootherFactory, Constructor) 
  {
    RCP<SmootherPrototype> smooProto = rcp( new FakeSmootherPrototype() );

    // TEST: Teuchos::null is a valid argument
    testConstructor1(Teuchos::null);

    // TEST: One argument == same pre and post smoother prototype
    testConstructor1(smooProto);
  }

  TEUCHOS_UNIT_TEST(SmootherFactory, Constructor_OneArg) 
  {
    RCP<SmootherPrototype> smooProtoA = rcp( new FakeSmootherPrototype() );
    RCP<SmootherPrototype> smooProtoB = rcp( new FakeSmootherPrototype() );

    testConstructor2(Teuchos::null, Teuchos::null);  
    testConstructor2(smooProtoA, smooProtoA);
    testConstructor2(smooProtoA, smooProtoB);
    testConstructor2(smooProtoA, Teuchos::null);
    testConstructor2(Teuchos::null, smooProtoB);
  }
  
  TEUCHOS_UNIT_TEST(SmootherFactory, SetSmootherPrototypes_OneArg)
  {
    
    RCP<SmootherPrototype> smooProto = rcp( new FakeSmootherPrototype() );
    
    testSet1(Teuchos::null);  
    testSet1(smooProto);
  }
  
  TEUCHOS_UNIT_TEST(SmootherFactory, SetSmootherPrototypes_TwoArg)
  {
    
    RCP<SmootherPrototype> smooProtoA = rcp( new FakeSmootherPrototype() );
    RCP<SmootherPrototype> smooProtoB = rcp( new FakeSmootherPrototype() );
    
    testSet2(Teuchos::null, Teuchos::null);  
    testSet2(smooProtoA, smooProtoA);
    testSet2(smooProtoA, smooProtoB);
    testSet2(smooProtoA, Teuchos::null);
    testSet2(Teuchos::null, smooProtoB);
  }

  TEUCHOS_UNIT_TEST(SmootherFactory, GetSmootherPrototypes)
  {
    // Already tested by previous tests.
  }

  TEUCHOS_UNIT_TEST(SmootherFactory, Build)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
      {
        out << "version: " << MueLu::Version() << std::endl;
        out << "Testing SmootherFactory::Build method" << std::endl;
#ifdef HAVE_MUELU_IFPACK
        //pre-smoother prototype
        RCP<SmootherPrototype> smooProto1 = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel");

        //post-smoother prototype
        RCP<SmootherPrototype>  smooProto2 = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Jacobi");

        RCP<SmootherFactory> smootherFactoryory = rcp(new SmootherFactory(smooProto1) );

        RCP<Level> aLevel = rcp(new Level() );
        MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSingleLevelHierarchy(*aLevel);
        RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99);

        //Check for exception if matrix is not set in Level.
        //FIXME once Level-template Get< RCP<Operator> >("A") doesn't throw an exception, this must be changed
        //TEST_THROW(smootherFactoryory->Build(aLevel,preSmoo,postSmoo),MueLu::Exceptions::RuntimeError);
        TEST_THROW(smootherFactoryory->Build(*aLevel),std::logic_error);

        aLevel->Set("A",A);

        //same prototype for pre and post smoothers
        smootherFactoryory->Build(*aLevel);
        {
          RCP<SmootherBase> preSmoo = aLevel->Get< RCP<SmootherBase> >("PreSmoother", smootherFactoryory);
          RCP<SmootherBase> postSmoo = aLevel->Get< RCP<SmootherBase> >("PostSmoother", smootherFactoryory);
          //JGTODO          TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel");
          //JGTODO          TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Gauss-Seidel");
        }

        //different prototypes for pre and post smoothers
        smootherFactoryory = rcp(new SmootherFactory(smooProto1,smooProto2) );
        smootherFactoryory->Build(*aLevel);
        {
          RCP<SmootherBase> preSmoo = aLevel->Get< RCP<SmootherBase> >("PreSmoother", smootherFactoryory);
          RCP<SmootherBase> postSmoo = aLevel->Get< RCP<SmootherBase> >("PostSmoother", smootherFactoryory);
          //JGTODO TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel");
          //JGTODO TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Jacobi");
        }
#endif
      }
  }

#undef testConstructor1
#undef testConstructor2
#undef testSet1
#undef testSet2

} // namespace MueLuTests
