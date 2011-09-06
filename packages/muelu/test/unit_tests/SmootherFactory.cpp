#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_Version.hpp"

#include "MueLu_TestHelpers.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  TEUCHOS_UNIT_TEST(SmootherFactory, DefaultCtor_Exception)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<SmootherPrototype> smoother = Teuchos::null;
    RCP<SmootherFactory> smooFact;

    TEST_THROW( smooFact = rcp(new SmootherFactory(smoother) ) , MueLu::Exceptions::RuntimeError );

  }

  TEUCHOS_UNIT_TEST(SmootherFactory, DefaultCtor_OneArg)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
      {
 
        out << "version: " << MueLu::Version() << std::endl;
#ifdef HAVE_MUELU_IFPACK
        RCP<SmootherPrototype> smoother = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype();
        RCP<SmootherFactory> smooFact = rcp(new SmootherFactory(smoother) );
        TEST_EQUALITY(smooFact != Teuchos::null, true);
#endif
      }
  }

  TEUCHOS_UNIT_TEST(SmootherFactory, DefaultCtor_TwoArgs)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
      {
 
        out << "version: " << MueLu::Version() << std::endl;
#ifdef HAVE_MUELU_IFPACK
        RCP<SmootherPrototype> smoother = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype();
        RCP<SmootherFactory> smooFact = rcp(new SmootherFactory(smoother,smoother) );
        TEST_EQUALITY(smooFact != Teuchos::null, true);
#endif
      }
  }

  TEUCHOS_UNIT_TEST(SmootherFactory, SetSmootherPrototypes)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
      {
 
        out << "version: " << MueLu::Version() << std::endl;
#ifdef HAVE_MUELU_IFPACK
        RCP<SmootherPrototype> smoother = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype();
        RCP<SmootherFactory> smooFact = rcp(new SmootherFactory(smoother) );

        RCP<SmootherPrototype>  newSmoo1 = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype();
        RCP<SmootherPrototype>  newSmoo2 = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype();
        smooFact->SetSmootherPrototypes(newSmoo1,newSmoo2);
        RCP<SmootherPrototype>  checkSmoo1;
        RCP<SmootherPrototype>  checkSmoo2;
        smooFact->GetSmootherPrototypes(checkSmoo1,checkSmoo2);

        TEST_EQUALITY(checkSmoo1 == newSmoo1, true);
        TEST_EQUALITY(checkSmoo2 == newSmoo2, true);
#endif
      }
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

        RCP<SmootherFactory> smooFactory = rcp(new SmootherFactory(smooProto1) );

        RCP<Level> aLevel = rcp(new Level() );
        MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSingleLevelHierarchy(*aLevel);
        RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99);

        //Check for exception if matrix is not set in Level.
        //FIXME once Level-template Get< RCP<Operator> >("A") doesn't throw an exception, this must be changed
        //TEST_THROW(smooFactory->Build(aLevel,preSmoo,postSmoo),MueLu::Exceptions::RuntimeError);
        TEST_THROW(smooFactory->Build(*aLevel),std::logic_error);

        aLevel->Set("A",A);

        //same prototype for pre and post smoothers
        smooFactory->Build(*aLevel);
        {
          RCP<SmootherBase> preSmoo = aLevel->Get< RCP<SmootherBase> >("PreSmoother", smooFactory);
          RCP<SmootherBase> postSmoo = aLevel->Get< RCP<SmootherBase> >("PostSmoother", smooFactory);
          //JGTODO          TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel");
          //JGTODO          TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Gauss-Seidel");
        }

        //different prototypes for pre and post smoothers
        smooFactory = rcp(new SmootherFactory(smooProto1,smooProto2) );
        smooFactory->Build(*aLevel);
        {
          RCP<SmootherBase> preSmoo = aLevel->Get< RCP<SmootherBase> >("PreSmoother", smooFactory);
          RCP<SmootherBase> postSmoo = aLevel->Get< RCP<SmootherBase> >("PostSmoother", smooFactory);
          //JGTODO TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel");
          //JGTODO TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Jacobi");
        }
#endif
      }
  }

}//namespace MueLuTests

