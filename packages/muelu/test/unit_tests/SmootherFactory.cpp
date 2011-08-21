#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_Version.hpp"

#include "MueLu_TestHelpers.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_IfpackSmoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

//TODO from JG: should be tested using a fakeSmoother

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
    //we are now in a class method declared by the above macro, and
    //that method has these input arguments:
    //Teuchos::FancyOStream& out, bool& success
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
      {
 
        out << "version: " << MueLu::Version() << std::endl;
#ifdef HAVE_MUELU_IFPACK
        Teuchos::ParameterList  ifpackList;
        ifpackList.set("relaxation: type", "Gauss-Seidel");
        ifpackList.set("relaxation: sweeps", (LO) 1);
        ifpackList.set("relaxation: damping factor", (SC) 1.0);
        RCP<SmootherPrototype>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
        RCP<SmootherFactory> smooFact = rcp(new SmootherFactory(smoother) );
        TEUCHOS_TEST_EQUALITY(smooFact != Teuchos::null, true, out, success);
#endif
      }
  }

  TEUCHOS_UNIT_TEST(SmootherFactory, DefaultCtor_TwoArgs)
  {
    //we are now in a class method declared by the above macro, and
    //that method has these input arguments:
    //Teuchos::FancyOStream& out, bool& success
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
      {
 
        out << "version: " << MueLu::Version() << std::endl;
#ifdef HAVE_MUELU_IFPACK
        Teuchos::ParameterList  ifpackList;
        ifpackList.set("relaxation: type", "Gauss-Seidel");
        ifpackList.set("relaxation: sweeps", (LO) 1);
        ifpackList.set("relaxation: damping factor", (SC) 1.0);
        RCP<SmootherPrototype>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
        RCP<SmootherFactory> smooFact = rcp(new SmootherFactory(smoother,smoother) );
        TEUCHOS_TEST_EQUALITY(smooFact != Teuchos::null, true, out, success);
#endif
      }
  }

  TEUCHOS_UNIT_TEST(SmootherFactory, SetSmootherPrototypes)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
      {
 
        out << "version: " << MueLu::Version() << std::endl;
#ifdef HAVE_MUELU_IFPACK
        Teuchos::ParameterList  ifpackList;
        ifpackList.set("relaxation: type", "Gauss-Seidel");
        ifpackList.set("relaxation: sweeps", (LO) 1);
        ifpackList.set("relaxation: damping factor", (SC) 1.0);
        RCP<SmootherPrototype>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
        RCP<SmootherFactory> smooFact = rcp(new SmootherFactory(smoother) );

        RCP<SmootherPrototype>  newSmoo1 = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
        RCP<SmootherPrototype>  newSmoo2 = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
        smooFact->SetSmootherPrototypes(newSmoo1,newSmoo2);
        RCP<SmootherPrototype>  checkSmoo1;
        RCP<SmootherPrototype>  checkSmoo2;
        smooFact->GetSmootherPrototypes(checkSmoo1,checkSmoo2);

        TEUCHOS_TEST_EQUALITY(checkSmoo1 == newSmoo1, true, out, success);
        TEUCHOS_TEST_EQUALITY(checkSmoo2 == newSmoo2, true, out, success);
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
        Teuchos::ParameterList  ifpackList1;
        ifpackList1.set("relaxation: type", "Gauss-Seidel");
        ifpackList1.set("relaxation: sweeps", (LO) 1);
        ifpackList1.set("relaxation: damping factor", (SC) 1.0);
        RCP<SmootherPrototype>  smooProto1 = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList1) );

        //post-smoother prototype
        Teuchos::ParameterList  ifpackList2;
        ifpackList2.set("relaxation: type", "Jacobi");
        ifpackList2.set("relaxation: sweeps", (LO) 1);
        ifpackList2.set("relaxation: damping factor", (SC) 1.0);
        RCP<SmootherPrototype>  smooProto2 = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList2) );
        //RCP<SmootherFactory> smooFactory = rcp(new SmootherFactory(smooProto1,smooProto2) );
        RCP<SmootherFactory> smooFactory = rcp(new SmootherFactory(smooProto1) );

        RCP<Level> aLevel = rcp(new Level() );
        aLevel->SetLevelID(1);
        RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99);

        RCP<SmootherPrototype>  preSmoo, postSmoo;
        //Check for exception if matrix is not set in Level.
        //FIXME once Level-template Get< RCP<Operator> >("A") doesn't throw an exception, this must be changed
        //TEST_THROW(smooFactory->Build(aLevel,preSmoo,postSmoo),MueLu::Exceptions::RuntimeError);
        TEST_THROW(smooFactory->Build(aLevel,preSmoo,postSmoo),std::logic_error);

        aLevel->Set("A",A);

        //same prototype for pre and post smoothers
        smooFactory->Build(aLevel,preSmoo,postSmoo);
        TEUCHOS_TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel",out,success);
        TEUCHOS_TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Gauss-Seidel",out,success);

        //different prototypes for pre and post smoothers
        smooFactory = rcp(new SmootherFactory(smooProto1,smooProto2) );
        smooFactory->Build(aLevel,preSmoo,postSmoo);
        TEUCHOS_TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel",out,success);
        TEUCHOS_TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Jacobi",out,success);
#endif
      }
  }

}//namespace MueLuTests

