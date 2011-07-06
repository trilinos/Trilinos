#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace MueLu::TestHelpers;

  //this macro declares the unit-test-class:
  TEUCHOS_UNIT_TEST(Amesos, NotSetup)
  {
    //we are now in a class method declared by the above macro, and
    //that method has these input arguments:
    //Teuchos::FancyOStream& out, bool& success

    MUELU_TEST_ONLY_FOR(Cthulhu::UseEpetra)
      {
        out << "version: " << MueLu::Version() << std::endl;

        GO nx,ny,nz;
        nx = ny = nz = 5;
        GO numGlobalElements = nx*ny*nz;
        RCP<const Map> map = MapFactory::Build(Parameters::getLib(), numGlobalElements, 0, Parameters::getDefaultComm());

        Teuchos::ParameterList amesosList;
        RCP<AmesosSmoother>    smoother = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );

        RCP<MultiVector> X   = MultiVectorFactory::Build(map,1);
        RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

        //try applying without setting up
        TEST_THROW( smoother->Apply(*X,*RHS), MueLu::Exceptions::RuntimeError );
        TEST_THROW( smoother->SetNIts(5),     MueLu::Exceptions::RuntimeError ); //JG: Why ?
      }
  }

  TEUCHOS_UNIT_TEST(Amesos, KLUSolve)
  {
    MUELU_TEST_ONLY_FOR(Cthulhu::UseEpetra)
      {
        out << "version: " << MueLu::Version() << std::endl;

        Teuchos::ParameterList amesosList;
        amesosList.set("PrintTiming",false);
        amesosList.set("OutputLevel",0);
        RCP<AmesosSmoother> smoother = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );

        RCP<Operator> Op = MueLu::TestHelpers::Factory<SC,LO,GO,NO,LMO>::Build1DPoisson(125);
        Level aLevel;
        aLevel.SetA(Op);

        smoother->Setup(aLevel);

        RCP<MultiVector> X   = MultiVectorFactory::Build(Op->getDomainMap(),1);
        RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getDomainMap(),1);

        X->setSeed(846930886);
        X->randomize();
        Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

        X->putScalar((SC) 0.0);
        Teuchos::Array<ST::magnitudeType> res;
        res = Utils::ResidualNorm(*Op,*X,*RHS);
        out << "||initial residual|| = " << res[0] << std::endl;

        smoother->Apply(*X,*RHS);
        res = Utils::ResidualNorm(*Op,*X,*RHS);
        out << "||final residual|| = " << res[0] << std::endl;
        TEUCHOS_TEST_EQUALITY(res[0] < 1e-12,true,out,success);
      }
  }
  
} // namespace <anonymous>
