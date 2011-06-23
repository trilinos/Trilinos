#include "Teuchos_UnitTestHarness.hpp"
#include "Cthulhu_DefaultPlatform.hpp"

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Utilities.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

  using Teuchos::rcp;
  using Teuchos::RCP;
  using namespace MueLu::TestHelpers;
  
  TEUCHOS_UNIT_TEST(Transpose, Correctness)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    RCP<Operator> Op = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(27*comm->getSize());
    Level fineLevel;
    fineLevel.SetA(Op);
    Level coarseLevel;

    SaPFactory sapFactory;
    sapFactory.BuildP(fineLevel,coarseLevel);

    RCP<Operator> P = coarseLevel.GetP();
    RCP<Operator> A = fineLevel.GetA();

    RCP<Operator> R = Utils::Transpose(P);

    RCP<MultiVector> result1 = MultiVectorFactory::Build(P->getDomainMap(),1);
    RCP<MultiVector> result2  = MultiVectorFactory::Build(R->getRangeMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(P->getRangeMap(),1);
    X->randomize();
    std::string filename="P.dat";
    Utils::Write(filename,P);
    filename="R.dat";
    Utils::Write(filename,R);

    out.precision(12);
    out.setOutputToRootOnly(-1);
    X->describe(out,Teuchos::VERB_EXTREME);

    //Calculate P^T * X
    P->apply(*X,*result1,Teuchos::TRANS,(SC)1.0,(SC)0.0);
    //Calculate R * X
    R->apply(*X,*result2,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

    std::cout << "P'*x" << std::endl;
    result1->describe(out,Teuchos::VERB_EXTREME);
    std::cout << "R*x" << std::endl;
    result2->describe(out,Teuchos::VERB_EXTREME);

    Teuchos::Array<ST::magnitudeType> normX(1), normResult1(1),normResult2(1);
    X->norm2(normX);
    out << "This test checks the correctness of the Transpose utility" << std::endl;
    out << "||X||_2 = " << normX << std::endl; 
    result1->norm2(normResult1);
    result2->norm2(normResult2);
    TEUCHOS_TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12, out, success);
  } //Correctness test

}//namespace <anonymous>

