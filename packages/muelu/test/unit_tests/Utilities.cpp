#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Utilities.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

// This file is intended to house all the tests for MueLu_Utilities.hpp.

namespace MueLuTests {
  
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_TPETRA)
  TEUCHOS_UNIT_TEST(Utilities,MatMatMult_EpetraVsTpetra)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "This test compares the matrix matrix multiply between Tpetra and Epetra" << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Teuchos::Array<ST::magnitudeType> normResult1(1);

    //Calculate result = (Op*Op)*X for Tpetra
    int nx = 37*comm->getSize();
    int ny=nx;
    RCP<Operator> Op = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build2DPoisson(nx,ny,Xpetra::UseEpetra);
    RCP<Operator> OpOp = Utils::TwoMatrixMultiply(Op,false,Op,false);
    RCP<MultiVector> result = MultiVectorFactory::Build(OpOp->getRangeMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(OpOp->getDomainMap(),1);
    Teuchos::Array<ST::magnitudeType> xnorm(1);
    X->setSeed(8675309);
    X->randomize(true);
    X->norm2(xnorm);
    OpOp->apply(*X,*result,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    result->norm2(normResult1);

    //Calculate result = (Op*Op)*X for Epetra
    Op = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build2DPoisson(nx,ny,Xpetra::UseTpetra);
    OpOp = Utils::TwoMatrixMultiply(Op,false,Op,false);
    result = MultiVectorFactory::Build(OpOp->getRangeMap(),1);
    X = MultiVectorFactory::Build(OpOp->getDomainMap(),1);
    X->setSeed(8675309);
    X->randomize(true);
    X->norm2(xnorm);
    OpOp->apply(*X,*result,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    Teuchos::Array<ST::magnitudeType> normResult2(1);
    result->norm2(normResult2);

    TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12);

  } //EpetraVersusTpetra
#endif

}//namespace MueLuTests

