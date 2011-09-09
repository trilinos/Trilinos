#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

/*
  Comments:
    1) Chebyshev smoothing must pass for any number of processors.
    2) Gauss-Seidel must pass for 1 and 4 processors.
    3) For any processor count except 1 and 4, the Gauss-Seidel test will
       report "passing", but this is only because the Teuchos test macro is skipped.
*/

namespace MueLuTests {

  // TODO added by JG: This test should be factorized for every smoothers, to avoid code duplication
  TEUCHOS_UNIT_TEST(IfpackSmoother, NotSetup)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)
      {
        out << "version: " << MueLu::Version() << std::endl;

        GO nx,ny,nz;
        nx = ny = nz = 5;
        GO numGlobalElements = nx*ny*nz;
        RCP<const Map> map = MapFactory::Build(Parameters::getLib(), numGlobalElements, 0, Parameters::getDefaultComm());

        Teuchos::ParameterList  ifpackList;
        RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );

        RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
        RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

        //try applying without setting up
        TEST_THROW( smoother->Apply(*X,*RHS) , MueLu::Exceptions::RuntimeError );
        TEST_THROW( smoother->SetNIts(5), MueLu::Exceptions::RuntimeError );
      }
  }

  TEUCHOS_UNIT_TEST(IfpackSmoother, GaussSeidel)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)
    {
      out << "version: " << MueLu::Version() << std::endl;
      out << "Tests interface to Ifpack's Gauss-Seidel preconditioner." << std::endl;
      RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();

      Teuchos::ParameterList  ifpackList;
      ifpackList.set("relaxation: type", "Gauss-Seidel");
      ifpackList.set("relaxation: sweeps", (int) 1);
      ifpackList.set("relaxation: damping factor", (double) 1.0);
      ifpackList.set("relaxation: zero starting solution", false);
      RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
      Level aLevel;
      MueLu::TestHelpers::Factory<SC,LO,GO,NO,LMO>::createSingleLevelHierarchy(aLevel);

      RCP<Operator> Op = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(125);
      aLevel.Set("A",Op);

      RCP<MultiVector> X = MultiVectorFactory::Build(Op->getDomainMap(),1);
      RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getRangeMap(),1);

      smoother->Setup(aLevel);

      Teuchos::Array<ST::magnitudeType> norms(1);

      X->putScalar( (SC) 1.0);
      RHS->putScalar( (SC) 0.0);

      out << "Applying one GS sweep" << std::endl;
      smoother->Apply(*X,*RHS);
      norms = Utils::ResidualNorm(*Op,*X,*RHS);
      switch (comm->getSize()) {
        case 1:
        case 4:
          TEST_FLOATING_EQUALITY(norms[0],5.773502691896257e-01,1e-12);
          break;
        default:
          out << "Pass/Fail is checked only for 1 and 4 processes." << std::endl;
          break;
      } //switch

      int numIts = 10;
      smoother->SetNIts(numIts);
      TEST_EQUALITY(smoother->GetNIts(),10);
      out << "Applying " << numIts << " GS sweeps" << std::endl;
      X->putScalar( (SC) 1.0);
      smoother->Apply(*X,*RHS);
      norms = Utils::ResidualNorm(*Op,*X,*RHS);
      switch (comm->getSize()) {
        case 1:
          TEST_FLOATING_EQUALITY(norms[0],8.326553652741774e-02,1e-12);
          break;
        case 4:
          TEST_FLOATING_EQUALITY(norms[0],8.326553653078517e-02,1e-12);
          break;
        default:
          out << "Pass/Fail is checked only for 1 and 4 processes." << std::endl;
          break;
      } //switch
    }
  } //GaussSeidel

  TEUCHOS_UNIT_TEST(IfpackSmoother, Chebyshev)
  {
    //TODO This test should really calculate the reduction analytically.

    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)
      {
        out << "version: " << MueLu::Version() << std::endl;

        out << "Tests interface to Ifpack's Chebyshev preconditioner." << std::endl;

        Teuchos::ParameterList  ifpackList;
        ifpackList.set("chebyshev: degree", (int) 3);
        ifpackList.set("chebyshev: max eigenvalue", (double) 1.98476);
        ifpackList.set("chebyshev: min eigenvalue", (double) 1.0);
        ifpackList.set("chebyshev: ratio eigenvalue", (double) 20);
        ifpackList.set("chebyshev: zero starting solution", false);
        RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("Chebyshev",ifpackList) );
        
        Level aLevel;
        MueLu::TestHelpers::Factory<SC,LO,GO,NO,LMO>::createSingleLevelHierarchy(aLevel);

        RCP<Operator> Op = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(125);
        aLevel.Set("A",Op);

        RCP<MultiVector> X = MultiVectorFactory::Build(Op->getDomainMap(),1);
        RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getRangeMap(),1);
        RCP<MultiVector> RES = MultiVectorFactory::Build(Op->getRangeMap(),1);
        RHS->putScalar( (SC) 0.0);
        X->putScalar(1.0);

        smoother->Setup(aLevel);

        Teuchos::Array<ST::magnitudeType> norms(1);
        X->norm2(norms);
        out << "||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

        int numIts = 3;
        smoother->SetNIts(numIts);
        TEST_EQUALITY(smoother->GetNIts(),numIts);
        out << "Applying degree " << numIts << " Chebyshev smoother" << std::endl;
        X->putScalar(1.0);
        out << "||X_initial|| = " << std::setiosflags(ios::fixed) << std::setprecision(25) << norms[0] << std::endl;
        smoother->Apply(*X,*RHS);
        Op->apply(*X,*RES,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
        RES->norm2(norms);
        out << "||RES|| = " << std::setiosflags(ios::fixed) <<
          std::setprecision(20) << norms[0] << std::endl;
        TEST_FLOATING_EQUALITY(norms[0],5.269156e-01,1e-7);  //Compare to residual reported by ML

      }
  } //Chebyshev

  TEUCHOS_UNIT_TEST(IfpackSmoother, ILU)
  {

    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)
      {
        //FIXME this will probably fail in parallel b/c it becomes block Jacobi

        out << "version: " << MueLu::Version() << std::endl;

        out << "Tests interface to Ifpack's ILU(0) preconditioner." << std::endl;

        RCP<const Teuchos::Comm<int> > comm = MueLu::TestHelpers::Parameters::getDefaultComm();

        Teuchos::ParameterList  ifpackList;
        RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("ILU",ifpackList) );

        Level aLevel;
        MueLu::TestHelpers::Factory<SC,LO,GO,NO,LMO>::createSingleLevelHierarchy(aLevel);

        RCP<Operator> Op = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(125);
        aLevel.Set("A",Op);

        RCP<MultiVector> Xtrue = MultiVectorFactory::Build(Op->getDomainMap(),1);
        RCP<MultiVector> X = MultiVectorFactory::Build(Op->getDomainMap(),1);
        RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getDomainMap(),1);

        smoother->Setup(aLevel);
        Xtrue->setSeed(846930886);
        Xtrue->randomize();
        Teuchos::Array<ST::magnitudeType> norms(1);
        Xtrue->norm2(norms);
        Xtrue->scale(1/norms[0]);
        Op->apply(*Xtrue,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
        RCP<Epetra_MultiVector> epRHS = Utils::MV2NonConstEpetraMV(RHS);
        RHS->norm2(norms);
        out << "||RHS|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
        X->putScalar( (SC) 0.0);


        smoother->Apply(*X,*RHS);
        X->norm2(norms);
        out << "||X_final|| = " << std::setiosflags(ios::fixed) << std::setprecision(25) << norms[0] << std::endl;
        norms = Utils::ResidualNorm(*Op,*X,*RHS);
        out << "||residual|| = " << norms[0] << std::endl;
        if (comm->getSize() == 1) {
          TEST_EQUALITY(norms[0]<1e-10,true);
        } else {
          out << "Pass/Fail is only checked in serial." << std::endl;
        }
        
      }
  } //ILU

}//namespace MueLuTests
