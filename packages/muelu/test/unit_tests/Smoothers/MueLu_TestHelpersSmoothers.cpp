#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_TestHelpersSmoothers.hpp"

#include "MueLu_Exceptions.hpp"

#include "MueLu_UseDefaultTypes.hpp"

// See .hpp file for description

namespace MueLuTests {

  namespace Smoother {
#include "MueLu_UseShortNames.hpp"

    // SmootherPrototype test
    void testApplyNoSetup(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success) {
      GO numGlobalElements = 125;
      RCP<const Map> map = MapFactory::Build(Parameters::getLib(), numGlobalElements, 0, Parameters::getDefaultComm());
      
      RCP<MultiVector> X   = MultiVectorFactory::Build(map,1);
      RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);
      
      TEST_THROW(smoother.Apply(*X,*RHS), MueLu::Exceptions::RuntimeError);
    }

    // SmootherBase test
    ST::magnitudeType testApply(const Operator& A, const SmootherBase & smoother, MultiVector & X, const MultiVector & RHS, Teuchos::FancyOStream & out, bool & success) {
      Array<ST::magnitudeType> norms(1);

      RHS.norm2(norms);
      out << "||RHS|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

      Teuchos::Array<ST::magnitudeType> initialNorms(1); X.norm2(initialNorms);
      out << "||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

      smoother.Apply(X, RHS); // TODO: bool const &InitialGuessIsZero=false

      Teuchos::Array<ST::magnitudeType> finalNorms(1); X.norm2(finalNorms);
      out << "||X_final|| = " << std::setiosflags(ios::fixed) << std::setprecision(25) << norms[0] << std::endl;

      Teuchos::Array<ST::magnitudeType> residualNorms = Utils::ResidualNorm(A, X, RHS);
      out << "||Residual|| = " << std::setiosflags(ios::fixed) << std::setprecision(20) << residualNorms[0] << std::endl;

      return residualNorms[0];
    }

    // SmootherBase test
    ST::magnitudeType testApply_X1_RHS0(const Operator& A, const SmootherBase & smoother, Teuchos::FancyOStream & out, bool & success) {
      RCP<MultiVector> X   = MultiVectorFactory::Build(A.getDomainMap(),1);
      RCP<MultiVector> RHS = MultiVectorFactory::Build(A.getRangeMap(),1);
      X->putScalar((SC) 1.0);
      RHS->putScalar((SC) 0.0);

      return testApply(A, smoother, *X, *RHS, out, success);
    }
    
    // SmootherBase test
    ST::magnitudeType testApply_X0_RandomRHS(const Operator& A, const SmootherBase & smoother, Teuchos::FancyOStream & out, bool & success) {
      RCP<MultiVector> X   = MultiVectorFactory::Build(A.getDomainMap(),1);
      RCP<MultiVector> RHS = MultiVectorFactory::Build(A.getRangeMap(),1);

      // Random X
      X->setSeed(846930886);
      X->randomize();
      
      // Normalize X
      Array<ST::magnitudeType> norms(1); X->norm2(norms);
      X->scale(1/norms[0]);
      
      // Compute RHS corresponding to X
      A.apply(*X,*RHS, Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
      
      // Reset X to 0
      X->putScalar((SC) 0.0);

      return testApply(A, smoother, *X, *RHS, out, success);
    }

    // SmootherPrototype helper function
    void setupSmoother(RCP<Operator>& A, SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success) {
      Level level; MueLu::TestHelpers::Factory<SC,LO,GO,NO,LMO>::createSingleLevelHierarchy(level);
      level.Set("A", A);
      
      smoother.Setup(level);
    }

    // SmootherPrototype test
    ST::magnitudeType testApply_A125_X1_RHS0(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success) {
      Teuchos::RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(125);

      setupSmoother(A, smoother, out, success);
      return testApply_X1_RHS0(*A, smoother, out, success); // in MueLuTests::SmootherBase
    }

    // SmootherPrototype test
    ST::magnitudeType testApply_A125_X0_RandomRHS(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success) {
      Teuchos::RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(125);
      
      setupSmoother(A, smoother, out, success);
      return testApply_X0_RandomRHS(*A, smoother, out, success); // in MueLuTests::SmootherBase
    }

  } // namespace Smoother
  
} // MueLuTests


// TODO: add a test to check if Apply() throw an exception if vectors are not compatible with the smoother object (ie: with A)
// TODO: test Vector and MultiVector
// TODO: testApplyReduceError
