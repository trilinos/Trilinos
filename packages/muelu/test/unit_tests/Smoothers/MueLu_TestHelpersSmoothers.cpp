// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_TestHelpersSmoothers.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Utilities.hpp"

#include "MueLu_Exceptions.hpp"

#include "MueLu_UseDefaultTypes.hpp"

// See .hpp file for description

namespace MueLuTests {
  namespace TestHelpers {
    namespace Smoothers {

#include "MueLu_UseShortNames.hpp"

    // SmootherPrototype test
    void testApplyNoSetup(SmootherPrototype const & smoother, Teuchos::FancyOStream & out, bool & success) {
      GO numGlobalElements = 125;
      RCP<const Map> map = MapFactory::Build(Parameters::getLib(), numGlobalElements, 0, Parameters::getDefaultComm());

      RCP<MultiVector> X   = MultiVectorFactory::Build(map,1);
      RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

      TEST_THROW(smoother.Apply(*X,*RHS), MueLu::Exceptions::RuntimeError);
    }

    // SmootherBase test
    ST::magnitudeType testApply(const Matrix& A, const SmootherBase & smoother, MultiVector & X, const MultiVector & RHS, Teuchos::FancyOStream & out, bool & success) {
      Array<ST::magnitudeType> norms(1);

      RHS.norm2(norms);
      out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

      Teuchos::Array<ST::magnitudeType> initialNorms(1); X.norm2(initialNorms);
      out << "||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

      smoother.Apply(X, RHS); // TODO: bool const &InitialGuessIsZero=false

      Teuchos::Array<ST::magnitudeType> finalNorms(1); X.norm2(finalNorms);
      out << "||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(25) << norms[0] << std::endl;

      Teuchos::Array<ST::magnitudeType> residualNorms = Utils::ResidualNorm(A, X, RHS);
      out << "||Residual|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorms[0] << std::endl;

      return residualNorms[0];
    }

    // SmootherBase test
    ST::magnitudeType testApply_X1_RHS0(const Matrix& A, const SmootherBase & smoother, Teuchos::FancyOStream & out, bool & success) {
      RCP<MultiVector> X   = MultiVectorFactory::Build(A.getDomainMap(),1);
      RCP<MultiVector> RHS = MultiVectorFactory::Build(A.getRangeMap(),1);
      X->putScalar((SC) 1.0);
      RHS->putScalar((SC) 0.0);

      return testApply(A, smoother, *X, *RHS, out, success);
    }

    // SmootherBase test
    ST::magnitudeType testApply_X0_RandomRHS(const Matrix& A, const SmootherBase & smoother, Teuchos::FancyOStream & out, bool & success) {
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
    void setupSmoother(RCP<Matrix>& A, SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success) {
      Level level; TestHelpers::Factory<SC,LO,GO,NO,LMO>::createSingleLevelHierarchy(level);
      level.Set("A", A);
      smoother.Setup(level);
    }

    // SmootherPrototype test
    ST::magnitudeType testApply_A125_X1_RHS0(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success) {
      Teuchos::RCP<Matrix> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(125);

      setupSmoother(A, smoother, out, success);
      return testApply_X1_RHS0(*A, smoother, out, success); // in MueLuTests::SmootherBase
    }

    // SmootherPrototype test
    ST::magnitudeType testApply_A125_X0_RandomRHS(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success) {
      Teuchos::RCP<Matrix> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(125);

      setupSmoother(A, smoother, out, success);
      return testApply_X0_RandomRHS(*A, smoother, out, success); // in MueLuTests::SmootherBase
    }

    void testDirectSolver(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success) {
      ST::magnitudeType residualNorms = testApply_A125_X0_RandomRHS(smoother, out, success);
      TEST_EQUALITY(residualNorms < 1e-12, true);
    }

  } // namespace Smoothers
  } // namespace  TestHelpers
} // namespace MueLuTests


// TODO: add a test to check if Apply() throw an exception if vectors are not compatible with the smoother object (ie: with A)
// TODO: test Vector and MultiVector
// TODO: testApplyReduceError
