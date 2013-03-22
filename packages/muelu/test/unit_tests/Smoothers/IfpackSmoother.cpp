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
#include "MueLu_TestHelpers.hpp"
#include "MueLu_TestHelpersSmoothers.hpp"

#include "MueLu_IfpackSmoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

/*
  Comments about tests with hard coded results:
  1) Chebyshev smoothing must pass for any number of processors.
  2) Gauss-Seidel must pass for 1 and 4 processors.
  3) For any processor count except 1 and 4, the Gauss-Seidel test will
  report "passing", but this is only because the Teuchos test macro is skipped.
*/

namespace MueLuTests {

  using namespace TestHelpers::Smoothers;

  TEUCHOS_UNIT_TEST(IfpackSmoother, NotSetup)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {

      IfpackSmoother smoother("point relaxation stand-alone", Teuchos::ParameterList());
      testApplyNoSetup(smoother, out, success);

    }
  }

  // Tests interface to Ifpack's Gauss-Seidel preconditioner.
  TEUCHOS_UNIT_TEST(IfpackSmoother, HardCodedResult_GaussSeidel1)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {

      RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

      Teuchos::ParameterList paramList;
      paramList.set("relaxation: type", "Gauss-Seidel");
      paramList.set("relaxation: sweeps", (int) 1);
      paramList.set("relaxation: damping factor", (double) 1.0);
      paramList.set("relaxation: zero starting solution", false);

      IfpackSmoother smoother("point relaxation stand-alone", paramList);

      ST::magnitudeType residualNorms = testApply_A125_X1_RHS0(smoother, out, success);

      switch (comm->getSize()) {
      case 1:
      case 4:
        TEST_FLOATING_EQUALITY(residualNorms, 5.773502691896257e-01,1e-12);
        break;
      default:
        out << "Pass/Fail is checked only for 1 and 4 processes." << std::endl;
        break;
      } // switch

    }
  }

  // Tests interface to Ifpack's Gauss-Seidel preconditioner.
  TEUCHOS_UNIT_TEST(IfpackSmoother, HardCodedResult_GaussSeidel2)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {

      RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

      Teuchos::ParameterList paramList;
      paramList.set("relaxation: type", "Gauss-Seidel");
      paramList.set("relaxation: sweeps", (int) 10);
      paramList.set("relaxation: damping factor", (double) 1.0);
      paramList.set("relaxation: zero starting solution", false);

      IfpackSmoother smoother("point relaxation stand-alone", paramList);

      ST::magnitudeType residualNorms = testApply_A125_X1_RHS0(smoother, out, success);

      switch (comm->getSize()) {
      case 1:
        TEST_FLOATING_EQUALITY(residualNorms,8.326553652741774e-02,1e-12);
        break;
      case 4:
        TEST_FLOATING_EQUALITY(residualNorms,8.326553653078517e-02,1e-12);
        break;
      default:
        out << "Pass/Fail is checked only for 1 and 4 processes." << std::endl;
        break;

      } // switch

    }
  } // GaussSeidel

    // Tests interface to Ifpack's Chebyshev preconditioner.
  TEUCHOS_UNIT_TEST(IfpackSmoother, HardCodedResult_Chebyshev)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {

      // TODO This test should really calculate the reduction analytically.

      Teuchos::ParameterList paramList;
      paramList.set("chebyshev: degree", (int) 3);
      paramList.set("chebyshev: max eigenvalue", (double) 1.98476);
      paramList.set("chebyshev: min eigenvalue", (double) 1.0);
      paramList.set("chebyshev: ratio eigenvalue", (double) 20);
      paramList.set("chebyshev: zero starting solution", false);

      IfpackSmoother smoother("Chebyshev", paramList);

      ST::magnitudeType residualNorms = testApply_A125_X1_RHS0(smoother, out, success);

      TEST_FLOATING_EQUALITY(residualNorms, 5.269156e-01, 1e-7); // Compare to residual reported by ML

    }
  } // Chebyshev

  // Tests interface to Ifpack's ILU(0) preconditioner.
  TEUCHOS_UNIT_TEST(IfpackSmoother, HardCodedResult_ILU)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {

      // FIXME this will probably fail in parallel b/c it becomes block Jacobi

      IfpackSmoother smoother("ILU", Teuchos::ParameterList());

      ST::magnitudeType residualNorms = testApply_A125_X0_RandomRHS(smoother, out, success);

      RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
      if (comm->getSize() == 1) {
        TEST_EQUALITY(residualNorms < 1e-10, true);
      } else {
        out << "Pass/Fail is only checked in serial." << std::endl;
      }

    }
  } // ILU

} // namespace MueLuTests

//TODO: test if 10 its of ifpack == 10 apply call
//TODO:     TEST_EQUALITY(smoother->GetNIts(),10);
//TODO: check failed if TPETRA
