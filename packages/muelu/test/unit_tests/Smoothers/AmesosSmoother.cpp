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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Amesos_config.h>
#include "MueLu_TestHelpers.hpp"
#include "MueLu_TestHelpersSmoothers.hpp"

#include "MueLu_AmesosSmoother.hpp"

namespace MueLuTests {

typedef double Scalar;
typedef int LocalOrdinal;
typedef int GlobalOrdinal;
typedef Xpetra::EpetraNode Node;
#include "MueLu_UseShortNames.hpp"

using namespace TestHelpers::Smoothers;

TEUCHOS_UNIT_TEST(AmesosSmoother, NotSetup) {
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {
    testApplyNoSetup(AmesosSmoother(), out, success);
  }
}

TEUCHOS_UNIT_TEST(AmesosSmoother, Apply_Correctness) {
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {
    Teuchos::RCP<AmesosSmoother> smoother;
#ifdef HAVE_AMESOS_KLU
    smoother = rcp(new AmesosSmoother("Klu"));
    testDirectSolver(*smoother, out, success);
#endif

#ifdef HAVE_AMESOS_SUPERLU
    smoother = rcp(new AmesosSmoother("Superlu"));
    testDirectSolver(*smoother, out, success);
#endif
  }
}

}  // namespace MueLuTests

// TODO
//  Test with invalid std::string
//  Test with invalid parameterList? (== a characterization test for Amesos)
//  Test if paramList takes into account
//  Check if all the defaults that are used by MueLu are tested
