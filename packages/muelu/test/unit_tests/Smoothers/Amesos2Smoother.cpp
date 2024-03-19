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
#include <Amesos2_config.h>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Amesos2Smoother.hpp"
#include "MueLu_TestHelpersSmoothers.hpp"

namespace MueLuTests {

// this namespace already has:  #include "MueLu_UseShortNames.hpp"
using namespace TestHelpers::Smoothers;

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Amesos2Smoother, NotSetup, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
#if defined HAVE_AMESOS2_KLU2 || defined HAVE_AMESOS2_SUPERLU
    testApplyNoSetup(Amesos2Smoother(), out, success);
#endif
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Amesos2Smoother, Apply_Correctness, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    Teuchos::RCP<Amesos2Smoother> smoother;
#ifdef HAVE_AMESOS2_KLU2
    smoother = Teuchos::rcp(new Amesos2Smoother("Klu"));
    testDirectSolver(*smoother, out, success);
#endif

#ifdef HAVE_AMESOS2_SUPERLU
    smoother = Teuchos::rcp(new Amesos2Smoother("Superlu"));
    testDirectSolver(*smoother, out, success);
#endif
  }
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Amesos2Smoother, NotSetup, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Amesos2Smoother, Apply_Correctness, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
