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
#ifndef MUELU_TEST_HELPERS_SMOOTHERS_H
#define MUELU_TEST_HELPERS_SMOOTHERS_H

#include <Teuchos_FancyOStream.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_UseDefaultTypes.hpp"

// Helper functions to test if derived classes conforms to the SmootherBase and SmootherPrototype interfaces

namespace MueLuTests {
  namespace TestHelpers {
    namespace Smoothers {

#include "MueLu_UseShortNames.hpp"

      //! Test if Apply() throw an exception when called before Setup()
      void testApplyNoSetup(SmootherPrototype const & smoother, Teuchos::FancyOStream & out, bool & success);

      //! Apply smoother with Poisson1D(125), RHS=0 and X initialized to 1
      //  This method calls Setup() internally
      ST::magnitudeType testApply_A125_X1_RHS0(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success);

      //!
      //! Apply smoother with Poisson1D(125) RHS choosed randomly and X initialized to 0
      ST::magnitudeType testApply_A125_X0_RandomRHS(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success);

      //! Test if a smoother reduces effectively the error on a simple problem
      // TODO

      //! Test direct solver
      void testDirectSolver(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success);

    } // namespace Smoothers
  } // namespace TestHelpers
} // MueLuTests

#endif // MUELU_TEST_HELPERS_SMOOTHERS_H
