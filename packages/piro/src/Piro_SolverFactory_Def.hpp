// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Piro_ConfigDefs.hpp"

#ifdef Piro_ENABLE_NOX
#include "Piro_NOXSolver.hpp"
#endif /* Piro_ENABLE_NOX */

#ifdef Piro_ENABLE_Rythmos
#include "Piro_RythmosSolver.hpp"
#endif /* Piro_ENABLE_Rythmos */

#include "Teuchos_TestForException.hpp"

#include <string>
#include <stdexcept>

namespace Piro {

template <typename Scalar>
Teuchos::RCP<Thyra::ResponseOnlyModelEvaluatorBase<Scalar> > SolverFactory::createSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model)
{
  Teuchos::RCP<Thyra::ResponseOnlyModelEvaluatorBase<Scalar> > result;

  const std::string &solverType = piroParams->get("Piro Solver", "NOX");

#ifdef Piro_ENABLE_NOX
  if (solverType == "NOX") {
    result = Teuchos::rcp(new NOXSolver<Scalar>(piroParams, model));
  } else
#endif /* Piro_ENABLE_NOX */
#ifdef Piro_ENABLE_Rythmos
  if (solverType == "Rythmos") {
    result = rcp(new RythmosSolver<Scalar>(piroParams, model));
  } else
#endif /* Piro_ENABLE_Rythmos */
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::runtime_error,
        "Error: Unknown Piro Solver Type: " << solverType);
  }

  return result;
}

} // namespace Piro
