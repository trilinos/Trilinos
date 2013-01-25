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

#include "Thyra_ModelEvaluatorHelpers.hpp"

#include "Teuchos_Assert.hpp"

namespace Piro {

template <typename Scalar>
void PerformSolveImpl(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities)
{
  // Empty function in/out arguments
  sensitivities.clear();
  responses.clear();

  // Parse parameters
  const bool computeSensitivities = solveParams.get("Compute Sensitivities", false);

  // Setup output arguments
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs = piroModel.createOutArgs();
  {
    const int num_g = outArgs.Ng(); // Number of vectors of responses (including state at index num_g - 1)
    responses.resize(num_g);
    sensitivities.resize(num_g);

    for (int j = 0; j < num_g; ++j) {
      {
        const Teuchos::RCP<Thyra::VectorBase<Scalar> > g = Thyra::createMember(*piroModel.get_g_space(j));
        outArgs.set_g(j, g);
        responses[j] = g;
      }

      const int num_p = outArgs.Np(); // Number of vectors of parameters
      sensitivities[j].resize(num_p);

      if (computeSensitivities) {
        for (int l = 0; l < num_p; ++l) {
          const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
            outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
          const Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdp_orient =
            Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
          if (dgdp_support.supports(dgdp_orient)) {
            const Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar> dgdp =
              Thyra::create_DgDp_mv(piroModel, j, l, dgdp_orient);
            outArgs.set_DgDp(j, l, dgdp);
            sensitivities[j][l] = dgdp.getMultiVector();
          }
        }
      }
    }
  }

  // Solve the problem using the nominal values for the parameters
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = piroModel.getNominalValues();
  piroModel.evalModel(inArgs, outArgs);
}

template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities)
{
  PerformSolveImpl(piroModel, solveParams, responses, sensitivities);
}

} // namespace Piro
