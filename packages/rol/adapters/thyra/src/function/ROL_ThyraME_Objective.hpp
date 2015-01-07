// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_THYRAME_OBJECTIVE_H
#define ROL_THYRAME_OBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include <iostream>

/** \class ROL::ThyraME_Objective
    \brief Implements the ROL::Objective interface for a Thyra Model Evaluator Objective.
*/

namespace ROL {

template <class Real>
class ThyraME_Objective : public Objective<Real> {
public:

  ThyraME_Objective(Thyra::ModelEvaluatorDefaultBase<double>& thyra_model_, int g_index_ = 0, int p_index_ = 0) : thyra_model(thyra_model_), g_index(g_index_), p_index(p_index_){};

  /** \brief Compute value.

      This function returns the objective function value.
      @param[in]          rol_x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  Real value( const Vector<Real> &rol_x, Real &tol ) {
    const ThyraVector<Real>  & thyra_p = Teuchos::dyn_cast<const ThyraVector<Real> >(rol_x);
    Teuchos::RCP< Thyra::VectorBase<Real> > g = Thyra::createMember<Real>(thyra_model.get_g_space(g_index));

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

    inArgs.set_p(p_index, thyra_p.getVector());
    outArgs.set_g(g_index, g);

    thyra_model.evalModel(inArgs, outArgs);

    return ::Thyra::get_ele(*g,0);
  };

  /** \brief Compute gradient.

      This function returns the objective function gradient.
      @param[out]         rol_g   is the gradient.
      @param[in]          rol_x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  void gradient( Vector<Real> &rol_g, const Vector<Real> &rol_x, Real &tol ) {
    const ThyraVector<Real>  & thyra_p = Teuchos::dyn_cast<const ThyraVector<Real> >(rol_x);
    ThyraVector<Real>  & thyra_dgdp = Teuchos::dyn_cast<ThyraVector<Real> >(rol_g);

    Teuchos::RCP<Thyra::MultiVectorBase<Real> > dgdp = thyra_dgdp.getVector();

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

    const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
          outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, g_index, p_index);

    Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdp_orient;
    if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM))
      dgdp_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
    else if(dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM))
      dgdp_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
       "ROL::ThyraME_Objective: DgDp does support neither DERIV_MV_JACOBIAN_FORM nor DERIV_MV_GRADIENT_FORM forms");
    }

    inArgs.set_p(p_index, thyra_p.getVector());
    outArgs.set_DgDp(g_index,p_index, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(dgdp, dgdp_orient));

    thyra_model.evalModel(inArgs, outArgs);

  };

private:
  Thyra::ModelEvaluatorDefaultBase<Real>& thyra_model;
  int g_index, p_index;

}; // class Objective

} // namespace ROL
#endif
