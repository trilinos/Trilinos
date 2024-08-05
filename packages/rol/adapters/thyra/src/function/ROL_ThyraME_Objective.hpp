// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(rol_x);
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
    const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(rol_x);
    ThyraVector<Real>  & thyra_dgdp = dynamic_cast<ThyraVector<Real>&>(rol_g);

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
      ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
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
