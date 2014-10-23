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

/** @ingroup func_group
    \class ROL::Objective
    \brief Provides the interface to evaluate objective functions.

    ROL's objective function interface is designed for Fr$eacute;chet differentiable 
    functionals \f$f:\mathcal{X}\to\mathbb{R}\f$, where \f$\mathcal{X}\f$ is a Banach
    space.  The basic operator interace, to be implemented by the user, requires:
    \li #value -- objective function evaluation.

    It is strongly recommended that the user additionally overload:
    \li #gradient -- the objective function gradient -- the default is a finite-difference approximation;
    \li #hessVec  -- the action of the Hessian -- the default is a finite-difference approximation.

    The user may also overload:
    \li #update     -- update the objective function at each new iteration;
    \li #dirDeriv   -- compute the directional derivative -- the default is a finite-difference approximation;
    \li #invHessVec -- the action of the inverse Hessian;
    \li #precond    -- the action of a preconditioner for the Hessian.

    ---
*/


namespace ROL {

template <class Real>
class ThyraME_Objective : public Objective<Real> {
public:

  ThyraME_Objective(Thyra::ModelEvaluatorDefaultBase<double>& _thyra_model) : thyra_model(_thyra_model), g_index(0), p_index(0){
    inArgs = thyra_model.createInArgs();
    outArgs = thyra_model.createOutArgs();

    const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
          outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, g_index, p_index);

    if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM))
      dgdp_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
    else if(dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM))
      dgdp_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
       "ROL::ThyraME_Objective: DgDp does support neither DERIV_MV_JACOBIAN_FORM nor DERIV_MV_GRADIENT_FORM forms");
    }
  };

  /** \brief Compute value.

      This function returns the objective function value.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  Real value( const Vector<Real> &rol_x, Real &tol ) {
    ROL::ThyraVector<Real>  & thyra_p = Teuchos::dyn_cast<ROL::ThyraVector<Real> >(const_cast <ROL::Vector<Real> &>(rol_x));
    g = Thyra::createMember<Real>(thyra_model.get_g_space(g_index));

    inArgs.set_p(p_index, thyra_p.getVector());
    outArgs.set_g(g_index, g);

    thyra_model.evalModel(inArgs, outArgs);

    Thyra::DetachedVectorView<Real> val_view(g);
    return  val_view[0];
  };

  /** \brief Compute gradient.

      This function returns the objective function gradient.
      @param[out]         g   is the gradient.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  void gradient( Vector<Real> &rol_g, const Vector<Real> &rol_x, Real &tol ) {
    ROL::ThyraVector<Real>  & thyra_p = Teuchos::dyn_cast<ROL::ThyraVector<Real> >(const_cast <ROL::Vector<Real> &>(rol_x));
    ROL::ThyraVector<Real>  & thyra_dgdp = Teuchos::dyn_cast<ROL::ThyraVector<Real> >(const_cast <ROL::Vector<Real> &>(rol_g));

    dgdp = thyra_dgdp.getNonConstVector();

    inArgs.set_p(p_index, thyra_p.getVector());
    outArgs.set_DgDp(g_index,p_index, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(dgdp, dgdp_orient));

    std::cout << "\n\nOrientation: " << dgdp_orient << "\n\n" <<std::endl;

    thyra_model.evalModel(inArgs, outArgs);
  };

private:
  Thyra::ModelEvaluatorBase::InArgs<Real> inArgs;
  Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs;
  Thyra::ModelEvaluatorDefaultBase<Real>& thyra_model;
  Teuchos::RCP<Thyra::MultiVectorBase<Real> > dgdp;
  Teuchos::RCP< Thyra::VectorBase<Real> > g;
  Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdp_orient;
  int g_index, p_index;

}; // class Objective

} // namespace ROL
#endif
