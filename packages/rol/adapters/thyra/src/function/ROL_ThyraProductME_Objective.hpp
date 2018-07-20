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

#ifndef ROL_THYRAPRODUCTME_OBJECTIVE_H
#define ROL_THYRAPRODUCTME_OBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
//#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include <vector>
#include <iostream>

/** \class ROL::ThyraProductME_Objective
    \brief Implements the ROL::Objective interface for a Thyra Model Evaluator Objective.
*/

namespace ROL {

template <class Real>
class ThyraProductME_Objective : public Objective<Real> {
public:

  ThyraProductME_Objective(Thyra::ModelEvaluatorDefaultBase<double>& thyra_model_, int g_index_, const std::vector<int>& p_indices_,Teuchos::RCP<Teuchos::ParameterList> params_ = Teuchos::null) :
    thyra_model(thyra_model_), g_index(g_index_), p_indices(p_indices_), params(params_) {
    valueUpdated = gradientUpdated = false;
    value_ = 0;
    x_ptr == Teuchos::null;
    grad_ptr = Teuchos::null;
    if(params != Teuchos::null)
      params->set<int>("Optimizer Iteration Number", -1);
  };

  /** \brief Compute value.

      This function returns the objective function value.
      @param[in]          rol_x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  Real value( const Vector<Real> &rol_x, Real &tol ) {

    if(!x_hasChanged(rol_x) &&  valueUpdated)
       return value_;

    const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(rol_x);
    Teuchos::RCP< Thyra::VectorBase<Real> > g = Thyra::createMember<Real>(thyra_model.get_g_space(g_index));
    Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();


    outArgs.set_g(g_index, g);
    for(std::size_t i=0; i<p_indices.size(); ++i)
      inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));

    thyra_model.evalModel(inArgs, outArgs);

    value_ = ::Thyra::get_ele(*g,0);

    valueUpdated = true;

    return value_;
  };

  /** \brief Compute gradient.

      This function returns the objective function gradient.
      @param[out]         rol_g   is the gradient.
      @param[in]          rol_x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  void gradient( Vector<Real> &rol_g, const Vector<Real> &rol_x, Real &tol ) {

    if( !x_hasChanged(rol_x) && gradientUpdated)
      return rol_g.set(*grad_ptr);

    if(params != Teuchos::null) {
      params->set<bool>("Update Functional", !valueUpdated);
      params->set<bool>("Update Functional Gradient", !gradientUpdated);
    }

    const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(rol_x);
    Teuchos::RCP<const  Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());
    ThyraVector<Real>  & thyra_dgdp = dynamic_cast<ThyraVector<Real>&>(rol_g);

    //Teuchos::RCP<Thyra::MultiVectorBase<Real> > dgdp = thyra_dgdp.getVector();
    Teuchos::RCP< Thyra::ProductMultiVectorBase<Real> > prodvec_dgdp_p = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Real>>(thyra_dgdp.getVector());

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();

    for(std::size_t i=0; i<p_indices.size(); ++i)
      inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

    Teuchos::RCP< Thyra::VectorBase<Real> > g;

    if(!valueUpdated) {
      g = Thyra::createMember<Real>(thyra_model.get_g_space(g_index));
      outArgs.set_g(g_index, g);
    }

    for(std::size_t i=0; i<p_indices.size(); ++i) {
      const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
            outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, g_index, p_indices[i]);
      Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdp_orient;
      if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM))
        dgdp_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
      else if(dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM))
        dgdp_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
      else {
        ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
         "ROL::ThyraProductME_Objective: DgDp does support neither DERIV_MV_JACOBIAN_FORM nor DERIV_MV_GRADIENT_FORM forms");
      }

      outArgs.set_DgDp(g_index,p_indices[i], Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(prodvec_dgdp_p->getNonconstMultiVectorBlock(i), dgdp_orient));
    }
    thyra_model.evalModel(inArgs, outArgs);

    if (grad_ptr == Teuchos::null)
      grad_ptr = rol_g.clone();
    grad_ptr->set(rol_g);

    if(!valueUpdated) {
      value_ = ::Thyra::get_ele(*g,0);
      valueUpdated = true;
    }
    gradientUpdated = true;
  };

  void update( const Vector<Real> & /*x*/, bool flag, int iter) {
     if(params != Teuchos::null)
       params->set<int>("Optimizer Iteration Number", iter);
  }

  bool x_hasChanged(const Vector<Real> &rol_x) {
    if (x_ptr == Teuchos::null) {
      x_ptr = rol_x.clone();
      x_ptr->set(rol_x);
      gradientUpdated = false;
      valueUpdated = false;
      return true;
    }
    else {
      x_ptr->axpy( -1.0, rol_x );
      Real norm = x_ptr->norm();
      x_ptr->set(rol_x);
      if (norm == 0) return false;
      else {
        gradientUpdated = false;
        valueUpdated = false;
        return true;
      }
    }
  }

public:
  bool gradientUpdated, valueUpdated;

private:
  Thyra::ModelEvaluatorDefaultBase<Real>& thyra_model;
  const int g_index;
  const std::vector<int> p_indices;
  Real value_;
  Teuchos::RCP<Vector<Real> > x_ptr, grad_ptr;
  Teuchos::RCP<Teuchos::ParameterList> params;

}; // class Objective

} // namespace ROL
#endif
