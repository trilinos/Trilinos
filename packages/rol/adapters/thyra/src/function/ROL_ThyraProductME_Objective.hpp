// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_THYRAPRODUCTME_OBJECTIVE_H
#define ROL_THYRAPRODUCTME_OBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
//#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include <vector>
#include "Teuchos_VerbosityLevel.hpp"

/** \class ROL::ThyraProductME_Objective
    \brief Implements the ROL::Objective interface for a Thyra Model Evaluator Objective.
*/

namespace ROL {

template <class Real>
class ThyraProductME_Objective : public Objective<Real> {
public:

  ThyraProductME_Objective(Thyra::ModelEvaluator<double>& thyra_model_, int g_index_, const std::vector<int>& p_indices_,
      Teuchos::RCP<Teuchos::ParameterList> params_ = Teuchos::null, Teuchos::EVerbosityLevel verbLevel= Teuchos::VERB_HIGH) :
    thyra_model(thyra_model_), g_index(g_index_), p_indices(p_indices_), params(params_),
    out(Teuchos::VerboseObjectBase::getDefaultOStream()),
    verbosityLevel(verbLevel) {
    computeValue = computeGradient = true;
    value_ = 0;
    if(Teuchos::nonnull(params)) {
      params->set<int>("Optimizer Iteration Number", -1);
      params->set<bool>("Compute State", true);
    }
  };

  /** \brief Compute value.

      This function returns the objective function value.
      @param[in]          rol_x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  Real value( const Vector<Real> &rol_x, Real &tol ) {

#ifdef  HAVE_ROL_DEBUG
    //x should be updated in the update function before calling value
    TEUCHOS_ASSERT(!x_hasChanged(rol_x));
#endif

    if(verbosityLevel >= Teuchos::VERB_MEDIUM)
      *out << "ROL::ThyraProductME_Objective::value" << std::endl;

    if(!computeValue) {
      if(verbosityLevel >= Teuchos::VERB_HIGH)
      *out << "ROL::ThyraProductME_Objective::value, Skipping Value Computation" << std::endl;
      return value_;
    }

    // Real norm = rol_x.norm();
    // std::cout << "Value norm: " << norm << std::endl;

    const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(rol_x);
    Teuchos::RCP< Thyra::VectorBase<Real> > g = Thyra::createMember<Real>(thyra_model.get_g_space(g_index));
    Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();


    outArgs.set_g(g_index, g);
    for(std::size_t i=0; i<p_indices.size(); ++i)
      inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));

    thyra_model.evalModel(inArgs, outArgs);

    if ((params == Teuchos::null) || !params->isParameter("State Solve Converged") || params->get<bool>("State Solve Converged"))
      value_ = ::Thyra::get_ele(*g,0);
    else
      value_ = 1.0e100;

    computeValue = false;

    return value_;
  };

  /** \brief Compute gradient.

      This function returns the objective function gradient.
      @param[out]         rol_g   is the gradient.
      @param[in]          rol_x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  void gradient( Vector<Real> &rol_g, const Vector<Real> &rol_x, Real &tol ) {

#ifdef  HAVE_ROL_DEBUG
    //x should be updated in the update function before calling value
    TEUCHOS_ASSERT(!x_hasChanged(rol_x));
#endif

    if(verbosityLevel >= Teuchos::VERB_MEDIUM)
      *out << "ROL::ThyraProductME_Objective::gradient" << std::endl;

    if(!computeGradient) {
      if(verbosityLevel >= Teuchos::VERB_HIGH)
        *out << "ROL::ThyraProductME_Objective::gradient, Skipping Gradient Computation" << std::endl;
      return rol_g.set(*grad_ptr_);
    }

    // Real norm = rol_x.norm();
    // std::cout << "In Gradient, Value norm: " << norm << std::endl;

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

    if(computeValue) {
      if(verbosityLevel >= Teuchos::VERB_HIGH)
        *out << "ROL::ThyraProductME_Objective::gradient, Computing Value" << std::endl;
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

    if(computeValue) {
      value_ = ::Thyra::get_ele(*g,0);
      computeValue = false;
    }
    
    if (grad_ptr_ == Teuchos::null)
      grad_ptr_ = rol_g.clone();
    grad_ptr_->set(rol_g);
    
    computeGradient = false;
  };

  void update( const Vector<Real> & x, bool flag = true, int iter = -1 ) {
    if(Teuchos::nonnull(params)) {
      params->set<int>("Optimizer Iteration Number", iter);
    }

    if(x_hasChanged(x)) {

      if(verbosityLevel >= Teuchos::VERB_HIGH)
        *out << "ROL::ThyraProductME_Objective::update, The Parameter Changed" << std::endl;
      computeValue = computeGradient = true;

      if(Teuchos::nonnull(params)) {
        params->set<bool>("Compute State", true);
        params->set<bool>("Optimization Variables Changed", true);
      }
    }
  }
  
  bool x_hasChanged(const Vector<Real> &rol_x) {
    bool changed = true;
    if (Teuchos::nonnull(rol_x_ptr)) {
      rol_x_ptr->axpy( -1.0, rol_x );
      Real norm = rol_x_ptr->norm();
      changed = (norm == 0) ? false : true;
    } else {
      rol_x_ptr = rol_x.clone();
    }
    rol_x_ptr->set(rol_x);
    return changed;
  }  

public:
  bool computeValue, computeGradient;

private:
  Thyra::ModelEvaluator<Real>& thyra_model;
  const int g_index;
  const std::vector<int> p_indices;
  Real value_;
  Teuchos::RCP<Vector<Real> > grad_ptr_;
  Teuchos::RCP<Vector<Real> > rol_x_ptr;
  Teuchos::RCP<Teuchos::ParameterList> params;
  Teuchos::RCP<Teuchos::FancyOStream> out;
  Teuchos::EVerbosityLevel verbosityLevel;

}; // class Objective

} // namespace ROL
#endif
