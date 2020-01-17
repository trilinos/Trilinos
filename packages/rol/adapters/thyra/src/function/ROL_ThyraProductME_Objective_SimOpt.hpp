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

#ifndef ROL_THYRAPRODUCTME_OBJECTIVE_SIMOPT
#define ROL_THYRAPRODUCTME_OBJECTIVE_SIMOPT

#include "Thyra_ProductVectorBase.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Types.hpp"

using namespace ROL;

template <class Real>
class ThyraProductME_Objective_SimOpt : public Objective_SimOpt<Real> {

public:


  ThyraProductME_Objective_SimOpt(Thyra::ModelEvaluatorDefaultBase<double>& thyra_model_, int g_index_, const std::vector<int>& p_indices_,Teuchos::RCP<Teuchos::ParameterList> params_ = Teuchos::null) :
    thyra_model(thyra_model_), g_index(g_index_), p_indices(p_indices_), params(params_) {
    updateValue = true;
    value_ = 0;
    x_ptr = Teuchos::null;
    if(params != Teuchos::null) {
      params->set<int>("Optimizer Iteration Number", -1);
    }
  };


  Real value(const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {

    if(updateValue) {
      const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
      const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);
      Teuchos::RCP< Thyra::VectorBase<Real> > g = Thyra::createMember<Real>(thyra_model.get_g_space(g_index));
      Teuchos::RCP<const Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();
      Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

      outArgs.set_g(g_index, g);
      for(std::size_t i=0; i<p_indices.size(); ++i)
        inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));
      inArgs.set_x(thyra_x.getVector());

      thyra_model.evalModel(inArgs, outArgs);

      value_ = ::Thyra::get_ele(*g,0);

      updateValue = false;
    }
    return value_;
  }

  void gradient_1(Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {

    const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
    const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);

    Teuchos::RCP<const  Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());
    ThyraVector<Real>  & thyra_dgdx = dynamic_cast<ThyraVector<Real>&>(g);

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();

    for(std::size_t i=0; i<p_indices.size(); ++i)
      inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));
    inArgs.set_x(thyra_x.getVector());

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

    Teuchos::RCP< Thyra::VectorBase<Real> > thyra_g;

    if(updateValue) {
      thyra_g = Thyra::createMember<Real>(thyra_model.get_g_space(g_index));
      outArgs.set_g(g_index, thyra_g);
    }

    for(std::size_t i=0; i<p_indices.size(); ++i) {
      const Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_support =
          outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, g_index);
      Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdx_orient;
      if (dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM))
        dgdx_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
      else if(dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM))
        dgdx_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
      else {
        ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
            "ROL::ThyraProductME_Objective: DgDx does support neither DERIV_MV_JACOBIAN_FORM nor DERIV_MV_GRADIENT_FORM forms");
      }

      outArgs.set_DgDx(g_index, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(thyra_dgdx.getVector(), dgdx_orient));
    }
    thyra_model.evalModel(inArgs, outArgs);

    if(updateValue) {
      value_ = ::Thyra::get_ele(*thyra_g,0);
      updateValue = false;
    }
  }

  void gradient_2(Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {

    const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
    const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);

    Teuchos::RCP<const  Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());
    ThyraVector<Real>  & thyra_dgdp = dynamic_cast<ThyraVector<Real>&>(g);

    //Teuchos::RCP<Thyra::MultiVectorBase<Real> > dgdp = thyra_dgdp.getVector();
    Teuchos::RCP< Thyra::ProductMultiVectorBase<Real> > prodvec_dgdp_p = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Real>>(thyra_dgdp.getVector());

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();

    for(std::size_t i=0; i<p_indices.size(); ++i)
      inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));
    inArgs.set_x(thyra_x.getVector());

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

    Teuchos::RCP< Thyra::VectorBase<Real> > thyra_g;

    if(updateValue) {
      thyra_g = Thyra::createMember<Real>(thyra_model.get_g_space(g_index));
      outArgs.set_g(g_index, thyra_g);
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

    if(updateValue) {
      value_ = ::Thyra::get_ele(*thyra_g,0);
      updateValue = false;
    }
  }

  void update( const Vector<Real> &/*u*/, const Vector<Real> &/*z*/, bool flag = true, int iter = -1) {
    updateValue = flag;
    if(params != Teuchos::null) {
      params->set<int>("Optimizer Iteration Number", iter);
      if(flag == true)
        params->set<bool>("Optimization Variables Changed",true);
    }
  }

public:
  bool updateValue;

private:
  Thyra::ModelEvaluatorDefaultBase<Real>& thyra_model;
  const int g_index;
  const std::vector<int> p_indices;
  Real value_;
  Teuchos::RCP<Vector<Real> > x_ptr;
  Teuchos::RCP<Teuchos::ParameterList> params;

};


#endif
