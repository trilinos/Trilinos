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
#include "Teuchos_VerbosityLevel.hpp"

using namespace ROL;

template <class Real>
class ThyraProductME_Objective_SimOpt : public Objective_SimOpt<Real> {

public:


  ThyraProductME_Objective_SimOpt(Thyra::ModelEvaluatorDefaultBase<double>& thyra_model_, int g_index_, const std::vector<int>& p_indices_,
      Teuchos::RCP<Teuchos::ParameterList> params_ = Teuchos::null, Teuchos::EVerbosityLevel verbLevel= Teuchos::VERB_HIGH) :
        thyra_model(thyra_model_), g_index(g_index_), p_indices(p_indices_), params(params_),
        out(Teuchos::VerboseObjectBase::getDefaultOStream()),
        verbosityLevel(verbLevel) {
    computeValue = computeGradient1 = computeGradient2 = true;
    value_ = 0;
    rol_u_ptr = rol_z_ptr = Teuchos::null;
    if(params != Teuchos::null) {
      params->set<int>("Optimizer Iteration Number", -1);
      params->set<Teuchos::RCP<Vector<Real> > >("Optimization Variable", Teuchos::null);
    }
  };


  Real value(const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {

#ifdef  HAVE_ROL_DEBUG
    //u and z should be updated in the update functions before calling applyAdjointJacobian_2
    TEUCHOS_ASSERT(!u_hasChanged(u));
    TEUCHOS_ASSERT(!z_hasChanged(z));
#endif

    if(verbosityLevel >= Teuchos::VERB_MEDIUM)
      *out << "ROL::ThyraProductME_Objective_SimOpt::value" << std::endl;

    if(!computeValue) {
      if(verbosityLevel >= Teuchos::VERB_HIGH)
        *out << "ROL::ThyraProductME_Objective_SimOpt::value, Skipping Computation of Value" << std::endl;
      return value_;
    }

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

    computeValue = false;

    return value_;
  }

  void gradient_1(Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {

#ifdef  HAVE_ROL_DEBUG
    //u and z should be updated in the update functions before calling gradient_1
    TEUCHOS_ASSERT(!u_hasChanged(u));
    TEUCHOS_ASSERT(!z_hasChanged(z));
#endif

    if(verbosityLevel >= Teuchos::VERB_MEDIUM)
      *out << "ROL::ThyraProductME_Objective_SimOpt::gradient_1" << std::endl;

    if(!computeGradient1) {
      if(verbosityLevel >= Teuchos::VERB_HIGH)
        *out << "ROL::ThyraProductME_Objective_SimOpt::gradient_1, Skipping Computation of Gradient 1" << std::endl;
      return g.set(*grad1_ptr_);
    }

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

    if(computeValue) {
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

    if(computeValue) {
      if(verbosityLevel >= Teuchos::VERB_HIGH)
        *out << "ROL::ThyraProductME_Objective_SimOpt::gradient_1, Computing Value" << std::endl;
      value_ = ::Thyra::get_ele(*thyra_g,0);
      computeValue = false;
    }

    if (Teuchos::is_null(grad1_ptr_))
      grad1_ptr_ = g.clone();
    grad1_ptr_->set(g);

    computeGradient1 = false;
  }

  void gradient_2(Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {

#ifdef  HAVE_ROL_DEBUG
    //u and z should be updated in the update functions before calling gradient_2
    TEUCHOS_ASSERT(!u_hasChanged(u));
    TEUCHOS_ASSERT(!z_hasChanged(z));
#endif

    if(verbosityLevel >= Teuchos::VERB_MEDIUM)
      *out << "ROL::ThyraProductME_Objective_SimOpt::gradient_2" << std::endl;


    if(!computeGradient2) {
      if(verbosityLevel >= Teuchos::VERB_HIGH)
        *out << "ROL::ThyraProductME_Objective_SimOpt::gradient_2, Skipping Computation of Gradient 2" << std::endl;
      return g.set(*grad2_ptr_);
    }

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

    if(computeValue) {
      if(verbosityLevel >= Teuchos::VERB_HIGH)
        *out << "ROL::ThyraProductME_Objective_SimOpt::gradient_2, Computing Value" << std::endl;
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

    if(computeValue) {
      value_ = ::Thyra::get_ele(*thyra_g,0);
      computeValue = false;
    }

    if (grad2_ptr_ == Teuchos::null)
      grad2_ptr_ = g.clone();
    grad2_ptr_->set(g);

    computeGradient2 = false;
  }


  void hessVec_11( Vector<Real> &hv, const Vector<Real> &v,
      const Vector<Real> &u,  const Vector<Real> &z, Real &/*tol*/ ) {

#ifdef  HAVE_ROL_DEBUG
    //u and z should be updated in the update functions before calling this function
    TEUCHOS_ASSERT(!u_hasChanged(u));
    TEUCHOS_ASSERT(!z_hasChanged(z));
#endif

    if(verbosityLevel >= Teuchos::VERB_MEDIUM)
      *out << "ROL::ThyraProductME_Objective_SimOpt::hessVec_11" << std::endl;

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();
    bool supports_deriv = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xx, g_index);

    if(supports_deriv) { //use derivatives computed by model evaluator

      const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
      const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);
      const ThyraVector<Real>  & thyra_v = dynamic_cast<const ThyraVector<Real>&>(v);

      Teuchos::RCP<const  Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());
      ThyraVector<Real>  & thyra_hv = dynamic_cast<ThyraVector<Real>&>(hv);

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();

      for(std::size_t i=0; i<p_indices.size(); ++i) {
        inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));
      }
      inArgs.set_x(thyra_x.getVector());
      inArgs.set_x_direction(thyra_v.getVector());

      Teuchos::RCP< Thyra::VectorBase<Real> > multiplier_g = Thyra::createMember<Real>(thyra_model.get_g_multiplier_space(g_index));
      Thyra::put_scalar(1.0, multiplier_g.ptr());
      inArgs.set_g_multiplier(g_index, multiplier_g);

      ROL_TEST_FOR_EXCEPTION( !supports_deriv, std::logic_error, "ROL::ThyraProductME_Objective: H_xx product vector is not supported");
      outArgs.set_hess_vec_prod_g_xx(g_index, thyra_hv.getVector());

      thyra_model.evalModel(inArgs, outArgs);

    } else { //compute derivatives with 2nd-order finite differences

      Real gtol = std::sqrt(ROL_EPSILON<Real>());
      // Compute step length
      Real h = std::cbrt(ROL_EPSILON<Real>());
      if (v.norm() > h) {
        h *= std::max(1.0,u.norm()/v.norm());
      }
      // Evaluate gradient of first component at (u+hv,z)
      ROL::Ptr<Vector<Real> > unew = u.clone();
      unew->set(u);
      unew->axpy(h,v);
      this->update(*unew,z);
      hv.zero();
      this->gradient_1(hv,*unew,z,gtol);
      // Evaluate gradient of first component at (u-hv,z)
      ROL::Ptr<Vector<Real> > g = hv.clone();
      unew->axpy(-2.0*h,v);
      this->update(*unew,z);
      this->gradient_1(*g,*unew,z,gtol);
      // Compute Newton quotient
      hv.axpy(-1.0,*g);
      hv.scale(0.5/h);
    }
  }

  void hessVec_12( Vector<Real> &hv, const Vector<Real> &v,
      const Vector<Real> &u, const Vector<Real> &z, Real &/*tol*/ ) {

#ifdef  HAVE_ROL_DEBUG
    //u and z should be updated in the update functions before calling this function
    TEUCHOS_ASSERT(!u_hasChanged(u));
    TEUCHOS_ASSERT(!z_hasChanged(z));
#endif

    if(verbosityLevel >= Teuchos::VERB_MEDIUM)
      *out << "ROL::ThyraProductME_Objective_SimOpt::hessVec_12" << std::endl;

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();
    bool supports_deriv = true;
    for(std::size_t j=0; j<p_indices.size(); ++j)
      supports_deriv =  supports_deriv && outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xp, g_index, p_indices[j]);

    if(supports_deriv) { //use derivatives computed by model evaluator
      const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
      const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);
      const ThyraVector<Real>  & thyra_v = dynamic_cast<const ThyraVector<Real>&>(v);

      Teuchos::RCP<const  Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());
      Teuchos::RCP<const  Thyra::ProductVectorBase<Real> > thyra_prodvec_v = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_v.getVector());
      ThyraVector<Real>  & thyra_hv = dynamic_cast<ThyraVector<Real>&>(hv);

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();

      for(std::size_t i=0; i<p_indices.size(); ++i) {
        inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));
        inArgs.set_p_direction(p_indices[i], thyra_prodvec_v->getVectorBlock(i));
      }
      inArgs.set_x(thyra_x.getVector());

      Teuchos::RCP< Thyra::VectorBase<Real> > multiplier_g = Thyra::createMember<Real>(thyra_model.get_g_multiplier_space(g_index));
      Thyra::put_scalar(1.0, multiplier_g.ptr());
      inArgs.set_g_multiplier(g_index, multiplier_g);

      std::vector<Teuchos::RCP< Thyra::MultiVectorBase<Real> > > hv_vec(p_indices.size());

      hv_vec[0] = thyra_hv.getVector();
      for(std::size_t j=1; j<p_indices.size(); ++j) {
        hv_vec[j] = thyra_hv.getVector()->clone_v();
      }

      for(std::size_t j=0; j<p_indices.size(); ++j) {
        bool supports_deriv =   outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xp, g_index, p_indices[j]);
        ROL_TEST_FOR_EXCEPTION( !supports_deriv, std::logic_error, "ROL::ThyraProductME_Objective_SimOpt: H_xp product vector is not supported");
        outArgs.set_hess_vec_prod_g_xp(g_index,p_indices[j], hv_vec[j]);
      }
      thyra_model.evalModel(inArgs, outArgs);

      for(std::size_t j=1; j<p_indices.size(); ++j)
        hv_vec[0]->update(1.0, *hv_vec[j]);

    } else { //compute derivatives with 2nd-order finite differences

      Real gtol = std::sqrt(ROL_EPSILON<Real>());
      // Compute step length
      Real h = std::cbrt(ROL_EPSILON<Real>());
      if (v.norm() > h) {
        h *= std::max(1.0,u.norm()/v.norm());
      }
      // Evaluate gradient of first component at (u,z+hv)
      ROL::Ptr<Vector<Real> > znew = z.clone();
      znew->set(z);
      znew->axpy(h,v);
      this->update(u,*znew);
      hv.zero();
      this->gradient_1(hv,u,*znew,gtol);
      // Evaluate gradient of first component at (u,z-hv)
      ROL::Ptr<Vector<Real> > g = hv.clone();
      znew->axpy(-2.0*h,v);
      this->update(u,*znew);
      this->gradient_1(*g,u,*znew,gtol);
      // Compute Newton quotient
      hv.axpy(-1.0,*g);
      hv.scale(0.5/h);
    }
  }

  void hessVec_21( Vector<Real> &hv, const Vector<Real> &v,
      const Vector<Real> &u, const Vector<Real> &z, Real &/*tol*/ ) {

#ifdef  HAVE_ROL_DEBUG
    //u and z should be updated in the update functions before calling this function
    TEUCHOS_ASSERT(!u_hasChanged(u));
    TEUCHOS_ASSERT(!z_hasChanged(z));
#endif

    if(verbosityLevel >= Teuchos::VERB_MEDIUM)
      *out << "ROL::ThyraProductME_Objective_SimOpt::hessVec_21" << std::endl;


    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

    bool supports_deriv = true;
    for(std::size_t i=0; i<p_indices.size(); ++i)
      supports_deriv = supports_deriv &&  outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_px, g_index, p_indices[i]);


    if(supports_deriv) { //use derivatives computed by model evaluator
      const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
      const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);
      const ThyraVector<Real>  & thyra_v = dynamic_cast<const ThyraVector<Real>&>(v);

      Teuchos::RCP<const  Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());
      ThyraVector<Real>  & thyra_hv = dynamic_cast<ThyraVector<Real>&>(hv);

      Teuchos::RCP< Thyra::ProductMultiVectorBase<Real> > prodvec_hv = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Real>>(thyra_hv.getVector());

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();

      for(std::size_t i=0; i<p_indices.size(); ++i) {
        inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));
      }
      inArgs.set_x(thyra_x.getVector());
      inArgs.set_x_direction(thyra_v.getVector());

      Teuchos::RCP< Thyra::VectorBase<Real> > multiplier_g = Thyra::createMember<Real>(thyra_model.get_g_multiplier_space(g_index));
      Thyra::put_scalar(1.0, multiplier_g.ptr());
      inArgs.set_g_multiplier(g_index, multiplier_g);

      Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();

      for(std::size_t i=0; i<p_indices.size(); ++i) {
        bool supports_deriv =   outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_px, g_index, p_indices[i]);
        ROL_TEST_FOR_EXCEPTION( !supports_deriv, std::logic_error, "ROL::ThyraProductME_Objective_SimOpt: H_px product vector is not supported");
        outArgs.set_hess_vec_prod_g_px(g_index,p_indices[i], prodvec_hv->getNonconstMultiVectorBlock(i));
      }
      thyra_model.evalModel(inArgs, outArgs);

    } else { //compute derivatives with 2nd-order finite differences

      Real gtol = std::sqrt(ROL_EPSILON<Real>());
      // Compute step length
      Real h = std::cbrt(ROL_EPSILON<Real>());;
      if (v.norm() > h) {
        h *= std::max(1.0,u.norm()/v.norm());
      }
      // Evaluate gradient of first component at (u+hv,z)
      ROL::Ptr<Vector<Real> > unew = u.clone();
      unew->set(u);
      unew->axpy(h,v);
      this->update(*unew,z);
      hv.zero();
      this->gradient_2(hv,*unew,z,gtol);
      // Evaluate gradient of first component at (u-hv,z)
      ROL::Ptr<Vector<Real> > g = hv.clone();
      unew->axpy(-2.0*h,v);
      this->update(*unew,z);
      this->gradient_2(*g,*unew,z,gtol);
      // Compute Newton quotient
      hv.axpy(-1.0,*g);
      hv.scale(0.5/h);
    }
  }

  void hessVec_22( Vector<Real> &hv, const Vector<Real> &v,
      const Vector<Real> &u,  const Vector<Real> &z, Real &/*tol*/ ) {

#ifdef  HAVE_ROL_DEBUG
    //u and z should be updated in the update functions before calling this function
    TEUCHOS_ASSERT(!u_hasChanged(u));
    TEUCHOS_ASSERT(!z_hasChanged(z));
#endif

    if(verbosityLevel >= Teuchos::VERB_MEDIUM)
      *out << "ROL::ThyraProductME_Objective_SimOpt::hessVec_22" << std::endl;

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model.createOutArgs();
    bool supports_deriv = true;
    for(std::size_t i=0; i<p_indices.size(); ++i)
      for(std::size_t j=0; j<p_indices.size(); ++j)
        supports_deriv =  supports_deriv && outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index, p_indices[i], p_indices[j]);

    if(supports_deriv) { //use derivatives computed by model evaluator
      const ThyraVector<Real>  & thyra_p = dynamic_cast<const ThyraVector<Real>&>(z);
      const ThyraVector<Real>  & thyra_x = dynamic_cast<const ThyraVector<Real>&>(u);
      const ThyraVector<Real>  & thyra_v = dynamic_cast<const ThyraVector<Real>&>(v);

      Teuchos::RCP<const  Thyra::ProductVectorBase<Real> > thyra_prodvec_p = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_p.getVector());
      Teuchos::RCP<const  Thyra::ProductVectorBase<Real> > thyra_prodvec_v = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_v.getVector());
      ThyraVector<Real>  & thyra_hv = dynamic_cast<ThyraVector<Real>&>(hv);

      Teuchos::RCP< Thyra::ProductMultiVectorBase<Real> > prodvec_hv = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Real>>(thyra_hv.getVector());

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model.createInArgs();

      for(std::size_t i=0; i<p_indices.size(); ++i) {
        inArgs.set_p(p_indices[i], thyra_prodvec_p->getVectorBlock(i));
        inArgs.set_p_direction(p_indices[i], thyra_prodvec_v->getVectorBlock(i));
      }
      inArgs.set_x(thyra_x.getVector());
      Teuchos::RCP< Thyra::VectorBase<Real> > multiplier_g = Thyra::createMember<Real>(thyra_model.get_g_multiplier_space(g_index));
      Thyra::put_scalar(1.0, multiplier_g.ptr());
      inArgs.set_g_multiplier(g_index, multiplier_g);

      std::vector<std::vector<Teuchos::RCP< Thyra::MultiVectorBase<Real> > > > hv_vec(p_indices.size());

      for(std::size_t i=0; i<p_indices.size(); ++i) {
        hv_vec[i].resize(p_indices.size());
        hv_vec[i][0] = prodvec_hv->getNonconstMultiVectorBlock(i);
        for(std::size_t j=1; j<p_indices.size(); ++j) {
          hv_vec[i][j] = hv_vec[i][0]->clone_mv();
        }
      }

      for(std::size_t i=0; i<p_indices.size(); ++i) {
        for(std::size_t j=0; j<p_indices.size(); ++j) {
          bool supports_deriv =   outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index, p_indices[i], p_indices[j]);
          ROL_TEST_FOR_EXCEPTION( !supports_deriv, std::logic_error, "ROL::ThyraProductME_Objective_SimOpt: H_pp product vector is not supported");

          outArgs.set_hess_vec_prod_g_pp(g_index,p_indices[i], p_indices[j], hv_vec[i][j]);
        }
      }
      thyra_model.evalModel(inArgs, outArgs);

      for(std::size_t i=0; i<p_indices.size(); ++i) {
        for(std::size_t j=1; j<p_indices.size(); ++j)
          hv_vec[i][0]->update(1.0, *hv_vec[i][j]);
      }
    } else { //compute derivatives with 2nd-order finite differences

      Real gtol = std::sqrt(ROL_EPSILON<Real>());
      // Compute step length
      Real h = std::cbrt(ROL_EPSILON<Real>());
      if (v.norm() > h) {
        h *= std::max(1.0,u.norm()/v.norm());
      }
      // Evaluate gradient of first component at (u,z+hv)
      ROL::Ptr<Vector<Real> > znew = z.clone();
      znew->set(z);
      znew->axpy(h,v);
      update(u,*znew);
      hv.zero();
      gradient_2(hv,u,*znew,gtol);
      // Evaluate gradient of first component at (u,z-hv)
      ROL::Ptr<Vector<Real> > g = hv.clone();
      znew->axpy(-2.0*h,v);
      update(u,*znew);
      gradient_2(*g,u,*znew,gtol);
      // Compute Newton quotient
      hv.axpy(-1.0,*g);
      hv.scale(0.5/h);
    }
  }

  void update( const Vector<Real> &u, const Vector<Real> &z, bool /*flag*/ = true, int iter = -1) {
    if(z_hasChanged(z) || u_hasChanged(u)) {
      if(verbosityLevel >= Teuchos::VERB_HIGH)
        *out << "ROL::ThyraProductME_Objective_SimOpt::update, Either The State Or The Parameters Changed" << std::endl;
      computeValue = computeGradient1 = computeGradient2 = true;

      if (Teuchos::is_null(rol_z_ptr))
        rol_z_ptr = z.clone();
      rol_z_ptr->set(z);

      if (Teuchos::is_null(rol_u_ptr))
        rol_u_ptr = u.clone();
      rol_u_ptr->set(u);
    }

    if(params != Teuchos::null) {
      auto& z_stored_ptr = params->get<Teuchos::RCP<Vector<Real> > >("Optimization Variable");
      if(Teuchos::is_null(z_stored_ptr) || z_hasChanged(*z_stored_ptr)) {
        if(verbosityLevel >= Teuchos::VERB_HIGH)
          *out << "ROL::ThyraProductME_Objective_SimOpt::update, Signaling That Parameter Changed" << std::endl;
        params->set<bool>("Optimization Variables Changed", true);
        if(Teuchos::is_null(z_stored_ptr))
          z_stored_ptr = z.clone();
        z_stored_ptr->set(z);
      }
      params->set<int>("Optimizer Iteration Number", iter);
    }
  }

  bool z_hasChanged(const Vector<Real> &rol_z) const {
    bool changed = true;
    if (Teuchos::nonnull(rol_z_ptr)) {
      auto diff = rol_z.clone();
      diff->set(*rol_z_ptr);
      diff->axpy( -1.0, rol_z );
      Real norm = diff->norm();
      changed = (norm == 0) ? false : true;
    }
    return changed;
  }

  bool u_hasChanged(const Vector<Real> &rol_u) const {
    bool changed = true;
    if (Teuchos::nonnull(rol_u_ptr)) {
      auto diff = rol_u.clone();
      diff->set(*rol_u_ptr);
      diff->axpy( -1.0, rol_u );
      Real norm = diff->norm();
      changed = (norm == 0) ? false : true;
    }
    return changed;
  }

public:
  bool computeValue, computeGradient1, computeGradient2;

private:
  Thyra::ModelEvaluatorDefaultBase<Real>& thyra_model;
  const int g_index;
  const std::vector<int> p_indices;
  Real value_;
  Teuchos::RCP<Vector<Real> > grad1_ptr_;
  Teuchos::RCP<Vector<Real> > grad2_ptr_;
  Teuchos::RCP<Vector<Real> > rol_z_ptr;
  Teuchos::RCP<Vector<Real> > rol_u_ptr;
  Teuchos::RCP<Teuchos::ParameterList> params;
  Teuchos::RCP<Teuchos::FancyOStream> out;
  Teuchos::EVerbosityLevel verbosityLevel;

};


#endif
