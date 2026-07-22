// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_THYRAPRODUCTME_OBJECTIVE_SIMOPT
#define PIRO_THYRAPRODUCTME_OBJECTIVE_SIMOPT

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Piro_ROL_ObserverBase.hpp"
#include "Piro_ProductModelEval.hpp"

namespace Piro {

template <class Real>
class ThyraProductME_Objective_SimOpt : public ROL::Objective_SimOpt<Real> {

public:


  ThyraProductME_Objective_SimOpt(const Teuchos::RCP<const Thyra::ModelEvaluator<Real>>& thyra_model, int g_index,
      Teuchos::ParameterList& piroParams, Teuchos::EVerbosityLevel verbLevel= Teuchos::VERB_HIGH,
      Teuchos::RCP<ROL_ObserverBase<Real>> observer = Teuchos::null) :
        thyra_model_(thyra_model), g_index_(g_index),
        optParams_(piroParams.sublist("Optimization Status")),
        out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
        verbosityLevel_(verbLevel), observer_(observer)  {
    write_interval_ = optParams_.get("Write Interval", 1);
    optParams_.set<int>("Optimizer Iteration Number", -1);
    useObjectiveRecoveryValue_ = optParams_.isParameter("Objective Recovery Value");
    objectiveRecoveryValue_ = useObjectiveRecoveryValue_ ? optParams_.get<Real>("Objective Recovery Value") : Real(0.0);
    updateType_ = ROL::UpdateType::Temp;
   };


  Real value(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Objective_SimOpt::value" << std::endl;

    if(objectiveStr_.isValueValid_) {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Objective_SimOpt::value, Skipping Computation of Value" << std::endl;
      return objectiveStr_.value_;
    }

    const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
    ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
    unew->set(u);
    const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);
    Teuchos::RCP< Thyra::VectorBase<Real> > g = Thyra::createMember<Real>(thyra_model_->get_g_space(g_index_));

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

    outArgs.set_g(g_index_, g);
    inArgs.set_p(0, thyra_p.getVector());

    inArgs.set_x(thyra_x.getVector());

    thyra_model_->evalModel(inArgs, outArgs);

    objectiveStr_.value_ = ::Thyra::get_ele(*g,0);

    //set value to (large) recovery value if solver did not converge
    if(useObjectiveRecoveryValue_ && optParams_.isParameter("State Solve Converged") && !optParams_.get<bool>("State Solve Converged")) {
      if(verbosityLevel_ >= Teuchos::VERB_LOW)
        *out_ << "Piro::ThyraProductME_Objective_SimOpt::value, Setting objective value to recovery value " << objectiveRecoveryValue_ << std::endl;
      objectiveStr_.value_ = objectiveRecoveryValue_;
    }

    objectiveStr_.isValueValid_ = true;

    if((updateType_ == ROL::UpdateType::Initial) && (observer_ != Teuchos::null)) {
      int iter = 0;
      const ROL::ThyraVector<Real>  & thyra_u = dynamic_cast<const ROL::ThyraVector<Real>&>(u);
      observer_->observeSolution(iter, *(thyra_u.getVector()), Teuchos::null, Teuchos::null, Teuchos::null);
      if(objectiveStr_.isValueValid_)
        observer_->observeResponse(iter);
    }

    return objectiveStr_.value_;
  }

  void gradient_1(ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Objective_SimOpt::gradient_1" << std::endl;

    TEUCHOS_ASSERT(objectiveStr_.areGradientsAllocated_);


    if(objectiveStr_.isGradient1Valid_ ) {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Objective_SimOpt::gradient_1, Skipping Computation of Gradient 1" << std::endl;
      return g.set(*objectiveStr_.gradient1_ptr_);
    }

    const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
    ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
    unew->set(u);
    const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);

    ROL::ThyraVector<Real>  & thyra_dgdx = dynamic_cast<ROL::ThyraVector<Real>&>(g);

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();

    inArgs.set_p(0, thyra_p.getVector());
    inArgs.set_x(thyra_x.getVector());

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

    Teuchos::RCP< Thyra::VectorBase<Real> > thyra_g;

    if(!objectiveStr_.isValueValid_) {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Objective_SimOpt::gradient_1, Computing Objective Value" << std::endl;
      thyra_g = Thyra::createMember<Real>(thyra_model_->get_g_space(g_index_));
      outArgs.set_g(g_index_, thyra_g);
    }

    const Thyra::ModelEvaluatorBase::DerivativeSupport dgdx_support =
        outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, g_index_);
    Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdx_orient;
    if (dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM))
      dgdx_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
    else if(dgdx_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM))
      dgdx_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
    else {
      ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Piro::ThyraProductME_Objective: DgDx does support neither DERIV_MV_JACOBIAN_FORM nor DERIV_MV_GRADIENT_FORM forms");
    }

    outArgs.set_DgDx(g_index_, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(thyra_dgdx.getVector(), dgdx_orient));

    thyra_model_->evalModel(inArgs, outArgs);

    if(!objectiveStr_.isValueValid_) {
      objectiveStr_.value_ = ::Thyra::get_ele(*thyra_g,0);
      
      //set value to (large) recovery value if solver did not converge
      if(useObjectiveRecoveryValue_ && optParams_.isParameter("State Solve Converged") && !optParams_.get<bool>("State Solve Converged")) {
        if(verbosityLevel_ >= Teuchos::VERB_LOW)
          *out_ << "Piro::ThyraProductME_Objective_SimOpt::gradient_1, Setting objective value to recovery value " << objectiveRecoveryValue_ << std::endl;
        objectiveStr_.value_ = objectiveRecoveryValue_;
      }
      
      objectiveStr_.isValueValid_ = true;
    }

    objectiveStr_.gradient1_ptr_->set(g);

    objectiveStr_.isGradient1Valid_ = true;
  }

  void gradient_2(ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Objective_SimOpt::gradient_2" << std::endl;

    TEUCHOS_ASSERT(objectiveStr_.areGradientsAllocated_);

    if(objectiveStr_.isGradient2Valid_) {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Objective_SimOpt::gradient_2, Skipping Computation of Gradient 2" << std::endl;
      return g.set(*objectiveStr_.gradient2_ptr_);
    }

    const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
    ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
    unew->set(u);
    const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);

    ROL::ThyraVector<Real>  & thyra_dgdp = dynamic_cast<ROL::ThyraVector<Real>&>(g);

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();

    inArgs.set_p(0, thyra_p.getVector());
    inArgs.set_x(thyra_x.getVector());

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

    Teuchos::RCP< Thyra::VectorBase<Real> > thyra_g;

    if(!objectiveStr_.isValueValid_) {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Objective_SimOpt::gradient_2, Computing Objective Value" << std::endl;
      thyra_g = Thyra::createMember<Real>(thyra_model_->get_g_space(g_index_));
      outArgs.set_g(g_index_, thyra_g);
    }

    Teuchos::RCP<Thyra::ProductMultiVectorBase<Real> > prodvec_dgdp =
        Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Real>>(thyra_dgdp.getVector());
    if ( !thyra_dgdp.getVector().is_null()) {
      if ( !prodvec_dgdp.is_null()) {
        Teuchos::RCP<const Piro::ProductModelEvaluator<Real>> model_PME = getProductModelEvaluator(thyra_model_);

        if ( !model_PME.is_null()) {
          Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Real>> dgdp_op =
              Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<Real>>(model_PME->create_DgDp_op(g_index_, 0, prodvec_dgdp));
          Thyra::ModelEvaluatorBase::Derivative<Real> dgdp_der(Teuchos::rcp_dynamic_cast<Thyra::LinearOpBase<Real>>(dgdp_op));
          outArgs.set_DgDp(g_index_, 0, dgdp_der);
        }
        else {
          ROL_TEST_FOR_EXCEPTION( true, std::logic_error, "Piro::ThyraProductME_Objective: dgdp is not supported for the used ModelEvaluator.");
        }
      }
      else {
        const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
            outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, g_index_, 0);
        Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdp_orient;
        if (dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM))
          dgdp_orient = Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM;
        else if(dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM))
          dgdp_orient = Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
        else {
          ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
              "Piro::ThyraProductME_Objective: DgDp does support neither DERIV_MV_JACOBIAN_FORM nor DERIV_MV_GRADIENT_FORM forms");
        }
        outArgs.set_DgDp(g_index_, 0, Thyra::ModelEvaluatorBase::DerivativeMultiVector<Real>(thyra_dgdp.getVector(), dgdp_orient));
      }
    }
    thyra_model_->evalModel(inArgs, outArgs);

    if(!objectiveStr_.isValueValid_) {
      objectiveStr_.value_ = ::Thyra::get_ele(*thyra_g,0);
      
      //set value to (large) recovery value if solver did not converge
      if(useObjectiveRecoveryValue_ && optParams_.isParameter("State Solve Converged") && !optParams_.get<bool>("State Solve Converged")) {
        if(verbosityLevel_ >= Teuchos::VERB_LOW)
          *out_ << "Piro::ThyraProductME_Objective_SimOpt::gradient_2, Setting objective value to recovery value " << objectiveRecoveryValue_ << std::endl;
        objectiveStr_.value_ = objectiveRecoveryValue_;
      }

      objectiveStr_.isValueValid_ = true;
    }

    objectiveStr_.gradient2_ptr_->set(g);
    objectiveStr_.isGradient2Valid_ = true;
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
      const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &/*tol*/ ) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Objective_SimOpt::hessVec_11" << std::endl;

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
    bool supports_deriv = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xx, g_index_);

    if(supports_deriv) { //use derivatives computed by model evaluator

      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
      ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
      unew->set(u);
      const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);
      const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);

      ROL::ThyraVector<Real>  & thyra_hv = dynamic_cast<ROL::ThyraVector<Real>&>(hv);

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();

      inArgs.set_p(0, thyra_p.getVector());
      inArgs.set_x(thyra_x.getVector());
      inArgs.set_x_direction(thyra_v.getVector());

      Teuchos::RCP< Thyra::VectorBase<Real> > multiplier_g = Thyra::createMember<Real>(thyra_model_->get_g_multiplier_space(g_index_));
      Thyra::put_scalar(1.0, multiplier_g.ptr());
      inArgs.set_g_multiplier(g_index_, multiplier_g);

      ROL_TEST_FOR_EXCEPTION( !supports_deriv, std::logic_error, "Piro::ThyraProductME_Objective: H_xx product vector is not supported");
      outArgs.set_hess_vec_prod_g_xx(g_index_, thyra_hv.getVector());

      thyra_model_->evalModel(inArgs, outArgs);

    } else { //compute derivatives with 2nd-order finite differences

      Real gtol = std::sqrt(ROL::ROL_EPSILON<Real>());
      // Compute step length
      Real h = std::cbrt(ROL::ROL_EPSILON<Real>());
      if (v.norm() > h) {
        h *= std::max(1.0,u.norm()/v.norm());
      }
      // Evaluate gradient of first component at (u+hv,z)
      ROL::Ptr<ROL::Vector<Real> > unew = u.clone();
      unew->set(u);
      unew->axpy(h,v);
      this->update(*unew,z,ROL::UpdateType::Temp);
      hv.zero();
      this->gradient_1(hv,*unew,z,gtol);
      // Evaluate gradient of first component at (u-hv,z)ROL::Ptr<ROL::Constraint_SimOpt<Scalar> > constr_ptr = ROL::makePtrFromRef(constr);
      ROL::Ptr<ROL::Vector<Real> > g = hv.clone();
      unew->axpy(-2.0*h,v);
      this->update(*unew,z,ROL::UpdateType::Temp);
      this->gradient_1(*g,*unew,z,gtol);
      // Compute Newton quotient
      hv.axpy(-1.0,*g);
      hv.scale(0.5/h);
      this->update(u,z,ROL::UpdateType::Temp);
    }
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
      const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &/*tol*/ ) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Objective_SimOpt::hessVec_12" << std::endl;

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
    bool supports_deriv = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xp, g_index_, 0);

    if(supports_deriv) { //use derivatives computed by model evaluator
      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
      ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
      unew->set(u);
      const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);
      const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);

      ROL::ThyraVector<Real>  & thyra_hv = dynamic_cast<ROL::ThyraVector<Real>&>(hv);

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();

      inArgs.set_p(0, thyra_p.getVector());
      inArgs.set_p_direction(0, thyra_v.getVector());
      inArgs.set_x(thyra_x.getVector());

      Teuchos::RCP< Thyra::VectorBase<Real> > multiplier_g = Thyra::createMember<Real>(thyra_model_->get_g_multiplier_space(g_index_));
      Thyra::put_scalar(1.0, multiplier_g.ptr());
      inArgs.set_g_multiplier(g_index_, multiplier_g);

      Teuchos::RCP< Thyra::MultiVectorBase<Real> >  hv_vec = thyra_hv.getVector();

      ROL_TEST_FOR_EXCEPTION( !outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_xp, g_index_, 0), 
        std::logic_error, "Piro::ThyraProductME_Objective_SimOpt: H_xp product vector is not supported");
      outArgs.set_hess_vec_prod_g_xp(g_index_,0, hv_vec);

      thyra_model_->evalModel(inArgs, outArgs);

    } else { //compute derivatives with 2nd-order finite differences

      Real gtol = std::sqrt(ROL::ROL_EPSILON<Real>());
      // Compute step length
      Real h = std::cbrt(ROL::ROL_EPSILON<Real>());
      if (v.norm() > h) {
        h *= std::max(1.0,u.norm()/v.norm());
      }
      // Evaluate gradient of first component at (u,z+hv)
      ROL::Ptr<ROL::Vector<Real> > znew = z.clone();
      znew->set(z);
      znew->axpy(h,v);
      this->update(u,*znew,ROL::UpdateType::Temp);
      hv.zero();
      this->gradient_1(hv,u,*znew,gtol);
      // Evaluate gradient of first component at (u,z-hv)
      ROL::Ptr<ROL::Vector<Real> > g = hv.clone();
      znew->axpy(-2.0*h,v);
      this->update(u,*znew,ROL::UpdateType::Temp);
      this->gradient_1(*g,u,*znew,gtol);
      // Compute Newton quotient
      hv.axpy(-1.0,*g);
      hv.scale(0.5/h);
      this->update(u,z,ROL::UpdateType::Temp);
    }
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
      const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &/*tol*/ ) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Objective_SimOpt::hessVec_21" << std::endl;


    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

    bool supports_deriv = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_px, g_index_, 0);


    if(supports_deriv) { //use derivatives computed by model evaluator
      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
      ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
      unew->set(u);
      const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);
      const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);

      ROL::ThyraVector<Real>  & thyra_hv = dynamic_cast<ROL::ThyraVector<Real>&>(hv);

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();

      inArgs.set_p(0, thyra_p.getVector());
      inArgs.set_x(thyra_x.getVector());
      inArgs.set_x_direction(thyra_v.getVector());

      Teuchos::RCP< Thyra::VectorBase<Real> > multiplier_g = Thyra::createMember<Real>(thyra_model_->get_g_multiplier_space(g_index_));
      Thyra::put_scalar(1.0, multiplier_g.ptr());
      inArgs.set_g_multiplier(g_index_, multiplier_g);

      ROL_TEST_FOR_EXCEPTION( !outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_px, g_index_, 0), 
        std::logic_error, "Piro::ThyraProductME_Objective_SimOpt: H_px product vector is not supported");
      outArgs.set_hess_vec_prod_g_px(g_index_,0, thyra_hv.getVector());

      thyra_model_->evalModel(inArgs, outArgs);

    } else { //compute derivatives with 2nd-order finite differences

      Real gtol = std::sqrt(ROL::ROL_EPSILON<Real>());
      // Compute step length
      Real h = std::cbrt(ROL::ROL_EPSILON<Real>());;
      if (v.norm() > h) {
        h *= std::max(1.0,u.norm()/v.norm());
      }
      // Evaluate gradient of first component at (u+hv,z)
      ROL::Ptr<ROL::Vector<Real> > unew = u.clone();
      unew->set(u);
      unew->axpy(h,v);
      this->update(*unew,z,ROL::UpdateType::Temp);
      hv.zero();
      this->gradient_2(hv,*unew,z,gtol);
      // Evaluate gradient of first component at (u-hv,z)
      ROL::Ptr<ROL::Vector<Real> > g = hv.clone();
      unew->axpy(-2.0*h,v);
      this->update(*unew,z,ROL::UpdateType::Temp);
      this->gradient_2(*g,*unew,z,gtol);
      // Compute Newton quotient
      hv.axpy(-1.0,*g);
      hv.scale(0.5/h);
      this->update(u,z,ROL::UpdateType::Temp);
    }
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v,
      const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &/*tol*/ ) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Objective_SimOpt::hessVec_22" << std::endl;

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
    bool supports_deriv = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index_, 0, 0);

    if(supports_deriv) { //use derivatives computed by model evaluator
      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
      ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
      unew->set(u);
      const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);
      const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);

      ROL::ThyraVector<Real>  & thyra_hv = dynamic_cast<ROL::ThyraVector<Real>&>(hv);

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();

      inArgs.set_p(0, thyra_p.getVector());
      inArgs.set_p_direction(0, thyra_v.getVector());
      inArgs.set_x(thyra_x.getVector());
      Teuchos::RCP< Thyra::VectorBase<Real> > multiplier_g = Thyra::createMember<Real>(thyra_model_->get_g_multiplier_space(g_index_));
      Thyra::put_scalar(1.0, multiplier_g.ptr());
      inArgs.set_g_multiplier(g_index_, multiplier_g);

      ROL_TEST_FOR_EXCEPTION( !outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_g_pp, g_index_, 0, 0), 
        std::logic_error, "Piro::ThyraProductME_Objective_SimOpt: H_pp product vector is not supported");

      outArgs.set_hess_vec_prod_g_pp(g_index_,0, 0, thyra_hv.getVector());

      thyra_model_->evalModel(inArgs, outArgs);

    } else { //compute derivatives with 2nd-order finite differences

      Real gtol = std::sqrt(ROL::ROL_EPSILON<Real>());
      // Compute step length
      Real h = std::cbrt(ROL::ROL_EPSILON<Real>());
      if (v.norm() > h) {
        h *= std::max(1.0,u.norm()/v.norm());
      }
      // Evaluate gradient of first component at (u,z+hv)
      ROL::Ptr<ROL::Vector<Real> > znew = z.clone();
      znew->set(z);
      znew->axpy(h,v);
      update(u,*znew,ROL::UpdateType::Temp);
      hv.zero();
      gradient_2(hv,u,*znew,gtol);
      // Evaluate gradient of first component at (u,z-hv)
      ROL::Ptr<ROL::Vector<Real> > g = hv.clone();
      znew->axpy(-2.0*h,v);
      update(u,*znew,ROL::UpdateType::Temp);
      gradient_2(*g,u,*znew,gtol);
      // Compute Newton quotient
      hv.axpy(-1.0,*g);
      hv.scale(0.5/h);
      this->update(u,z,ROL::UpdateType::Temp);
    }
  }

  void update( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, ROL::UpdateType type, int iter = -1) {
    
    if(verbosityLevel_ >= Teuchos::VERB_HIGH)
      *out_ << "Piro::ThyraProductME_Objective_SimOpt::update, UpdateType::" << ROL::UpdateTypeToString(type) << " Iter " << iter << std::endl;

    updateType_ = type;

    switch(type) {
      case ROL::UpdateType::Initial:    
        cached_objectiveStr_.allocateGradients(u, z);
        if(!tmp_objectiveStr_.areGradientsAllocated_)
          tmp_objectiveStr_.allocateGradients(u, z);
        objectiveStr_ = cached_objectiveStr_;
        break;
      case ROL::UpdateType::Accept: {
        tmp_objectiveStr_ = cached_objectiveStr_;
        cached_objectiveStr_ = objectiveStr_;

        TEUCHOS_ASSERT(!tmp_objectiveStr_.shareGradients(cached_objectiveStr_));
        if((observer_ != Teuchos::null) && (iter%write_interval_ == 0)) {
          const ROL::ThyraVector<Real>  & thyra_u = dynamic_cast<const ROL::ThyraVector<Real>&>(u);
          observer_->observeSolution(iter, *(thyra_u.getVector()), Teuchos::null, Teuchos::null, Teuchos::null);
          if(objectiveStr_.isValueValid_)
            observer_->observeResponse(iter);
        }
        break;
      }
      case ROL::UpdateType::Revert:
        objectiveStr_ = cached_objectiveStr_;
        break;
      case  ROL::UpdateType::Trial:
        objectiveStr_ = tmp_objectiveStr_;
        objectiveStr_.markAsNotValid();
        break;
      case ROL::UpdateType::Temp:
        if(!tmp_objectiveStr_.areGradientsAllocated_)
          tmp_objectiveStr_.allocateGradients(u, z);
        objectiveStr_ = tmp_objectiveStr_;
        objectiveStr_.markAsNotValid();
        break;
    }

    optParams_.set<int>("Optimizer Iteration Number", iter);
  }

  void update( const ROL::Vector<Real> &, const ROL::Vector<Real> &z, bool, int  = -1 ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl << "Piro::ThyraProductME_Objective_SimOpt::update:  Deprecated Update function, it should not be called." << std::endl);
  }

private:

  struct ObjectiveStruct {
    Real value_;
    Teuchos::RCP<ROL::Vector<Real> > gradient1_ptr_;
    Teuchos::RCP<ROL::Vector<Real> > gradient2_ptr_;
    bool isValueValid_;
    bool isGradient1Valid_;
    bool isGradient2Valid_;
    bool areGradientsAllocated_;

    ObjectiveStruct() : 
      value_(0), gradient1_ptr_(Teuchos::null), gradient2_ptr_(Teuchos::null),
      isValueValid_(false), isGradient1Valid_(false), isGradient2Valid_(false),
      areGradientsAllocated_(false) {}

    void allocateGradients(const ROL::Vector<Real>& gradient1, const ROL::Vector<Real>& gradient2) {
      gradient1_ptr_ = gradient1.clone();
      gradient2_ptr_ = gradient2.clone();
      areGradientsAllocated_ = true;
      markAsNotValid();
    }

    bool shareGradients(const ObjectiveStruct& objectiveStruct) {
      return (objectiveStruct.gradient1_ptr_.ptr() == gradient1_ptr_.ptr()) ||
        (objectiveStruct.gradient2_ptr_.ptr() == gradient2_ptr_.ptr());
    }

    void markAsNotValid() {
    isValueValid_ = isGradient1Valid_ = isGradient2Valid_ = false;
    }
  };



  const Teuchos::RCP<const Thyra::ModelEvaluator<Real>> thyra_model_;
  const int g_index_;
  Real objectiveRecoveryValue_;
  bool useObjectiveRecoveryValue_;
  ROL::UpdateType updateType_;

  ObjectiveStruct objectiveStr_, cached_objectiveStr_,  tmp_objectiveStr_;

  Teuchos::ParameterList& optParams_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::EVerbosityLevel verbosityLevel_;
  Teuchos::RCP<ROL_ObserverBase<Real>> observer_;
  int write_interval_;

};



}
#endif
