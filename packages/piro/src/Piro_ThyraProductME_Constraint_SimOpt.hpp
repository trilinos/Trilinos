// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_THYRAPRODUCTME_EQUALITYCONSTRAINT_SIMOPT
#define PIRO_THYRAPRODUCTME_EQUALITYCONSTRAINT_SIMOPT

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_ThyraVector.hpp"
#include "ROL_Constraint_SimOpt.hpp"

#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Piro_ROL_ObserverBase.hpp"

namespace Piro {

//! \brief ROL interface wrapper for Sacado SimOpt Constraint
/// @tparam Real 
template<class Real>
class ThyraProductME_Constraint_SimOpt : public ROL::Constraint_SimOpt<Real> {

public:

  ThyraProductME_Constraint_SimOpt(const Teuchos::RCP<const Thyra::ModelEvaluator<Real>>& thyra_model, 
      const Teuchos::RCP<const Thyra::ModelEvaluator<Real>>& thyra_adjointModel, 
      Teuchos::ParameterList& piroParams, Teuchos::EVerbosityLevel verbLevel= Teuchos::VERB_HIGH,
      Teuchos::RCP<ROL_ObserverBase<Real>> observer = Teuchos::null) :
        thyra_model_(thyra_model), thyra_adjointModel_(thyra_adjointModel),
        optParams_(piroParams.sublist("Optimization Status")),
        out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
        verbosityLevel_(verbLevel), observer_(observer) {
    thyra_solver_ = Teuchos::null;
    computeJacobian1_ = true;
    computeAdjointJacobian1_ = true;
    num_responses_ = -1;
    jacobian1_ = Teuchos::null;
    adjointJacobian1_ = Teuchos::null;
    isJacobianPrecValid_ = isJacobianSolverValid_ = false;
    isTransposedJacobianSolverValid_ = isAdjointJacobianSolverValid_ = false;

    availableAdjointModel_ = Teuchos::nonnull(thyra_adjointModel);
    optParams_.set<int>("Optimizer Iteration Number", -1);
    print_ = false;
  };

  void setExternalSolver(Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_solver) {
    thyra_solver_ = thyra_solver;
  }

  void setNumResponses(int num_responses) {
    num_responses_ = num_responses;
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::value" << std::endl;


    const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
    ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
    unew->set(u);
    const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);
    ROL::ThyraVector<Real>  & thyra_f = dynamic_cast<ROL::ThyraVector<Real>&>(c);

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

    outArgs.set_f(thyra_f.getVector());
    inArgs.set_p(0, thyra_p.getVector());
    inArgs.set_x(thyra_x.getVector());

    thyra_model_->evalModel(inArgs, outArgs);

    /*
    {
      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
      const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(u);
      const ROL::ThyraVector<Real>  & thyra_f = dynamic_cast<const ROL::ThyraVector<Real>&>(c);

      Thyra::ConstDetachedVectorView<Real> x_view(thyra_x.getVector());
      Thyra::ConstDetachedVectorView<Real> p_view(thyra_p.getVector());
      Thyra::ConstDetachedVectorView<Real> f_view(thyra_f.getVector());

      std::cout << "\nEnd of value... x:" << " ";
      for (std::size_t i=0; i<x_view.subDim(); ++i)
        std::cout << x_view(i) << " ";
      std::cout << "\np:" << " ";
      for (std::size_t i=0; i<p_view.subDim(); ++i)
        std::cout << p_view(i) << " ";
      std::cout << "\nf:" << " ";
      for (std::size_t i=0; i<f_view.subDim(); ++i)
        std::cout << f_view(i) << " ";

      std::cout << "Norm: " << c.norm() <<std::endl;
    }
    */
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
      const ROL::Vector<Real> &z, Real &tol) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyJacobian_1" << std::endl;

    computeJacobian1(u, z, "applyJacobian_1");

    const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);
    ROL::ThyraVector<Real>  & thyra_jv = dynamic_cast<ROL::ThyraVector<Real>&>(jv);
    jacobian1_->apply(Thyra::NOTRANS, *thyra_v.getVector(), thyra_jv.getVector().ptr(), 1.0, 0.0);

  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
      const ROL::Vector<Real> &z, Real &tol) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyJacobian_2" << std::endl;

    const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
    const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(u);
    const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);
    ROL::ThyraVector<Real>  & thyra_jv = dynamic_cast<ROL::ThyraVector<Real>&>(jv);

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

    inArgs.set_p(0, thyra_p.getVector());
    inArgs.set_x(thyra_x.getVector());

    Thyra::ModelEvaluatorBase::DerivativeSupport ds =  outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,0);

    if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
      auto dfdp_op = thyra_model_->create_DfDp_op(0);
      TEUCHOS_TEST_FOR_EXCEPTION(
          dfdp_op == Teuchos::null, std::logic_error,
          std::endl << "Piro::ThyraProductME_Constraint_SimOpt::applyJacobian_2():  " <<
          "Needed df/dp operator (" << 0 << ") is null!" << std::endl);
      outArgs.set_DfDp(0,dfdp_op);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(!ds.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP),
          std::logic_error,
          std::endl <<
          "Piro::ThyraProductME_Constraint_SimOpt::applyJacobian_2():  " <<
          "The code related to df/dp multivector has been commented out because never tested.  " <<
          std::endl);

    }

    thyra_model_->evalModel(inArgs, outArgs);
    thyra_jv.zero();

    Teuchos::RCP<const  Thyra::ProductVectorBase<Real> > thyra_prodvec_v = 
      Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Real>>(thyra_v.getVector());
    
    Thyra::ModelEvaluatorBase::Derivative<Real> dfdp_dv = outArgs.get_DfDp(0);

    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Real>> dfdp_block_op = 
      Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<Real>>(dfdp_dv.getLinearOp());    
    if (thyra_prodvec_v != Teuchos::null && dfdp_block_op != Teuchos::null) {
      for(int i=0; i<thyra_prodvec_v->productSpace()->numBlocks(); ++i) {
        auto dfdp_op = dfdp_block_op->getBlock(0, i);
        if (dfdp_op != Teuchos::null) {
          auto temp_jv_ptr = Teuchos::rcp_dynamic_cast<ROL::ThyraVector<Real>>(thyra_jv.clone());
          temp_jv_ptr->zero();
          dfdp_op->apply(Thyra::NOTRANS,*thyra_prodvec_v->getVectorBlock(i), temp_jv_ptr->getVector().ptr(),1.0, 0.0);
          thyra_jv.axpy(1.0, *temp_jv_ptr);
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(
              dfdp_op == Teuchos::null,
              std::logic_error,
              std::endl <<
              "Piro::ThyraProductME_Constraint_SimOpt::applyJacobian_2():  " <<
              "The code related to df/dp multivector has been commented out because never tested.  " <<
              std::endl);
        }
      }
    }
    else {
      auto dfdp_op = dfdp_dv.getLinearOp();
      auto dfdp = dfdp_dv.getMultiVector();

      if (dfdp_op != Teuchos::null) {
        auto temp_jv_ptr = Teuchos::rcp_dynamic_cast<ROL::ThyraVector<Real>>(thyra_jv.clone());
        temp_jv_ptr->zero();
        dfdp_op->apply(Thyra::NOTRANS,*thyra_v.getVector(), temp_jv_ptr->getVector().ptr(),1.0, 0.0);
        thyra_jv.axpy(1.0, *temp_jv_ptr);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(
            dfdp_op == Teuchos::null,
            std::logic_error,
            std::endl <<
            "Piro::ThyraProductME_Constraint_SimOpt::applyJacobian_2():  " <<
            "The code related to df/dp multivector has been commented out because never tested.  " <<
            std::endl);
      }
    }
  }


  void applyInverseJacobian_1(ROL::Vector<Real> &ijv,
      const ROL::Vector<Real> &v,
      const ROL::Vector<Real> &u,
      const ROL::Vector<Real> &z,
      Real &tol) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyInverseJacobian_1" << std::endl;

    initializeJacobianSolver(u, z, "applyInverseJacobian_1");
    solveLinearSystem(*jacobianSolver_, v, ijv, "applyInverseJacobian_1");
  }



  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
      const ROL::Vector<Real> &z, Real &tol) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyAdjointJacobian_1" << std::endl;

    computeJacobian1(u, z, "applyAdjointJacobian_1");

    const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);
    ROL::ThyraVector<Real>  & thyra_ajv = dynamic_cast<ROL::ThyraVector<Real>&>(ajv);
    jacobian1_->apply(Thyra::TRANS, *thyra_v.getVector(), thyra_ajv.getVector().ptr(), 1.0, 0.0);

  }

  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv,
      const ROL::Vector<Real> &v,
      const ROL::Vector<Real> &u,
      const ROL::Vector<Real> &z,
      Real &tol) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyInverseAdjointJacobian_1" << std::endl;

    // Solve with the Jacobian computed by the adjoint model, if available, or with the transpose of the Jacobian
    Teuchos::RCP<const Thyra::LinearOpWithSolveBase<Real>> solver;
    if(availableAdjointModel_) {
      initializeAdjointJacobianSolver(u, z, "applyInverseAdjointJacobian_1");
      solver = adjointJacobianSolver_;
    } else {
      initializeTransposedJacobianSolver(u, z, "applyInverseAdjointJacobian_1");
      solver = transposedJacobianSolver_;
    }

    ROL::ThyraVector<Real>  & thyra_iajv = dynamic_cast<ROL::ThyraVector<Real>&>(iajv);
    thyra_iajv.getVector()->assign(0);
    solveLinearSystem(*solver, v, iajv, "applyInverseAdjointJacobian_1");
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
      const ROL::Vector<Real> &z, Real &tol) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyAdjointJacobian_2" << std::endl;

    const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
    const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(u);
    const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);

    ROL::ThyraVector<Real>  & thyra_ajv = dynamic_cast<ROL::ThyraVector<Real>&>(ajv);
    thyra_ajv.zero();

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

    inArgs.set_x(thyra_x.getVector());
    inArgs.set_p(0, thyra_p.getVector());

    {
      // df/dp

      Thyra::ModelEvaluatorBase::DerivativeSupport ds =  outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, 0);
      // Determine which layout to use for df/dp.

      if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
        auto dfdp_op = thyra_model_->create_DfDp_op(0);
        TEUCHOS_TEST_FOR_EXCEPTION(
            dfdp_op == Teuchos::null, std::logic_error,
            std::endl << "Piro::ThyraProductME_Constraint_SimOpt::applyAdjointJacobian_2:  " <<
            "Needed df/dp operator (" << 0 << ") is null!" << std::endl);
        outArgs.set_DfDp(0,dfdp_op);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(
            !ds.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP),
            std::logic_error,
            std::endl <<
            "Piro::ThyraProductME_Constraint_SimOpt::applyAdjointJacobian_2():  " <<
            "The code related to df/dp multivector has been commented out because never tested.  " <<
            std::endl);
      }
    }

    thyra_model_->evalModel(inArgs, outArgs);

    Teuchos::RCP<Thyra::ProductVectorBase<Real> > thyra_prodvec_ajv = 
      Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<Real>>(thyra_ajv.getVector());
    
    Thyra::ModelEvaluatorBase::Derivative<Real> dfdp_dv = outArgs.get_DfDp(0);

    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Real>> dfdp_block_op = 
      Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<Real>>(dfdp_dv.getLinearOp());    
    if (thyra_prodvec_ajv != Teuchos::null && dfdp_block_op != Teuchos::null) {
      for(int i=0; i<thyra_prodvec_ajv->productSpace()->numBlocks(); ++i) {
        auto dfdp_op = dfdp_block_op->getBlock(0, i);
        if (dfdp_op != Teuchos::null) {
          dfdp_op->apply(Thyra::TRANS,*thyra_v.getVector(), thyra_prodvec_ajv->getNonconstVectorBlock(i).ptr(),1.0, 0.0);
          // Thyra::update(1.0,  *tmp, thyra_ajv.getMultiVector().ptr());
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(
              dfdp_op == Teuchos::null,
              std::logic_error,
              std::endl <<
              "Piro::ThyraProductME_Constraint_SimOpt::applyAdjointJacobian_2():  " <<
              "The code related to df/dp multivector has been commented out because never tested.  " <<
              std::endl);
        }
      }
    }
    else {
      auto dfdp_op = dfdp_dv.getLinearOp();      
      if (dfdp_op != Teuchos::null) {
        dfdp_op->apply(Thyra::TRANS,*thyra_v.getVector(), thyra_ajv.getMultiVector().ptr(),1.0, 0.0);
        // Thyra::update(1.0,  *tmp, thyra_ajv.getMultiVector().ptr());
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(
            dfdp_op == Teuchos::null,
            std::logic_error,
            std::endl <<
            "Piro::ThyraProductME_Constraint_SimOpt::applyAdjointJacobian_2():  " <<
            "The code related to df/dp multivector has been commented out because never tested.  " <<
            std::endl);
      }
    }
  }

  void solve(ROL::Vector<Real> &c,
      ROL::Vector<Real> &u,
      const ROL::Vector<Real> &z,
      Real &tol) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::solve" << std::endl;

    if(thyra_solver_.is_null()) {
      bool print = ROL::Constraint_SimOpt<Real>::print_;
      ROL::Constraint_SimOpt<Real>::print_ = false;  //unfortunately solve uses cout instead of user-provided stream, so temporarily disabling print
      ROL::Constraint_SimOpt<Real>::solve(c,u,z,tol);
      ROL::Constraint_SimOpt<Real>::print_ = print; 
    }
    else {
      if(this->zero_) { //SimOpt.Solve.Zero Initial Guess
        u.zero();
      }

      solve_update(u,z,ROL::UpdateType::Initial,0);

      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
      ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<ROL::ThyraVector<Real>&>(u);
      ROL::ThyraVector<Real>  & thyra_f = dynamic_cast<ROL::ThyraVector<Real>&>(c);

      //the last response will contain the solution
      Teuchos::RCP< Thyra::VectorBase<Real> > gx = Thyra::createMember<Real>(thyra_solver_->get_g_space(num_responses_));

      Teuchos::ParameterList origTestList;

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_solver_->createInArgs();
      Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_solver_->createOutArgs();

      outArgs.set_f(thyra_f.getVector());      
      outArgs.set_g(num_responses_, gx); //will contain the solution
      inArgs.set_p(0, thyra_p.getVector());

      inArgs.set_x(thyra_x.getVector());

      optParams_.set<bool>("Compute State", true);

      thyra_solver_->evalModel(inArgs, outArgs);

      //copy the solution into thyra_x  (u).
      Teuchos::RCP<const Thyra::VectorBase<Real> > gx_out = outArgs.get_g(num_responses_);
      if (Teuchos::nonnull(gx_out)) {
        Thyra::copy(*gx_out, thyra_x.getVector().ptr());
      }

      //Trial is always expected before accept. Accept does not reset computation flags.
      solve_update(u,z,ROL::UpdateType::Trial); 
      solve_update(u,z,ROL::UpdateType::Accept);
    }
  }


  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv,
      const ROL::Vector<Real> &w,
      const ROL::Vector<Real> &v,
      const ROL::Vector<Real> &u,
      const ROL::Vector<Real> &z,
      Real &tol) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyAdjointHessian_11" << std::endl;

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
    bool supports_deriv = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_xx);

    if(supports_deriv) { //use derivatives computed by model evaluator
      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
      ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
      unew->set(u);
      const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);
      const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);
      const ROL::ThyraVector<Real>  & thyra_w = dynamic_cast<const ROL::ThyraVector<Real>&>(w);

      ROL::ThyraVector<Real>  & thyra_ahwv = dynamic_cast<ROL::ThyraVector<Real>&>(ahwv);

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();

      inArgs.set_p(0, thyra_p.getVector());
      inArgs.set_x(thyra_x.getVector());
      inArgs.set_x_direction(thyra_v.getVector());

      inArgs.set_f_multiplier(thyra_w.getVector());

      ROL_TEST_FOR_EXCEPTION( !supports_deriv, std::logic_error, "Piro::ThyraProductME_Constraint: H_xx product vector is not supported");
      outArgs.set_hess_vec_prod_f_xx(thyra_ahwv.getVector());

      thyra_model_->evalModel(inArgs, outArgs);

    } else {  //compute derivatives with 2nd-order finite differences

      Real jtol = std::sqrt(ROL::ROL_EPSILON<Real>());
      // Compute step size
      Real h = std::cbrt(ROL::ROL_EPSILON<Real>());;
      if (v.norm() > h) {
        h *= std::max(1.0,u.norm()/v.norm());
      }
      // Evaluate Jacobian at (u+hv,z)
      ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
      unew->set(u);
      unew->axpy(h,v);
      this->update(*unew,z,ROL::UpdateType::Temp);
      applyAdjointJacobian_1(ahwv,w,*unew,z,jtol);
      // Evaluate Jacobian at (u-hv,z)
      ROL::Ptr<ROL::Vector<Real>> jv = ahwv.clone();
      unew->axpy(-2.*h,v);
      this->update(*unew,z,ROL::UpdateType::Temp);
      applyAdjointJacobian_1(*jv,w,*unew,z,jtol);
      // Compute Newton quotient
      ahwv.axpy(-1.0,*jv);
      ahwv.scale(0.5/h);
      this->update(u,z,ROL::UpdateType::Temp);
    }
  }


  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv,
      const ROL::Vector<Real> &w,
      const ROL::Vector<Real> &v,
      const ROL::Vector<Real> &u,
      const ROL::Vector<Real> &z,
      Real &/*tol*/) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyAdjointHessian_12" << std::endl;

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

    bool supports_deriv = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_px, 0);

    if(supports_deriv) {  //use derivatives computed by model evaluator

      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
      ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
      unew->set(u);
      const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);
      const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);
      const ROL::ThyraVector<Real>  & thyra_w = dynamic_cast<const ROL::ThyraVector<Real>&>(w);

      ROL::ThyraVector<Real>  & thyra_ahwv = dynamic_cast<ROL::ThyraVector<Real>&>(ahwv);

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();

      inArgs.set_p(0, thyra_p.getVector());
      inArgs.set_x(thyra_x.getVector());
      inArgs.set_x_direction(thyra_v.getVector());
      inArgs.set_f_multiplier(thyra_w.getVector());

      outArgs.set_hess_vec_prod_f_px(0, thyra_ahwv.getVector());

      thyra_model_->evalModel(inArgs, outArgs);

    } else {  //compute derivatives with 2nd-order finite differences


      Real jtol = std::sqrt(ROL::ROL_EPSILON<Real>());
      // Compute step size
      Real h = std::cbrt(ROL::ROL_EPSILON<Real>());
      if (v.norm() > h) {
        h *= std::max(1.0,u.norm()/v.norm());
      }
      // Evaluate Jacobian at (u+hv,z)
      ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
      unew->set(u);
      unew->axpy(h,v);
      this->update(*unew,z,ROL::UpdateType::Temp);
      applyAdjointJacobian_2(ahwv,w,*unew,z,jtol);
      // Evaluate Jacobian at (u - hv,z)
      ROL::Ptr<ROL::Vector<Real>> jv = ahwv.clone();
      unew->axpy(-2.0*h,v);
      this->update(*unew,z,ROL::UpdateType::Temp);
      applyAdjointJacobian_2(*jv,w,*unew,z,jtol);
      // Compute Newton quotient
      ahwv.axpy(-1.0,*jv);
      ahwv.scale(0.5/h);
      this->update(u,z,ROL::UpdateType::Temp);
    }
  }


  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv,
      const ROL::Vector<Real> &w,
      const ROL::Vector<Real> &v,
      const ROL::Vector<Real> &u,
      const ROL::Vector<Real> &z,
      Real &/*tol*/) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyAdjointHessian_21" << std::endl;

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
    bool supports_deriv = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_xp, 0);

    if(supports_deriv) { //use derivatives computed by model evaluator

      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
      ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
      unew->set(u);
      const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);
      const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);
      const ROL::ThyraVector<Real>  & thyra_w = dynamic_cast<const ROL::ThyraVector<Real>&>(w);

      ROL::ThyraVector<Real>  & thyra_ahwv = dynamic_cast<ROL::ThyraVector<Real>&>(ahwv);

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();

      inArgs.set_p(0, thyra_p.getVector());
      inArgs.set_p_direction(0, thyra_v.getVector());

      inArgs.set_x(thyra_x.getVector());
      inArgs.set_f_multiplier(thyra_w.getVector());

      outArgs.set_hess_vec_prod_f_xp(0, thyra_ahwv.getVector());

      thyra_model_->evalModel(inArgs, outArgs);

    } else {  //compute derivatives with 2nd-order finite differences

      Real jtol = std::sqrt(ROL::ROL_EPSILON<Real>());
      // Compute step size
      Real h = std::cbrt(ROL::ROL_EPSILON<Real>());
      if (v.norm() > h) {
        h *= std::max(1.0,u.norm()/v.norm());
      }
      // Evaluate Jacobian at (u,z+hv)
      ROL::Ptr<ROL::Vector<Real>> znew = z.clone();
      znew->set(z);
      znew->axpy(h,v);
      this->update(u,*znew,ROL::UpdateType::Temp);
      applyAdjointJacobian_1(ahwv,w,u,*znew,jtol);
      // Evaluate Jacobian at (u,z-hv)
      ROL::Ptr<ROL::Vector<Real>> jv = ahwv.clone();
      znew->axpy(-2.0*h,v);
      this->update(u,*znew,ROL::UpdateType::Temp);
      applyAdjointJacobian_1(*jv,w,u,*znew,jtol);
      // Compute Newton quotient
      ahwv.axpy(-1.0,*jv);
      ahwv.scale(0.5/h);
      this->update(u,z,ROL::UpdateType::Temp);
    }
  }

  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv,
      const ROL::Vector<Real> &w,
      const ROL::Vector<Real> &v,
      const ROL::Vector<Real> &u,
      const ROL::Vector<Real> &z,
      Real &/*tol*/) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyAdjointHessian_22" << std::endl;

    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();
    bool supports_deriv = outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_pp, 0, 0);

    if(supports_deriv) {  //use derivatives computed by model evaluator

      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
      ROL::Ptr<ROL::Vector<Real>> unew = u.clone();
      unew->set(u);
      const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(*unew);
      const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);
      const ROL::ThyraVector<Real>  & thyra_w = dynamic_cast<const ROL::ThyraVector<Real>&>(w);

      ROL::ThyraVector<Real>  & thyra_ahwv = dynamic_cast<ROL::ThyraVector<Real>&>(ahwv);

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();

      inArgs.set_p(0, thyra_p.getVector());
      inArgs.set_p_direction(0, thyra_v.getVector());
      inArgs.set_x(thyra_x.getVector());
      inArgs.set_f_multiplier(thyra_w.getVector());

      ROL_TEST_FOR_EXCEPTION( !outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_hess_vec_prod_f_pp, 0, 0),
        std::logic_error, "Piro::ThyraProductME_Constraint_SimOpt: H_pp product vector is not supported");

      outArgs.set_hess_vec_prod_f_pp(0, 0, thyra_ahwv.getVector());

      thyra_model_->evalModel(inArgs, outArgs);

    } else {  //compute derivatives with 2nd-order finite differences

      Real jtol = std::sqrt(ROL::ROL_EPSILON<Real>());
      // Compute step size
      Real h = std::cbrt(ROL::ROL_EPSILON<Real>());
      if (v.norm() > h) {
        h *= std::max(1.0,u.norm()/v.norm());
      }
      // Evaluate Jacobian at (u,z+hv)
      ROL::Ptr<ROL::Vector<Real>> znew = z.clone();
      znew->set(z);
      znew->axpy(h,v);
      this->update(u,*znew,ROL::UpdateType::Temp);
      applyAdjointJacobian_2(ahwv,w,u,*znew,jtol);
      // Evaluate Jacobian at (u,z-hv)
      ROL::Ptr<ROL::Vector<Real>> jv = ahwv.clone();
      znew->axpy(-2.0*h,v);
      this->update(u,*znew,ROL::UpdateType::Temp);
      applyAdjointJacobian_2(*jv,w,u,*znew,jtol);
      // Compute Newton quotient
      ahwv.axpy(-1.0,*jv);
      ahwv.scale(0.5/h);
      this->update(u,z,ROL::UpdateType::Temp);
    }
  }

  void solve_update(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, ROL::UpdateType type, int iter = -1) {
    if(verbosityLevel_ >= Teuchos::VERB_HIGH) {
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::solve_update, UpdateType::" << ROL::UpdateTypeToString(type) << " Iter " << iter << std::endl;
    }
    if(type != ROL::UpdateType::Accept)
      computeJacobian1_ = computeAdjointJacobian1_  = true;
  }

  void update_1( const ROL::Vector<Real> &u, ROL::UpdateType type, int iter = -1 ) {

    if(verbosityLevel_ >= Teuchos::VERB_HIGH) {
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::update_1, UpdateType::" << ROL::UpdateTypeToString(type) << " Iter " << iter << std::endl;
    }

    if(type != ROL::UpdateType::Accept)
      computeJacobian1_ = computeAdjointJacobian1_  = true;

    optParams_.set<int>("Optimizer Iteration Number", iter);
  }


  void update_2( const ROL::Vector<Real> &z, ROL::UpdateType type, int iter = -1 ) {

    if(verbosityLevel_ >= Teuchos::VERB_HIGH) {
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::update_2, UpdateType::" << ROL::UpdateTypeToString(type) << " Iter " << iter << std::endl;
    }
    
    if(type != ROL::UpdateType::Accept) {
      computeJacobian1_ = computeAdjointJacobian1_ = true;
      if(observer_ != Teuchos::null)
        observer_->parametersChanged();
    }
    
    optParams_.set<int>("Optimizer Iteration Number", iter);
  }


  void update_1( const ROL::Vector<Real> & , bool , int  = -1 ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl << "Piro::ThyraProductME_Constraint_SimOpt::update_1:  Deprecated Update function, it should not be called." << std::endl);
  }

  void update_2( const ROL::Vector<Real> & , bool , int  = -1 ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl << "Piro::ThyraProductME_Constraint_SimOpt::update_2:  Deprecated Update function, it should not be called." << std::endl);
  }

private:

  // The Jacobians and the linear solvers built from them are reused until the update functions
  // flag the Jacobians for recomputation (computeJacobian1_, computeAdjointJacobian1_).
  // Recomputing a Jacobian invalidates the preconditioners and the solvers built from it.

  //! Evaluates the Jacobian df/dx of the model at (u,z), unless it is up to date
  void computeJacobian1(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, const std::string& caller) {
    if(!computeJacobian1_) {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Constraint_SimOpt::" << caller << ", Skipping Jacobian Computation" << std::endl;
      return;
    }

    if(Teuchos::is_null(jacobian1_))
      jacobian1_ = thyra_model_->create_W_op();

    evalW(*thyra_model_, u, z, jacobian1_);

    computeJacobian1_ = false;
    isJacobianPrecValid_ = isJacobianSolverValid_ = isTransposedJacobianSolverValid_ = false;
  }

  //! Evaluates the Jacobian of the adjoint model, (df/dx)^T, at (u,z), unless it is up to date.
  //! The W operator of the adjoint model, and its preconditioner and solver, are created anew at every evaluation,
  //! because adjoint models can replace the matrix of W instead of filling it (e.g., they assign to it the transpose
  //! of the Jacobian, with a different graph), so that W could not be filled again.
  void computeAdjointJacobian1(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, const std::string& caller) {
    if(!computeAdjointJacobian1_) {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Constraint_SimOpt::" << caller << ", Skipping Adjoint Jacobian Computation" << std::endl;
      return;
    }

    adjointJacobian1_ = thyra_adjointModel_->create_W_op();
    adjointJacobianPrec_ = Teuchos::null;
    adjointJacobianSolver_ = Teuchos::null;

    evalW(*thyra_adjointModel_, u, z, adjointJacobian1_);

    computeAdjointJacobian1_ = false;
    isAdjointJacobianSolverValid_ = false;
  }

  //! Initializes the solver for the Jacobian df/dx at (u,z)
  void initializeJacobianSolver(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, const std::string& caller) {
    computeJacobian1(u, z, caller);
    if(isJacobianSolverValid_) {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Constraint_SimOpt::" << caller << ", Reusing Linear Solver" << std::endl;
      return;
    }

    const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Real>> lows_factory = thyra_model_->get_W_factory();
    TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));
    initializeJacobianPrec(u, z, *lows_factory);

    if(Teuchos::is_null(jacobianSolver_))
      jacobianSolver_ = lows_factory->createOp();
    if(Teuchos::nonnull(jacobianPrec_))
      Thyra::initializePreconditionedOp<Real>(*lows_factory, jacobian1_, jacobianPrec_, jacobianSolver_.ptr());
    else
      Thyra::initializeOp<Real>(*lows_factory, jacobian1_, jacobianSolver_.ptr());
    isJacobianSolverValid_ = true;
  }

  //! Initializes the solver for the transpose of the Jacobian df/dx at (u,z), used when the adjoint model is not available.
  //! The preconditioner is the transpose of the preconditioner of the Jacobian, so it must support the transposed apply.
  void initializeTransposedJacobianSolver(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, const std::string& caller) {
    computeJacobian1(u, z, caller);
    if(isTransposedJacobianSolverValid_) {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Constraint_SimOpt::" << caller << ", Reusing Linear Solver" << std::endl;
      return;
    }

    const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Real>> lows_factory = thyra_model_->get_W_factory();
    TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));
    initializeJacobianPrec(u, z, *lows_factory);

    if(Teuchos::is_null(transposedJacobianSolver_))
      transposedJacobianSolver_ = lows_factory->createOp();
    if(Teuchos::nonnull(jacobianPrec_))
      Thyra::initializePreconditionedOp<Real>(*lows_factory,
          Thyra::transpose<Real>(jacobian1_),
          Thyra::unspecifiedPrec<Real>(Thyra::transpose<Real>(jacobianPrec_->getUnspecifiedPrecOp())),
          transposedJacobianSolver_.ptr());
    else
      Thyra::initializeOp<Real>(*lows_factory, Thyra::transpose<Real>(jacobian1_), transposedJacobianSolver_.ptr());
    isTransposedJacobianSolverValid_ = true;
  }

  //! Initializes the solver for the Jacobian of the adjoint model at (u,z)
  void initializeAdjointJacobianSolver(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, const std::string& caller) {
    computeAdjointJacobian1(u, z, caller);
    if(isAdjointJacobianSolverValid_) {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Constraint_SimOpt::" << caller << ", Reusing Linear Solver" << std::endl;
      return;
    }

    const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Real>> lows_factory = thyra_adjointModel_->get_W_factory();
    TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));

    if(Teuchos::is_null(adjointJacobianPrec_))
      adjointJacobianPrec_ = createPrec(*thyra_adjointModel_, *lows_factory);
    initializePrec(*thyra_adjointModel_, *lows_factory, adjointJacobian1_, adjointJacobianPrec_, u, z);

    if(Teuchos::is_null(adjointJacobianSolver_))
      adjointJacobianSolver_ = lows_factory->createOp();
    if(Teuchos::nonnull(adjointJacobianPrec_))
      Thyra::initializePreconditionedOp<Real>(*lows_factory, adjointJacobian1_, adjointJacobianPrec_, adjointJacobianSolver_.ptr());
    else
      Thyra::initializeOp<Real>(*lows_factory, adjointJacobian1_, adjointJacobianSolver_.ptr());
    isAdjointJacobianSolverValid_ = true;
  }

  //! Initializes the preconditioner of the Jacobian df/dx, shared by the solvers of the Jacobian and of its transpose
  void initializeJacobianPrec(const ROL::Vector<Real> &u, const ROL::Vector<Real> &z,
      const Thyra::LinearOpWithSolveFactoryBase<Real>& lows_factory) {
    if(isJacobianPrecValid_) return;
    if(Teuchos::is_null(jacobianPrec_))
      jacobianPrec_ = createPrec(*thyra_model_, lows_factory);
    initializePrec(*thyra_model_, lows_factory, jacobian1_, jacobianPrec_, u, z);
    isJacobianPrecValid_ = true;
  }

  //! Creates the preconditioner with the preconditioner factory of the linear solver factory, if available,
  //! or with the model (W_prec). Returns null if neither is available.
  static Teuchos::RCP<Thyra::PreconditionerBase<Real>> createPrec(const Thyra::ModelEvaluator<Real>& model,
      const Thyra::LinearOpWithSolveFactoryBase<Real>& lows_factory) {
    const Teuchos::RCP<Thyra::PreconditionerFactoryBase<Real>> prec_factory = lows_factory.getPreconditionerFactory();
    if (Teuchos::nonnull(prec_factory))
      return prec_factory->createPrec();
    if (model.createOutArgs().supports(Thyra::ModelEvaluatorBase::OUT_ARG_W_prec))
      return model.create_W_prec();
    return Teuchos::null;
  }

  //! Computes the preconditioner prec of the operator op, the W operator of the model at (u,z)
  void initializePrec(const Thyra::ModelEvaluator<Real>& model, const Thyra::LinearOpWithSolveFactoryBase<Real>& lows_factory,
      const Teuchos::RCP<Thyra::LinearOpBase<Real>>& op, const Teuchos::RCP<Thyra::PreconditionerBase<Real>>& prec,
      const ROL::Vector<Real> &u, const ROL::Vector<Real> &z) const {
    if(Teuchos::is_null(prec)) return;
    const Teuchos::RCP<Thyra::PreconditionerFactoryBase<Real>> prec_factory = lows_factory.getPreconditionerFactory();
    if (Teuchos::nonnull(prec_factory)) {
      prec_factory->initializePrec(Teuchos::rcp(new Thyra::DefaultLinearOpSource<Real>(op)), prec.get());
    } else {
      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
      const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(u);
      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = model.createInArgs();
      Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = model.createOutArgs();
      inArgs.set_p(0, thyra_p.getVector());
      inArgs.set_x(thyra_x.getVector());
      outArgs.set_W_prec(prec);
      model.evalModel(inArgs, outArgs);
    }
  }

  //! Evaluates the W operator of the model at (u,z)
  static void evalW(const Thyra::ModelEvaluator<Real>& model, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z,
      const Teuchos::RCP<Thyra::LinearOpBase<Real>>& W) {
    const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
    const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(u);
    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = model.createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = model.createOutArgs();
    inArgs.set_p(0, thyra_p.getVector());
    inArgs.set_x(thyra_x.getVector());
    outArgs.set_W_op(W);
    model.evalModel(inArgs, outArgs);
  }

  //! Solves solver*sol = rhs, using the content of sol as initial guess. Prints a warning if the solver does not converge.
  //! The solve uses the default tolerance of the linear solver.
  void solveLinearSystem(const Thyra::LinearOpWithSolveBase<Real>& solver, const ROL::Vector<Real> &rhs,
      ROL::Vector<Real> &sol, const std::string& caller) const {
    const ROL::ThyraVector<Real>  & thyra_rhs = dynamic_cast<const ROL::ThyraVector<Real>&>(rhs);
    ROL::ThyraVector<Real>  & thyra_sol = dynamic_cast<ROL::ThyraVector<Real>&>(sol);
    const Thyra::SolveCriteria<Real> solve_criteria;
    const Thyra::SolveStatus<Real> solve_status = Thyra::solve<Real>(solver, Thyra::NOTRANS, *thyra_rhs.getVector(),
        thyra_sol.getVector().ptr(), Teuchos::ptr(&solve_criteria));
    if((solve_status.solveStatus == Thyra::SOLVE_STATUS_UNCONVERGED) && (verbosityLevel_ >= Teuchos::VERB_LOW))
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::" << caller << ", WARNING: the linear solver did not converge. "
            << solve_status.message << std::endl;
  }

  Teuchos::RCP<Thyra::PreconditionerBase<Real>> jacobianPrec_, adjointJacobianPrec_;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Real>> jacobianSolver_, transposedJacobianSolver_, adjointJacobianSolver_;
  bool isJacobianPrecValid_, isJacobianSolverValid_, isTransposedJacobianSolverValid_, isAdjointJacobianSolverValid_;

public:

  bool computeJacobian1_, computeAdjointJacobian1_;
  Teuchos::RCP<Thyra::ModelEvaluator<Real>> thyra_solver_;
  const Teuchos::RCP<const Thyra::ModelEvaluator<Real>> thyra_model_, thyra_adjointModel_;
  int num_responses_;
  Teuchos::ParameterList& optParams_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::RCP< Thyra::LinearOpBase<Real> > jacobian1_;
  Teuchos::RCP< Thyra::LinearOpBase<Real> > adjointJacobian1_;
  Teuchos::EVerbosityLevel verbosityLevel_;
  Teuchos::RCP<ROL_ObserverBase<Real>> observer_;
  bool availableAdjointModel_;
  bool print_;

};

}
#endif
