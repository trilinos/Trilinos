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
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"
#include "ROL_StdVector.hpp"
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

    if(computeJacobian1_) {
      // Create Jacobian
      const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(u);
      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();
      Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

      Teuchos::RCP< const Thyra::LinearOpWithSolveFactoryBase<Real> > lows_factory = thyra_model_->get_W_factory();
      TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));

      if(Teuchos::is_null(jacobian1_))
        jacobian1_ = thyra_model_->create_W_op();

      inArgs.set_p(0, thyra_p.getVector());
      inArgs.set_x(thyra_x.getVector());

      outArgs.set_W_op(jacobian1_);
      thyra_model_->evalModel(inArgs, outArgs);

      computeJacobian1_ = false;
    } else {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyJacobian_1, Skipping Jacobian Computation" << std::endl;
    }

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

    const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
    const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(u);
    const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);
    ROL::ThyraVector<Real>  & thyra_ijv = dynamic_cast<ROL::ThyraVector<Real>&>(ijv);

    Teuchos::RCP<Thyra::MultiVectorBase<Real> > thyra_ijv_ptr = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Real>>(thyra_ijv.getVector());

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

    inArgs.set_p(0, thyra_p.getVector());
    inArgs.set_x(thyra_x.getVector());

    // Create Jacobian and preconditioner
    Teuchos::RCP< const Thyra::LinearOpWithSolveFactoryBase<Real> > lows_factory = thyra_model_->get_W_factory();
    TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));

    if(Teuchos::is_null(jacobian1_))
      jacobian1_ = thyra_model_->create_W_op();

    Teuchos::RCP< const ::Thyra::DefaultLinearOpSource<Real> > losb = Teuchos::rcp(new ::Thyra::DefaultLinearOpSource<Real>(jacobian1_));
    Teuchos::RCP< ::Thyra::PreconditionerBase<Real> > prec;

    Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<Real> > prec_factory =  lows_factory->getPreconditionerFactory();
    if (Teuchos::nonnull(prec_factory)) {
      prec = prec_factory->createPrec();
    } else if (outArgs.supports( Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
      prec = thyra_model_->create_W_prec();
    }
    const Teuchos::RCP<Thyra::LinearOpWithSolveBase<Real> > jacobian = lows_factory->createOp();

    if(computeJacobian1_)
    {
      outArgs.set_W_op(jacobian1_);
      thyra_model_->evalModel(inArgs, outArgs);
      outArgs.set_W_op(Teuchos::null);
      inArgs.set_x(Teuchos::null);

      computeJacobian1_ = false;
    } else {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyInverseJacobian_1, Skipping Jacobian Computation" << std::endl;
    }

    if (Teuchos::nonnull(prec_factory))
      prec_factory->initializePrec(losb, prec.get());
    else if ( Teuchos::nonnull(prec) && (outArgs.supports( Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) ) {
      inArgs.set_x(thyra_x.getVector());
      outArgs.set_W_prec(prec);
      thyra_model_->evalModel(inArgs, outArgs);
    }

    if(Teuchos::nonnull(prec))
      Thyra::initializePreconditionedOp<Real>(*lows_factory,
          jacobian1_,
          prec,
          jacobian.ptr());
    else
      Thyra::initializeOp<Real>(*lows_factory, jacobian1_, jacobian.ptr());

    const Thyra::SolveCriteria<Real> solve_criteria;

    Thyra::solve(
        *jacobian,
        Thyra::NOTRANS,
        *thyra_v.getVector(),
        thyra_ijv_ptr.ptr(),
        Teuchos::ptr(&solve_criteria));
  }



  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
      const ROL::Vector<Real> &z, Real &tol) {

    if(verbosityLevel_ >= Teuchos::VERB_MEDIUM)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyAdjointJacobian_1" << std::endl;


    Teuchos::RCP< Thyra::LinearOpBase<Real> > lop;
    if(computeJacobian1_){
      if(Teuchos::is_null(jacobian1_))
        jacobian1_ = thyra_model_->create_W_op();

      const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
      const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(u);

      //Teuchos::RCP< Thyra::VectorBase<Real> > thyra_f = Thyra::createMember<Real>(thyra_model_->get_f_space());

      Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = thyra_model_->createInArgs();
      Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = thyra_model_->createOutArgs();

      inArgs.set_p(0, thyra_p.getVector());
      inArgs.set_x(thyra_x.getVector());

      // Create implicitly transpose Jacobian and preconditioner
      Teuchos::RCP< const Thyra::LinearOpWithSolveFactoryBase<Real> > lows_factory = thyra_model_->get_W_factory();
      TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));

      outArgs.set_W_op(jacobian1_);
      thyra_model_->evalModel(inArgs, outArgs);

      computeJacobian1_ = false;
    }
    else {
      if(verbosityLevel_ >= Teuchos::VERB_HIGH)
        *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyAdjointJacobian_1, Skipping Jacobian Computation" << std::endl;
    }


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

    const ROL::ThyraVector<Real>  & thyra_p = dynamic_cast<const ROL::ThyraVector<Real>&>(z);
    const ROL::ThyraVector<Real>  & thyra_x = dynamic_cast<const ROL::ThyraVector<Real>&>(u);
    const ROL::ThyraVector<Real>  & thyra_v = dynamic_cast<const ROL::ThyraVector<Real>&>(v);
    ROL::ThyraVector<Real>  & thyra_iajv = dynamic_cast<ROL::ThyraVector<Real>&>(iajv);
    Teuchos::RCP<Thyra::MultiVectorBase<Real> > thyra_iajv_ptr = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Real>>(thyra_iajv.getVector());

    Thyra::ModelEvaluatorBase::InArgs<Real> inArgs = availableAdjointModel_ ? thyra_adjointModel_->createInArgs() : thyra_model_->createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<Real> outArgs = availableAdjointModel_ ? thyra_adjointModel_->createOutArgs() : thyra_model_->createOutArgs();

    inArgs.set_p(0, thyra_p.getVector());
    inArgs.set_x(thyra_x.getVector());

    // Create implicitly transpose Jacobian and preconditioner
    Teuchos::RCP< const Thyra::LinearOpWithSolveFactoryBase<Real> > lows_factory = availableAdjointModel_ ? thyra_adjointModel_->get_W_factory() : thyra_model_->get_W_factory();
    TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));
    Teuchos::RCP< Thyra::LinearOpBase<Real> > lop;

    if(availableAdjointModel_) {
      if(Teuchos::is_null(adjointJacobian1_) || computeAdjointJacobian1_)
        adjointJacobian1_ = thyra_adjointModel_->create_W_op();
      lop = adjointJacobian1_;
    } else {
      if(Teuchos::is_null(jacobian1_))
        jacobian1_ = thyra_model_->create_W_op();
      lop = jacobian1_;
    }


    Teuchos::RCP< const ::Thyra::DefaultLinearOpSource<Real> > losb = Teuchos::rcp(new ::Thyra::DefaultLinearOpSource<Real>(lop));
    Teuchos::RCP< ::Thyra::PreconditionerBase<Real> > prec;

    Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<Real> > prec_factory =  lows_factory->getPreconditionerFactory();
    if (Teuchos::nonnull(prec_factory)) {
      prec = prec_factory->createPrec();
    } else if (outArgs.supports( Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
      prec = availableAdjointModel_ ? thyra_adjointModel_->create_W_prec() :  thyra_model_->create_W_prec();
    }
    const Teuchos::RCP<Thyra::LinearOpWithSolveBase<Real> > jacobian = lows_factory->createOp();

    bool compute_lop = (computeJacobian1_ && !availableAdjointModel_) || (computeAdjointJacobian1_ && availableAdjointModel_);

    if(compute_lop)
    {
      outArgs.set_W_op(lop);

      if(availableAdjointModel_) {
        thyra_adjointModel_->evalModel(inArgs, outArgs);
        computeAdjointJacobian1_ = false;
      }
      else {
        thyra_model_->evalModel(inArgs, outArgs);
        computeJacobian1_ = false;
      }
      outArgs.set_W_op(Teuchos::null);

    } else if (verbosityLevel_ >= Teuchos::VERB_HIGH)
      *out_ << "Piro::ThyraProductME_Constraint_SimOpt::applyInverseAdjointJacobian_1, Skipping Jacobian Computation" << std::endl;

    if (Teuchos::nonnull(prec_factory))
      prec_factory->initializePrec(losb, prec.get());
    else if ( Teuchos::nonnull(prec) && (outArgs.supports( Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) ) {
      outArgs.set_W_prec(prec);
      if(availableAdjointModel_)
        thyra_adjointModel_->evalModel(inArgs, outArgs);
      else
        thyra_model_->evalModel(inArgs, outArgs);
    }

    if(Teuchos::nonnull(prec)) {
      if(availableAdjointModel_) {
        Thyra::initializePreconditionedOp<Real>(*lows_factory,
          lop,
          prec,
          jacobian.ptr());
      } else {
        Thyra::initializePreconditionedOp<Real>(*lows_factory,
            Thyra::transpose<Real>(lop),
            Thyra::unspecifiedPrec<Real>(::Thyra::transpose<Real>(prec->getUnspecifiedPrecOp())),
            jacobian.ptr());
      }
    }
    else {
      if(availableAdjointModel_)
        Thyra::initializeOp<Real>(*lows_factory, lop, jacobian.ptr());
      else
        Thyra::initializeOp<Real>(*lows_factory, Thyra::transpose<Real>(lop), jacobian.ptr());
    }
    const Thyra::SolveCriteria<Real> solve_criteria;

    lop->apply(Thyra::NOTRANS, *thyra_v.getVector(), thyra_iajv.getVector().ptr(), 1.0, 0.0);

    thyra_iajv.getVector()->assign(0);
    Thyra::solve(
        *jacobian,
        Thyra::NOTRANS,
        *thyra_v.getVector(),
        thyra_iajv_ptr.ptr(),
        Teuchos::ptr(&solve_criteria));
  };

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
