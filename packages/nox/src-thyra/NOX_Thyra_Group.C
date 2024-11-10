// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_CompilerCodeTweakMacros.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_PreconditionerBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "NOX_Common.H"
#include "NOX_Thyra_Group.H"    // class definition
#include "NOX_Abstract_MultiVector.H"
#include "NOX_Thyra_MultiVector.H"
#include "NOX_Assert.H"
#include "NOX_SolverStats.hpp"

NOX::Thyra::Group::
Group(const NOX::Thyra::Vector& initial_guess,
      const Teuchos::RCP< const ::Thyra::ModelEvaluator<double> >& model,
      const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weight_vector,
      const Teuchos::RCP<const ::Thyra::VectorBase<double> >& right_weight_vector,
      const Teuchos::RCP<::Thyra::VectorBase<double> >& inv_right_weight_vector,
      const bool rightScalingFirst):
  model_(model),
  rightScalingFirst_(rightScalingFirst),
  updatePreconditioner_(true),
  last_linear_solve_status_(NOX::Abstract::Group::NotConverged),
  last_linear_solve_num_iters_(0),
  last_linear_solve_achieved_tol_(0.0),
  use_pseudo_transient_terms_(false),
  x_dot_(Teuchos::null),
  alpha_(0.0),
  beta_(1.0),
  t_(-1.0),
  use_base_point_(false)
{
  x_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(initial_guess, DeepCopy));

  // To support implicit function scaling, all vectors must be copy
  // constructed/cloned from a NOX::Thyra::Vector that already has the
  // weight vector set or you must manually set the weighting vector.
  // Here we set the x_vec_ to have the weighting and clone that for
  // everything.
  if (nonnull(weight_vector)) {
    weight_vec_ = weight_vector;
    x_vec_->setWeightVector(weight_vec_);
  }

  if (nonnull(right_weight_vector)) {
    right_weight_vec_ = right_weight_vector;
    if (nonnull(inv_right_weight_vector)) {
      inv_right_weight_vec_ = inv_right_weight_vector;
    } else {
      inv_right_weight_vec_ = right_weight_vec_->clone_v();
      ::Thyra::reciprocal(*right_weight_vec_, inv_right_weight_vec_.ptr());
    }
    scaled_x_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));
  }
  else
    rightScalingFirst_ = false;

  f_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));
  newton_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));
  gradient_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));

  lop_ = model->create_W_op();

  // create jacobian operator
  lows_factory_ = model->get_W_factory();

  TEUCHOS_ASSERT(nonnull(lows_factory_));

  // Create jacobian with solver
  shared_jacobian_ = Teuchos::rcp(new NOX::SharedObject< ::Thyra::LinearOpWithSolveBase<double>, NOX::Thyra::Group >(lows_factory_->createOp()));

  losb_ = Teuchos::rcp(new ::Thyra::DefaultLinearOpSource<double>(lop_));

  // create preconditioner
  prec_factory_ = lows_factory_->getPreconditionerFactory();

  if (Teuchos::nonnull(prec_factory_)) {
    prec_ = prec_factory_->createPrec();
  } else if (model->createOutArgs().supports( ::Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
    prec_ = model->create_W_prec();
  }

  resetIsValidFlags();
}

NOX::Thyra::Group::
Group(const NOX::Thyra::Vector& initial_guess,
      const Teuchos::RCP< const ::Thyra::ModelEvaluator<double> >& model,
      const Teuchos::RCP< ::Thyra::LinearOpBase<double> >& linear_op,
      const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> >& lows_factory,
      const Teuchos::RCP< ::Thyra::PreconditionerBase<double> >& prec_op,
      const Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double> >& prec_factory,
      const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weight_vector,
      const Teuchos::RCP<const ::Thyra::VectorBase<double> >& right_weight_vector,
      const Teuchos::RCP<::Thyra::VectorBase<double> >& inv_right_weight_vector,
      const bool rightScalingFirst,
      const bool updatePreconditioner,
      const bool jacobianIsEvaluated):
  model_(model),
  lop_(linear_op),
  lows_factory_(lows_factory),
  prec_(prec_op),
  prec_factory_(prec_factory),
  rightScalingFirst_(rightScalingFirst),
  updatePreconditioner_(updatePreconditioner),
  last_linear_solve_status_(NOX::Abstract::Group::NotConverged),
  last_linear_solve_num_iters_(0),
  last_linear_solve_achieved_tol_(0.0),
  use_pseudo_transient_terms_(false),
  x_dot_(Teuchos::null),
  alpha_(0.0),
  beta_(1.0),
  t_(-1.0),
  use_base_point_(false)
{
  TEUCHOS_TEST_FOR_EXCEPTION(jacobianIsEvaluated && Teuchos::is_null(linear_op),std::runtime_error,
                             "ERROR - NOX::Thyra::Group(...) - linear_op is null but JacobianIsEvaluated is true. Impossible combination!");

  x_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(initial_guess, DeepCopy));

  // To support implicit function scaling, all vectors must be copy
  // constructed/cloned from a NOX::Thyra::Vector that already has the
  // weight vector set or you must manually set the weighting vector.
  // Here we set the x_vec_ to have the weighting and clone that for
  // everything.
  if (nonnull(weight_vector)) {
    weight_vec_ = weight_vector;
    x_vec_->setWeightVector(weight_vec_);
  }

  if (nonnull(right_weight_vector)) {
    right_weight_vec_ = right_weight_vector;
    if (nonnull(inv_right_weight_vector)) {
      inv_right_weight_vec_ = inv_right_weight_vector;
    } else {
      inv_right_weight_vec_ = right_weight_vec_->clone_v();
      ::Thyra::reciprocal(*right_weight_vec_, inv_right_weight_vec_.ptr());
    }
    scaled_x_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));
  }
  else
    rightScalingFirst_ = false;

  f_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));
  gradient_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));

  // Create jacobian with solver
  if (nonnull(lop_) && nonnull(lows_factory_)) {
    newton_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));
    shared_jacobian_ = Teuchos::rcp(new NOX::SharedObject< ::Thyra::LinearOpWithSolveBase<double>, NOX::Thyra::Group >(lows_factory_->createOp()));
    losb_ = Teuchos::rcp(new ::Thyra::DefaultLinearOpSource<double>(lop_));
  }

  if ( nonnull(prec_factory_) && is_null(prec_) )
    prec_ = prec_factory_->createPrec();

  resetIsValidFlags();

  if (jacobianIsEvaluated) {
    is_valid_jacobian_ = true;
    shared_jacobian_->getObject(this);
  }
}

NOX::Thyra::Group::Group(const NOX::Thyra::Group& source, NOX::CopyType type) :
  model_(source.model_),
  shared_jacobian_(source.shared_jacobian_),
  lop_(source.lop_),
  lows_factory_(source.lows_factory_),
  losb_(source.losb_),
  prec_(source.prec_),
  prec_factory_(source.prec_factory_),
  right_weight_vec_(source.right_weight_vec_),
  inv_right_weight_vec_(source.inv_right_weight_vec_),
  rightScalingFirst_(source.rightScalingFirst_),
  updatePreconditioner_(source.updatePreconditioner_),
  last_linear_solve_status_(source.last_linear_solve_status_),
  last_linear_solve_num_iters_(source.last_linear_solve_num_iters_),
  last_linear_solve_achieved_tol_(source.last_linear_solve_achieved_tol_),
  use_pseudo_transient_terms_(source.use_pseudo_transient_terms_),
  x_dot_(source.x_dot_),
  alpha_(source.alpha_),
  beta_(source.beta_),
  t_(source.t_),
  use_base_point_(false),
  base_point_(source.base_point_)
{

  x_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*source.x_vec_, type));
  f_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*source.f_vec_, type));
  if (nonnull(source.newton_vec_))
  newton_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*source.newton_vec_, type));
  gradient_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*source.gradient_vec_, type));

  if (nonnull(source.weight_vec_))
    weight_vec_ = source.weight_vec_;

  if (nonnull(source.right_weight_vec_))
    scaled_x_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*source.scaled_x_vec_, type));

  if (type == NOX::DeepCopy) {
    is_valid_f_ = source.is_valid_f_;
    is_valid_jacobian_ = source.is_valid_jacobian_;
    is_valid_newton_dir_ = source.is_valid_newton_dir_;
    is_valid_gradient_dir_ = source.is_valid_gradient_dir_;
    is_valid_lows_ = source.is_valid_lows_;

    // New copy takes ownership of the shared Jacobian for DeepCopy
    if (nonnull(shared_jacobian_))
      if (this->isJacobian())
    shared_jacobian_->getObject(this);
  }
  else if (type == NOX::ShapeCopy) {
    resetIsValidFlags();
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
               "NOX Error - Copy type is invalid!");
  }
}

NOX::Thyra::Group::~Group()
{ }

void NOX::Thyra::Group::resetIsValidFlags()
{
  is_valid_f_ = false;
  is_valid_jacobian_ = false;
  is_valid_newton_dir_ = false;
  is_valid_gradient_dir_ = false;
  is_valid_lows_ = false;
}

Teuchos::RCP<NOX::Abstract::Group> NOX::Thyra::Group::
clone(NOX::CopyType type) const
{
  Teuchos::RCP<NOX::Abstract::Group> newgrp =
    Teuchos::rcp(new NOX::Thyra::Group(*this, type));
  return newgrp;
}

NOX::Abstract::Group& NOX::Thyra::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const NOX::Thyra::Group&> (source));
}

NOX::Abstract::Group& NOX::Thyra::Group::operator=(const Group& source)
{

  // Copy the xVector
  *x_vec_ = *source.x_vec_;

  is_valid_f_ = source.is_valid_f_;
  is_valid_jacobian_ = source.is_valid_jacobian_;
  is_valid_newton_dir_ = source.is_valid_newton_dir_;
  is_valid_gradient_dir_ = source.is_valid_gradient_dir_;
  is_valid_lows_ = source.is_valid_lows_;

  if (this->isF())
    *f_vec_ = *(source.f_vec_);

  // Jacobian is shared, always assign shared object
  shared_jacobian_ = source.shared_jacobian_;
  lows_factory_ = source.lows_factory_;
  lop_ = source.lop_;
  losb_ = source.losb_;
  prec_factory_ = source.prec_factory_;
  prec_ = source.prec_;

  if (this->isNewton())
    *newton_vec_ = *(source.newton_vec_);

  if (this->isGradient())
    *gradient_vec_ = *(source.gradient_vec_);

  if (nonnull(source.weight_vec_))
    weight_vec_ = source.weight_vec_;

  if (nonnull(source.right_weight_vec_)) {
    right_weight_vec_ = source.right_weight_vec_;
    inv_right_weight_vec_ = source.inv_right_weight_vec_;
    *scaled_x_vec_ = *source.scaled_x_vec_;    
  }

  // If valid, this takes ownership of the shared Jacobian
  if (nonnull(shared_jacobian_))
    if (this->isJacobian())
      shared_jacobian_->getObject(this);

  rightScalingFirst_ = source.rightScalingFirst_;
  updatePreconditioner_ = source.updatePreconditioner_;

  last_linear_solve_status_ = source.last_linear_solve_status_;
  last_linear_solve_num_iters_ = source.last_linear_solve_num_iters_;
  last_linear_solve_achieved_tol_ = source.last_linear_solve_achieved_tol_;

  use_pseudo_transient_terms_ = source.use_pseudo_transient_terms_;
  x_dot_ = source.x_dot_;
  alpha_ = source.alpha_;
  beta_ = source.beta_;
  t_ = source.t_;

  use_base_point_ = source.use_base_point_;
  base_point_ = source.base_point_;

  return *this;
}


Teuchos::RCP<const ::Thyra::VectorBase<double> >
NOX::Thyra::Group::get_current_x() const
{
  if (is_null(x_vec_))
    return Teuchos::null;
  return x_vec_->getThyraRCPVector();
}


Teuchos::RCP< ::Thyra::LinearOpBase<double> >
NOX::Thyra::Group::getNonconstJacobianOperator()
{
  NOX_ASSERT(nonnull(shared_jacobian_));
  NOX_ASSERT(nonnull(lop_));
  shared_jacobian_->getObject(this);
  return lop_;
}


Teuchos::RCP<const ::Thyra::LinearOpBase<double> >
NOX::Thyra::Group::getJacobianOperator() const
{
  NOX_ASSERT(nonnull(lop_));
  return lop_;
}


Teuchos::RCP<const ::Thyra::LinearOpBase<double> >
NOX::Thyra::Group::getScaledJacobianOperator() const
{
  NOX_ASSERT(nonnull(lop_));
  if(rightScalingFirst_){
    const Teuchos::RCP< ::Thyra::ScaledLinearOpBase<double> > lop_scaled =
      Teuchos::rcp_dynamic_cast< ::Thyra::ScaledLinearOpBase<double> >(lop_, true);
    lop_scaled->scaleRight(*right_weight_vec_);
  }
  return lop_;
}


void
NOX::Thyra::Group::unscaleJacobianOperator() const
{
  NOX_ASSERT(nonnull(lop_));
  if(rightScalingFirst_){
    const Teuchos::RCP< ::Thyra::ScaledLinearOpBase<double> > lop_scaled =
      Teuchos::rcp_dynamic_cast< ::Thyra::ScaledLinearOpBase<double> >(lop_, true);
    lop_scaled->scaleRight(*inv_right_weight_vec_);
  }
}


Teuchos::RCP< ::Thyra::LinearOpWithSolveBase<double> >
NOX::Thyra::Group::getNonconstJacobian()
{
  this->updateLOWS();
  return shared_jacobian_->getObject(this);
}


Teuchos::RCP<const ::Thyra::LinearOpWithSolveBase<double> >
NOX::Thyra::Group::getJacobian() const
{
  this->updateLOWS();
  return shared_jacobian_->getObject();
}


Teuchos::RCP< ::Thyra::PreconditionerBase<double> >
NOX::Thyra::Group::getNonconstPreconditioner()
{
  return prec_;
}


Teuchos::RCP<const ::Thyra::PreconditionerBase<double> >
NOX::Thyra::Group::getPreconditioner() const
{
  return prec_;
}

void NOX::Thyra::Group::setJacobianOperator(const Teuchos::RCP<::Thyra::LinearOpBase<double>>& op)
{
  lop_ = op;
}

void NOX::Thyra::Group::setPreconditionerMatrix(const Teuchos::RCP<const ::Thyra::DefaultLinearOpSource<double>>& op)
{
  losb_ = op;
}

void NOX::Thyra::Group::setX(const NOX::Abstract::Vector& y)
{
  setX(dynamic_cast<const NOX::Thyra::Vector&> (y));
}


void NOX::Thyra::Group::setX(const NOX::Thyra::Vector& y)
{
  resetIsValidFlags();
  *x_vec_ = y;

  if (nonnull(right_weight_vec_))
      computeScaledSolution();
}


void NOX::Thyra::Group::computeX(const NOX::Abstract::Group& grp,
                 const NOX::Abstract::Vector& d,
                 double step)
{
  const NOX::Thyra::Group& thyra_grp =
    dynamic_cast<const NOX::Thyra::Group&> (grp);

  const NOX::Thyra::Vector& thyra_d =
    dynamic_cast<const NOX::Thyra::Vector&> (d);

  this->computeX(thyra_grp, thyra_d, step);
}

void NOX::Thyra::Group::computeX(const NOX::Thyra::Group& grp,
                 const NOX::Thyra::Vector& d,
                 double step)
{
  this->resetIsValidFlags();
  x_vec_->update(1.0, *(grp.x_vec_), step, d);

  if (nonnull(right_weight_vec_))
      computeScaledSolution();
}

NOX::Abstract::Group::ReturnType NOX::Thyra::Group::computeF()
{
  if (this->isF())
    return NOX::Abstract::Group::Ok;

  auto in_args = model_->createInArgs();
  auto out_args = model_->createOutArgs();

  if (use_base_point_)
    in_args = base_point_;

  if (use_pseudo_transient_terms_) {
    in_args.set_x_dot(x_dot_);
    in_args.set_alpha(alpha_);
    in_args.set_beta(beta_);
    in_args.set_t(t_);
  }

  in_args.set_x(x_vec_->getThyraRCPVector().assert_not_null());
  out_args.set_f(f_vec_->getThyraRCPVector().assert_not_null());

  //model_->setVerbLevel(Teuchos::VERB_EXTREME);
  model_->evalModel(in_args, out_args);

  is_valid_f_ = true;

  if (out_args.isFailed())
    return NOX::Abstract::Group::Failed;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType NOX::Thyra::Group::computeJacobian()
{
  if (this->isJacobian())
    return NOX::Abstract::Group::Ok;

  NOX_ASSERT(nonnull(lop_));
  NOX_ASSERT(nonnull(shared_jacobian_));

  shared_jacobian_->getObject(this);

  auto in_args = model_->createInArgs();
  auto out_args = model_->createOutArgs();

  if (use_base_point_)
    in_args = base_point_;

  if (use_pseudo_transient_terms_ == true) {
    in_args.set_x_dot(x_dot_);
    in_args.set_alpha(alpha_);
    in_args.set_beta(beta_);
    in_args.set_t(t_);
  }

  in_args.set_x(x_vec_->getThyraRCPVector());
  out_args.set_W_op(lop_);

  model_->evalModel(in_args, out_args);

  is_valid_jacobian_ = true;

  if (out_args.isFailed())
    return NOX::Abstract::Group::Failed;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType NOX::Thyra::Group::computeGradient()
{
  NOX_ASSERT(nonnull(lop_));
  NOX_ASSERT(nonnull(lows_factory_));
  if ( ::Thyra::opSupported(*shared_jacobian_->getObject(), ::Thyra::TRANS) ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true,  std::logic_error,
               "NOX Error - compute gradient not implemented yet!");
    //return NOX::Abstract::Group::Ok;
  }
  return NOX::Abstract::Group::Failed;
}

NOX::Abstract::Group::ReturnType NOX::Thyra::Group::
computeNewton(Teuchos::ParameterList& p)
{
  NOX::Abstract::Group::ReturnType status =
    this->applyJacobianInverse(p, *f_vec_, *newton_vec_);
  newton_vec_->scale(-1.0);

  return status;
}

NOX::Abstract::Group::ReturnType
NOX::Thyra::Group::applyJacobian(const Abstract::Vector& input,
                 NOX::Abstract::Vector& result) const
{
  using Teuchos::dyn_cast;
  return applyJacobian(dyn_cast<const NOX::Thyra::Vector>(input),
               dyn_cast<NOX::Thyra::Vector>(result));

}

NOX::Abstract::Group::ReturnType
NOX::Thyra::Group::applyJacobian(const Vector& input, Vector& result) const
{
  NOX_ASSERT(nonnull(lop_));
  if ( !(this->isJacobian()) ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
               "NOX Error - Jacobian is not valid.  " <<
               "Call computeJacobian before calling applyJacobian!");
  }

  ::Thyra::apply(*lop_, ::Thyra::NOTRANS,
         input.getThyraVector(), result.getThyraRCPVector().ptr());

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
NOX::Thyra::Group::applyJacobianMultiVector(
                 const NOX::Abstract::MultiVector& input,
                 NOX::Abstract::MultiVector& result) const
{
  if ( !(this->isJacobian()) ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
               "NOX Error - Jacobian is not valid.  " <<
               "Call computeJacobian before calling applyJacobian!");
  }

  NOX_ASSERT(nonnull(lop_));

  const NOX::Thyra::MultiVector& nt_input =
    Teuchos::dyn_cast<const NOX::Thyra::MultiVector>(input);
  NOX::Thyra::MultiVector& nt_result =
    Teuchos::dyn_cast<NOX::Thyra::MultiVector>(result);

  ::Thyra::apply(*lop_,
         ::Thyra::NOTRANS,
         *nt_input.getThyraMultiVector(),
         nt_result.getThyraMultiVector().ptr());

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
NOX::Thyra::Group::applyJacobianTranspose(const NOX::Abstract::Vector& input,
                       NOX::Abstract::Vector& result) const
{
  using Teuchos::dyn_cast;
  return applyJacobianTranspose(dyn_cast<const NOX::Thyra::Vector>(input),
                dyn_cast<NOX::Thyra::Vector>(result));
}

NOX::Abstract::Group::ReturnType
NOX::Thyra::Group::applyJacobianTranspose(const NOX::Thyra::Vector& input,
                      NOX::Thyra::Vector& result) const
{
  if ( !(this->isJacobian()) ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
               "NOX Error - Jacobian is not valid.  " <<
               "Call computeJacobian before calling applyJacobian!");
  }

  NOX_ASSERT(nonnull(lop_));
  NOX_ASSERT(nonnull(shared_jacobian_));

  if ( ::Thyra::opSupported(*lop_, ::Thyra::TRANS) ) {
    ::Thyra::apply(*shared_jacobian_->getObject(), ::Thyra::TRANS,
           input.getThyraVector(), result.getThyraRCPVector().ptr());
    return NOX::Abstract::Group::Ok;
  }
  return NOX::Abstract::Group::Failed;
}

NOX::Abstract::Group::ReturnType
NOX::Thyra::Group::applyJacobianTransposeMultiVector(
                 const NOX::Abstract::MultiVector& input,
                 NOX::Abstract::MultiVector& result) const
{
  if ( !(this->isJacobian()) ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
               "NOX Error - Jacobian is not valid.  " <<
               "Call computeJacobian before calling applyJacobian!");
  }

  NOX_ASSERT(nonnull(lop_));
  NOX_ASSERT(nonnull(shared_jacobian_));

  if (! ::Thyra::opSupported(*shared_jacobian_->getObject(), ::Thyra::TRANS) )
    return NOX::Abstract::Group::Failed;

  const NOX::Thyra::MultiVector& nt_input =
    Teuchos::dyn_cast<const NOX::Thyra::MultiVector>(input);
  NOX::Thyra::MultiVector& nt_result =
    Teuchos::dyn_cast<NOX::Thyra::MultiVector>(result);

  ::Thyra::apply(*lop_,
         ::Thyra::TRANS,
         *nt_input.getThyraMultiVector(),
         nt_result.getThyraMultiVector().ptr());

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
NOX::Thyra::Group::applyJacobianInverse(Teuchos::ParameterList& p,
                    const Abstract::Vector& input,
                    NOX::Abstract::Vector& result) const
{
  using Teuchos::dyn_cast;
  return applyJacobianInverse(p, dyn_cast<const NOX::Thyra::Vector>(input),
                  dyn_cast<NOX::Thyra::Vector>(result));
}

NOX::Abstract::Group::ReturnType
NOX::Thyra::Group::applyJacobianInverse(Teuchos::ParameterList& p,
                    const NOX::Thyra::Vector& input,
                    NOX::Thyra::Vector& result) const
{
  return applyJacobianInverseMultiVector( p, input.getThyraVector(),
                      result.getThyraVector() );
}

NOX::Abstract::Group::ReturnType NOX::Thyra::Group::
applyJacobianInverseMultiVector(Teuchos::ParameterList& p,
                const NOX::Abstract::MultiVector& input,
                NOX::Abstract::MultiVector& result) const
{
  const NOX::Thyra::MultiVector& nt_input =
    Teuchos::dyn_cast<const NOX::Thyra::MultiVector>(input);
  NOX::Thyra::MultiVector& nt_result =
    Teuchos::dyn_cast<NOX::Thyra::MultiVector>(result);

  return applyJacobianInverseMultiVector(p,
                     *nt_input.getThyraMultiVector(),
                     *nt_result.getThyraMultiVector());
}

bool NOX::Thyra::Group::isF() const
{
  return is_valid_f_;
}

bool NOX::Thyra::Group::isJacobian() const
{
  NOX_ASSERT(nonnull(shared_jacobian_));
  return ((shared_jacobian_->isOwner(this)) && (is_valid_jacobian_));
}

bool NOX::Thyra::Group::isNewton() const
{
  return is_valid_newton_dir_;
}

bool NOX::Thyra::Group::isGradient() const
{
  return is_valid_gradient_dir_;
}

const NOX::Abstract::Vector& NOX::Thyra::Group::getX() const
{
  return *x_vec_;
}

const NOX::Abstract::Vector& NOX::Thyra::Group::getScaledX() const
{ 
  if(nonnull(inv_right_weight_vec_))
    return *scaled_x_vec_;
  else
    return *x_vec_;
}

const NOX::Abstract::Vector& NOX::Thyra::Group::getF() const
{
  return *f_vec_;
}

double NOX::Thyra::Group::getNormF() const
{
  if ( !(this->isF()) ) {
    std::cerr << "ERROR: NOX::Thyra::Group::getNormF() "
	      << "- F is not up to date.  Please call computeF()!" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return f_vec_->norm();
}

const NOX::Abstract::Vector& NOX::Thyra::Group::getNewton() const
{
  return *newton_vec_;
}

const NOX::Abstract::Vector& NOX::Thyra::Group::getGradient() const
{
  return *gradient_vec_;
}

Teuchos::RCP< const NOX::Abstract::Vector > NOX::Thyra::Group::getXPtr() const
{
  return x_vec_;
}

Teuchos::RCP< const NOX::Abstract::Vector > NOX::Thyra::Group::getFPtr() const
{
  return f_vec_;
}

Teuchos::RCP< const NOX::Abstract::Vector > NOX::Thyra::Group::getNewtonPtr() const
{
  NOX_ASSERT(nonnull(newton_vec_));
  return newton_vec_;
}

Teuchos::RCP< const NOX::Abstract::Vector > NOX::Thyra::Group::getGradientPtr() const
{
  return gradient_vec_;
}

void NOX::Thyra::Group::logLastLinearSolveStats(NOX::SolverStats& stats) const
{
  stats.linearSolve.logLinearSolve(last_linear_solve_status_ == NOX::Abstract::Group::Ok,
                                   last_linear_solve_num_iters_,
                                   last_linear_solve_achieved_tol_,
                                   0.0,
                                   0.0);
}

void NOX::Thyra::Group::print() const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
}

// protected
NOX::Abstract::Group::ReturnType NOX::Thyra::Group::
applyJacobianInverseMultiVector(Teuchos::ParameterList& p,
                const ::Thyra::MultiVectorBase<double>& input,
                ::Thyra::MultiVectorBase<double>& result) const
{
  this->updateLOWS();

  // Create solve criteria
  ::Thyra::SolveCriteria<double> solveCriteria;
  solveCriteria.requestedTol = p.get("Tolerance", 1.0e-6);

  std::string numer_measure = p.get("Solve Measure Numerator",
                    "Norm Residual");
  std::string denom_measure = p.get("Solve Measure Denominator",
                    "Norm Initial Residual");
  solveCriteria.solveMeasureType =
    ::Thyra::SolveMeasureType(getThyraNormType(numer_measure),
                  getThyraNormType(denom_measure));

  // Initialize result to zero to remove possible NaNs
  ::Thyra::assign(Teuchos::ptrFromRef(result), 0.0);

  this->scaleResidualAndJacobian();

  ::Thyra::SolveStatus<double> solve_status;
  {
    NOX_FUNC_TIME_MONITOR("NOX Total Linear Solve");

    solve_status = ::Thyra::solve(*shared_jacobian_->getObject(),
                  ::Thyra::NOTRANS, input,
                  Teuchos::ptrFromRef(result),
                  Teuchos::constPtr(solveCriteria));
  }

  // If the linear solver left us an iteration count, stash it on the output list
  if(!solve_status.extraParameters.is_null()) {
    int current_iters    = solve_status.extraParameters->get("Iteration Count",0);
    int cumulative_iters = p.sublist("Output").get("Cumulative Iteration Count",0);

    p.sublist("Output").set("Last Iteration Count",current_iters);
    p.sublist("Output").set("Cumulative Iteration Count",cumulative_iters + current_iters);
    last_linear_solve_num_iters_ = current_iters;
  }

  this->unscaleResidualAndJacobian();

  if (nonnull(right_weight_vec_)){
  
    const Teuchos::RCP< ::Thyra::ScaledLinearOpBase<double> > result_scaled =
      Teuchos::rcp_dynamic_cast< ::Thyra::ScaledLinearOpBase<double> >(Teuchos::rcpFromRef(result), true);
    result_scaled->scaleLeft(*right_weight_vec_);

  }

  last_linear_solve_achieved_tol_ = solve_status.achievedTol;

  last_linear_solve_status_ = NOX::Abstract::Group::Failed;
  if (solve_status.solveStatus == ::Thyra::SOLVE_STATUS_CONVERGED)
    last_linear_solve_status_ = NOX::Abstract::Group::Ok;
  else if (solve_status.solveStatus == ::Thyra::SOLVE_STATUS_UNCONVERGED)
    last_linear_solve_status_ = NOX::Abstract::Group::NotConverged;

  return last_linear_solve_status_;
}

NOX::Abstract::Group::ReturnType
NOX::Thyra::Group::applyRightPreconditioning(bool /* useTranspose */,
                         Teuchos::ParameterList& /* params */,
                         const NOX::Abstract::Vector& input,
                         NOX::Abstract::Vector& result) const
{
  NOX_ASSERT(nonnull(prec_));

  // A nonnull prec_factory_ means we are using Jacobian for M and use
  // the prec_factory_ to produce M^{-1}.  Otherwise, we assume that
  // M^{-1} is provided directly by the user from the model evaluator.
  // Finally if the model evauator does not support a prec then just
  // use the operator the user supplied prec_ object and assume they
  // know when to update it externally.

  if (nonnull(prec_factory_)) {
    NOX_ASSERT(nonnull(losb_));

    if (!this->isJacobian())
      const_cast<NOX::Thyra::Group*>(this)->computeJacobian();

    this->scaleResidualAndJacobian();
    prec_factory_->initializePrec(losb_, prec_.get());
  }
  else if (model_->createOutArgs().supports( ::Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
    auto in_args = model_->createInArgs();
    auto out_args = model_->createOutArgs();

    if (use_base_point_)
      in_args = base_point_;

    if (use_pseudo_transient_terms_ == true) {
      in_args.set_x_dot(x_dot_);
      in_args.set_alpha(alpha_);
      in_args.set_beta(beta_);
      in_args.set_t(t_);
    }

    in_args.set_x(x_vec_->getThyraRCPVector().assert_not_null());
    out_args.set_W_prec(prec_);

    model_->evalModel(in_args, out_args);
  }

  const NOX::Thyra::Vector & inputThyraVector = dynamic_cast<const NOX::Thyra::Vector&>(input);
  NOX::Thyra::Vector & resultThyraVector = dynamic_cast<NOX::Thyra::Vector&>(result);

  // Could be left, right or unspecified
  Teuchos::RCP<const ::Thyra::LinearOpBase<double> > tmp_prec_ = prec_->getRightPrecOp();
  if (is_null(tmp_prec_))
    tmp_prec_ = prec_->getUnspecifiedPrecOp();
  NOX_ASSERT(nonnull(tmp_prec_));

  ::Thyra::apply(*tmp_prec_,
         ::Thyra::NOTRANS,
         *inputThyraVector.getThyraRCPVector(),
         outArg(*resultThyraVector.getThyraRCPVector().ptr()));

  if (nonnull(prec_factory_))
    this->unscaleResidualAndJacobian();

  return NOX::Abstract::Group::Ok;
}

::Thyra::ESolveMeasureNormType
NOX::Thyra::Group::getThyraNormType(const std::string& name) const
{
  if (name == "None")
    return ::Thyra::SOLVE_MEASURE_ONE;
  else if (name == "Norm Residual")
    return ::Thyra::SOLVE_MEASURE_NORM_RESIDUAL;
  else if (name == "Norm Solution")
    return ::Thyra::SOLVE_MEASURE_NORM_SOLUTION;
  else if (name == "Norm Initial Residual")
    return ::Thyra::SOLVE_MEASURE_NORM_INIT_RESIDUAL;
  else if (name == "Norm RHS")
    return ::Thyra::SOLVE_MEASURE_NORM_RHS;
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,  std::logic_error,
               "NOX Error - unknown solve measure " << name);
  }
  TEUCHOS_UNREACHABLE_RETURN(::Thyra::SOLVE_MEASURE_ONE);
}

void NOX::Thyra::Group::updateLOWS() const
{
  if (is_valid_lows_)
    return;

  NOX_ASSERT(nonnull(lop_));
  NOX_ASSERT(nonnull(lows_factory_));

  this->scaleResidualAndJacobian();

  {
    NOX_FUNC_TIME_MONITOR("NOX Total Preconditioner Construction");

    if (nonnull(prec_) && (!updatePreconditioner_)) {
      // Use the preconditioner. If the matrix values need to be
      // updated, it must be handled by user manually
      ::Thyra::initializePreconditionedOp<double>(*lows_factory_,
                          lop_,
                          prec_,
                          shared_jacobian_->getObject(this).ptr());
    }
    else if (nonnull(prec_factory_) && updatePreconditioner_) {
      // Automatically update using the user supplied prec factory
      prec_factory_->initializePrec(losb_, prec_.get());

      ::Thyra::initializePreconditionedOp<double>(*lows_factory_,
                          lop_,
                          prec_,
                          shared_jacobian_->getObject(this).ptr());
    }
    else if ( nonnull(prec_) && (model_->createOutArgs().supports( ::Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) && updatePreconditioner_) {
      // Automatically update using the ModelEvaluator
      auto in_args = model_->createInArgs();
      auto out_args = model_->createOutArgs();

      if (use_base_point_)
        in_args = base_point_;

      if (use_pseudo_transient_terms_ == true) {
        in_args.set_x_dot(x_dot_);
        in_args.set_alpha(alpha_);
        in_args.set_beta(beta_);
        in_args.set_t(t_);
      }

      in_args.set_x(x_vec_->getThyraRCPVector().assert_not_null());
      out_args.set_W_prec(prec_);

      model_->evalModel(in_args, out_args);

      ::Thyra::initializePreconditionedOp<double>(*lows_factory_,
                          lop_,
                          prec_,
                          shared_jacobian_->getObject(this).ptr());

    }
    else {
      // No preconditioner
      ::Thyra::initializeOp<double>(*lows_factory_,
                    lop_,
                    shared_jacobian_->getObject(this).ptr());
    }

  }

  this->unscaleResidualAndJacobian();

  is_valid_lows_ = true;
}

void NOX::Thyra::Group::scaleResidualAndJacobian() const
{
  if (nonnull(weight_vec_)) {
    ::Thyra::ele_wise_scale(*weight_vec_, f_vec_->getThyraRCPVector().ptr());

    const Teuchos::RCP< ::Thyra::ScaledLinearOpBase<double> > W_scaled =
      Teuchos::rcp_dynamic_cast< ::Thyra::ScaledLinearOpBase<double> >(lop_, true);
    W_scaled->scaleLeft(*weight_vec_);
  }

  if (nonnull(right_weight_vec_)) {
    const Teuchos::RCP< ::Thyra::ScaledLinearOpBase<double> > W_scaled =
      Teuchos::rcp_dynamic_cast< ::Thyra::ScaledLinearOpBase<double> >(lop_, true);
    W_scaled->scaleRight(*right_weight_vec_);
  }
}

void NOX::Thyra::Group::unscaleResidualAndJacobian() const
{
  if (nonnull(weight_vec_)) {

    if (is_null(inv_weight_vec_))
      inv_weight_vec_ = weight_vec_->clone_v();

    // Always recompute inverse scaling.  We don't know when a user
    // might update the scaling vector in his code
    ::Thyra::reciprocal(*weight_vec_, inv_weight_vec_.ptr());

    ::Thyra::ele_wise_scale(*inv_weight_vec_, f_vec_->getThyraRCPVector().ptr());

    const Teuchos::RCP< ::Thyra::ScaledLinearOpBase<double> > W_scaled =
      Teuchos::rcp_dynamic_cast< ::Thyra::ScaledLinearOpBase<double> >(lop_, true);
    W_scaled->scaleLeft(*inv_weight_vec_);
  }

  if (nonnull(right_weight_vec_)) {

    const Teuchos::RCP< ::Thyra::ScaledLinearOpBase<double> > W_scaled =
      Teuchos::rcp_dynamic_cast< ::Thyra::ScaledLinearOpBase<double> >(lop_, true);
    W_scaled->scaleRight(*inv_right_weight_vec_);

  }
}

void NOX::Thyra::Group::computeScaledSolution()
{
  *scaled_x_vec_ = *x_vec_;
  NOX::Thyra::Vector scaling(inv_right_weight_vec_);
  scaled_x_vec_->scale(scaling);  
}

Teuchos::RCP< const ::Thyra::ModelEvaluator<double> >
NOX::Thyra::Group::getModel() const
{
  return model_;
}

void NOX::Thyra::Group::enablePseudoTransientTerms(const Teuchos::RCP<const ::Thyra::VectorBase<double>>& x_dot,
                                                   const double alpha, const double beta, const double t)
{
  use_pseudo_transient_terms_ = true;
  x_dot_ = x_dot;
  alpha_ = alpha;
  beta_ = beta;
  t_ = t;
}

void NOX::Thyra::Group::disablePseudoTransientTerms()
{
  use_pseudo_transient_terms_ = false;
  x_dot_ = Teuchos::null;
  alpha_ = 0.0;
  beta_ = 1.0;
  t_ = -1.0;
}

bool NOX::Thyra::Group::usingPseudoTransientTerms() const
{ return use_pseudo_transient_terms_; }

void NOX::Thyra::Group::setBasePoint(const ::Thyra::ModelEvaluatorBase::InArgs<double>& base_point_params)
{
  use_base_point_ = true;
  base_point_ = base_point_params;
}

void NOX::Thyra::Group::unsetBasePoint()
{
  use_base_point_ = false;
  // Set to default to reset RCPs to free memory
  base_point_ = model_->createInArgs();
}

bool NOX::Thyra::Group::usingBasePoint() const
{ return use_base_point_; }

Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double>>
NOX::Thyra::Group::getLinearOpWithSolveFactory() const
{ return lows_factory_; }

Teuchos::RCP<::Thyra::PreconditionerFactoryBase<double>>
NOX::Thyra::Group::getPreconditionerFactory() const
{ return prec_factory_; }

void
NOX::Thyra::Group::
takeControlOfPreconditionerUpdates(const Teuchos::RCP< ::Thyra::PreconditionerBase<double>>& prec)
{
  prec_ = prec;
  prec_factory_ = Teuchos::null;
  updatePreconditioner_ = false;

  // Unset factory so solver doesn't automatically update preconditioner
  Teuchos::rcp_const_cast<::Thyra::LinearOpWithSolveFactoryBase<double>>(lows_factory_)->unsetPreconditionerFactory(nullptr,nullptr);
}
