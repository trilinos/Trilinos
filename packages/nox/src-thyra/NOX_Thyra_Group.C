// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "Teuchos_Assert.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_TimeMonitor.hpp"
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
#include "NOX_Thyra_Group.H"	// class definition
#include "NOX_Abstract_MultiVector.H"
#include "NOX_Thyra_MultiVector.H"

NOX::Thyra::Group::
Group(const NOX::Thyra::Vector& initial_guess,
      const Teuchos::RCP< const ::Thyra::ModelEvaluator<double> >& model,
      const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weight_vector):
  model_(model)
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

  f_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));
  newton_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));
  gradient_vec_ = 
    Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));
  
  lop_ = model->create_W_op();

  // create jacobian operator
  lows_factory_ = model->get_W_factory();

  // Create jacobian with solver
  shared_jacobian_ = Teuchos::rcp(new NOX::SharedObject< ::Thyra::LinearOpWithSolveBase<double>, NOX::Thyra::Group >(lows_factory_->createOp()));

  losb_ = Teuchos::rcp(new ::Thyra::DefaultLinearOpSource<double>(lop_));

  // create preconditioner
  prec_factory_ = lows_factory_->getPreconditionerFactory();
  
  if (Teuchos::nonnull(prec_factory_))
    prec_ = prec_factory_->createPrec();

  // Create in/out args
  in_args_ = model_->createInArgs();
  out_args_ = model_->createOutArgs();

  resetIsValidFlags();
}

NOX::Thyra::Group::
Group(const NOX::Thyra::Vector& initial_guess,
      const Teuchos::RCP< const ::Thyra::ModelEvaluator<double> >& model,
      const Teuchos::RCP< ::Thyra::LinearOpBase<double> >& linear_op,
      const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> >& lows_factory,
      const Teuchos::RCP< ::Thyra::PreconditionerBase<double> >& prec_op,
      const Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double> >& prec_factory,
      const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weight_vector):
  model_(model),
  lop_(linear_op),
  lows_factory_(lows_factory),
  prec_(prec_op),
  prec_factory_(prec_factory)
{
  TEUCHOS_ASSERT(nonnull(lop_));
  TEUCHOS_ASSERT(nonnull(lows_factory_));

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

  f_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));
  newton_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));
  gradient_vec_ = 
    Teuchos::rcp(new NOX::Thyra::Vector(*x_vec_, ShapeCopy));
  
  // Create jacobian with solver
  shared_jacobian_ = Teuchos::rcp(new NOX::SharedObject< ::Thyra::LinearOpWithSolveBase<double>, NOX::Thyra::Group >(lows_factory_->createOp()));

  losb_ = Teuchos::rcp(new ::Thyra::DefaultLinearOpSource<double>(lop_));

  if ( nonnull(prec_factory_) && is_null(prec_) )
    prec_ = prec_factory_->createPrec();

  // Create in/out args
  in_args_ = model_->createInArgs();
  out_args_ = model_->createOutArgs();

  resetIsValidFlags();
}

NOX::Thyra::Group::Group(const NOX::Thyra::Group& source, NOX::CopyType type) :
  model_(source.model_),
  shared_jacobian_(source.shared_jacobian_),
  lop_(source.lop_),
  lows_factory_(source.lows_factory_),
  losb_(source.losb_),
  prec_(source.prec_),
  prec_factory_(source.prec_factory_)
{

  x_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*source.x_vec_, type));
  f_vec_ = Teuchos::rcp(new NOX::Thyra::Vector(*source.f_vec_, type));
  newton_vec_ = 
    Teuchos::rcp(new NOX::Thyra::Vector(*source.newton_vec_, type));
  gradient_vec_ = 
    Teuchos::rcp(new NOX::Thyra::Vector(*source.gradient_vec_, type));

  if (nonnull(source.weight_vec_))
    weight_vec_ = source.weight_vec_;

  in_args_ = model_->createInArgs();
  out_args_ = model_->createOutArgs();
  
  if (type == NOX::DeepCopy) {
    is_valid_f_ = source.is_valid_f_;
    is_valid_jacobian_ = source.is_valid_jacobian_;
    is_valid_newton_dir_ = source.is_valid_newton_dir_;
    is_valid_gradient_dir_ = source.is_valid_gradient_dir_;
    is_valid_lows_ = source.is_valid_lows_;

    // New copy takes ownership of the shared Jacobian for DeepCopy
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

  // If valid, this takes ownership of the shared Jacobian
  if (this->isJacobian())
    shared_jacobian_->getObject(this);

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
  shared_jacobian_->getObject(this);
  return lop_;
}


Teuchos::RCP<const ::Thyra::LinearOpBase<double> >
NOX::Thyra::Group::getJacobianOperator() const
{
  return lop_;
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


void NOX::Thyra::Group::setX(const NOX::Abstract::Vector& y) 
{
  setX(dynamic_cast<const NOX::Thyra::Vector&> (y));
}


void NOX::Thyra::Group::setX(const NOX::Thyra::Vector& y) 
{
  resetIsValidFlags();
  *x_vec_ = y;
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
}

NOX::Abstract::Group::ReturnType NOX::Thyra::Group::computeF() 
{
  if (this->isF())
    return NOX::Abstract::Group::Ok;

  in_args_.set_x(x_vec_->getThyraRCPVector().assert_not_null());
  out_args_.set_f(f_vec_->getThyraRCPVector().assert_not_null());
  //model_->setVerbLevel(Teuchos::VERB_EXTREME);
  model_->evalModel(in_args_, out_args_);
  in_args_.set_x(Teuchos::null);
  out_args_.set_f(Teuchos::null);

  is_valid_f_ = true;

  if (out_args_.isFailed())
    return NOX::Abstract::Group::Failed;
  
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType NOX::Thyra::Group::computeJacobian() 
{
  if (this->isJacobian())
    return NOX::Abstract::Group::Ok;

  shared_jacobian_->getObject(this);

  in_args_.set_x(x_vec_->getThyraRCPVector());
  out_args_.set_W_op(lop_);
  model_->evalModel(in_args_, out_args_);
  in_args_.set_x(Teuchos::null);
  out_args_.set_W_op(Teuchos::null);

  is_valid_jacobian_ = true;

  if (out_args_.isFailed())
    return NOX::Abstract::Group::Failed;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType NOX::Thyra::Group::computeGradient() 
{
  if ( ::Thyra::opSupported(*shared_jacobian_->getObject(), ::Thyra::TRANS) ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true,  std::logic_error, 
		       "NOX Error - compute gradient not implemented yet!");
    return NOX::Abstract::Group::Ok;
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
  if ( !(this->isJacobian()) ) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
		       "NOX Error - Jacobian is not valid.  " <<
		       "Call computeJacobian before calling applyJacobian!");
  }
  
  //::Thyra::apply(*shared_jacobian_->getObject(), ::Thyra::NOTRANS,
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

const NOX::Abstract::Vector& NOX::Thyra::Group::getF() const 
{  
  return *f_vec_;
}

double NOX::Thyra::Group::getNormF() const
{
  if ( this->isF() )
    return f_vec_->norm();

  cerr << "ERROR: NOX::Thyra::Group::getNormF() "
       << "- F is not up to date.  Please call computeF()!" << endl;
  throw "NOX Error";
  
  return 0.0;
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
  return newton_vec_;
}

Teuchos::RCP< const NOX::Abstract::Vector > NOX::Thyra::Group::getGradientPtr() const 
{ 
  return gradient_vec_;
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

  this->unscaleResidualAndJacobian();

  // ToDo: Get the output statistics and achieved tolerance to pass
  // back ...
  
  if (solve_status.solveStatus == ::Thyra::SOLVE_STATUS_CONVERGED)
    return NOX::Abstract::Group::Ok;
  else if (solve_status.solveStatus == ::Thyra::SOLVE_STATUS_UNCONVERGED)
    return NOX::Abstract::Group::NotConverged;
  
  return NOX::Abstract::Group::Failed;
}

::Thyra::ESolveMeasureNormType 
NOX::Thyra::Group::getThyraNormType(const string& name) const
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
    return ::Thyra::SOLVE_MEASURE_ONE;
  }
}

void NOX::Thyra::Group::updateLOWS() const
{ 
  if (is_valid_lows_)
    return;

  this->scaleResidualAndJacobian();

  {
    NOX_FUNC_TIME_MONITOR("NOX Total Preconditioner Construction");
    
    if (nonnull(prec_factory_)) {
      prec_factory_->initializePrec(losb_, prec_.get());
      
      ::Thyra::initializePreconditionedOp<double>(*lows_factory_,
						  lop_,
						  prec_,
						  shared_jacobian_->getObject(this).ptr());
    }
    else if ( nonnull(prec_) && (out_args_.supports( ::Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) ) {
      in_args_.set_x(x_vec_->getThyraRCPVector().assert_not_null());
      out_args_.set_W_prec(prec_);
      model_->evalModel(in_args_, out_args_);
      in_args_.set_x(Teuchos::null);
      out_args_.set_W_prec(Teuchos::null);

      ::Thyra::initializePreconditionedOp<double>(*lows_factory_,
						  lop_,
						  prec_,
						  shared_jacobian_->getObject(this).ptr());
      
    }
    else {
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
}

Teuchos::RCP< const ::Thyra::ModelEvaluator<double> > 
NOX::Thyra::Group::getModel() const
{
  return model_;
}

::Thyra::ModelEvaluatorBase::InArgs<double>& 
NOX::Thyra::Group::getNonconstInArgs()
{
  return in_args_;
}

const ::Thyra::ModelEvaluatorBase::InArgs<double>& 
NOX::Thyra::Group::getInArgs() const
{
  return in_args_;
}
