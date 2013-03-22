// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#include "LOCA_Thyra_Group.H"	          // class definition
#include "NOX_Thyra_MultiVector.H"
#include "Teuchos_Assert.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Thyra::Group::Group(
	    const Teuchos::RCP<LOCA::GlobalData>& global_data,
	    const NOX::Thyra::Vector& initial_guess,
	    const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model,
	    const LOCA::ParameterVector& p,
	    int p_index,
	    bool impl_dfdp) :
  NOX::Thyra::Group(initial_guess, model),
  LOCA::Abstract::Group(global_data),
  globalData(global_data),
  params(p),
  param_index(p_index),
  saveDataStrategy(),
  implement_dfdp(impl_dfdp)
{
  // Create thyra vector to store parameters that is a view of LOCA
  // parameter vector
  RTOpPack::SubVectorView<double> pv(0,params.length(),
				     Teuchos::ArrayRCP<double>(params.getDoubleArrayPointer(),0,params.length(),false),1);
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > ps = 
    model_->get_p_space(param_index);
  param_thyra_vec = ::Thyra::createMemberView(ps,pv);

  // Create x_dot vector of zeros
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > xs = 
    model_->get_x_space();
  x_dot_vec = ::Thyra::createMember(xs);
  ::Thyra::put_scalar(0.0, x_dot_vec.ptr());
}

LOCA::Thyra::Group::Group(const LOCA::Thyra::Group& source, 
			   NOX::CopyType type) :
  NOX::Thyra::Group(source, type),
  LOCA::Abstract::Group(source, type),
  globalData(source.globalData),
  params(source.params),
  param_index(source.param_index),
  saveDataStrategy(source.saveDataStrategy),
  implement_dfdp(source.implement_dfdp)
{
  // Create thyra vector to store parameters that is a view of LOCA
  // parameter vector
  RTOpPack::SubVectorView<double> pv(0,params.length(),
				     Teuchos::ArrayRCP<double>(params.getDoubleArrayPointer(),0,params.length(),false),1);
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > ps = 
    model_->get_p_space(param_index);
  param_thyra_vec = ::Thyra::createMemberView(ps,pv);

  // Create x_dot vector of zeros
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > xs = 
    model_->get_x_space();
  x_dot_vec = ::Thyra::createMember(xs);
  ::Thyra::put_scalar(0.0, x_dot_vec.ptr());
}

LOCA::Thyra::Group::~Group() 
{
}

LOCA::Thyra::Group& 
LOCA::Thyra::Group::operator=(const LOCA::Thyra::Group& source)
{
  if (this != &source) {
    NOX::Thyra::Group::operator=(source);
    LOCA::Abstract::Group::copy(source);
    params = source.params;
    param_index = source.param_index;
    saveDataStrategy = source.saveDataStrategy;
    implement_dfdp = source.implement_dfdp;

    // Because param_thyra_vec is a view of params, we don't need to copy
  }
  return *this;
}

NOX::Abstract::Group& 
LOCA::Thyra::Group::operator=(const NOX::Abstract::Group& source)
{
  operator=(dynamic_cast<const Group&> (source));
  return *this;
}

NOX::Abstract::Group& 
LOCA::Thyra::Group::operator=(const NOX::Thyra::Group& source)
{
  operator=(dynamic_cast<const Group&> (source));
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::Thyra::Group::clone(NOX::CopyType type) const 
{
  return Teuchos::rcp(new LOCA::Thyra::Group(*this, type));
}

NOX::Abstract::Group::ReturnType
LOCA::Thyra::Group::computeF() 
{
  if (this->isF())
    return NOX::Abstract::Group::Ok;

  in_args_.set_x(x_vec_->getThyraRCPVector().assert_not_null());
  if (in_args_.supports(::Thyra::ModelEvaluatorBase::IN_ARG_x_dot))
    in_args_.set_x_dot(x_dot_vec);
  in_args_.set_p(param_index, param_thyra_vec);
  out_args_.set_f(f_vec_->getThyraRCPVector().assert_not_null());

  model_->evalModel(in_args_, out_args_);

  in_args_.set_x(Teuchos::null);
  in_args_.set_p(param_index, Teuchos::null);
  out_args_.set_f(Teuchos::null);

  is_valid_f_ = true;

  if (out_args_.isFailed())
    return NOX::Abstract::Group::Failed;
  
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::Thyra::Group::computeJacobian() 
{
  if (this->isJacobian())
    return NOX::Abstract::Group::Ok;

  shared_jacobian_->getObject(this);

  in_args_.set_x(x_vec_->getThyraRCPVector());
  if (in_args_.supports(::Thyra::ModelEvaluatorBase::IN_ARG_x_dot))
    in_args_.set_x_dot(x_dot_vec);
  if (in_args_.supports(::Thyra::ModelEvaluatorBase::IN_ARG_alpha))
    in_args_.set_alpha(0.0);
  if (in_args_.supports(::Thyra::ModelEvaluatorBase::IN_ARG_beta))
    in_args_.set_beta(1.0);
  in_args_.set_p(param_index, param_thyra_vec);
  out_args_.set_W_op(lop_);

  model_->evalModel(in_args_, out_args_);

  in_args_.set_x(Teuchos::null);
  in_args_.set_p(param_index, Teuchos::null);
  out_args_.set_W_op(Teuchos::null);

  is_valid_jacobian_ = true;

  if (out_args_.isFailed())
    return NOX::Abstract::Group::Failed;

  return NOX::Abstract::Group::Ok;
}

void
LOCA::Thyra::Group::copy(const NOX::Abstract::Group& source)
{
  *this = source;
}

void 
LOCA::Thyra::Group::setParams(const LOCA::ParameterVector& p)
{
  this->resetIsValidFlags();
  params = p;
}

void
LOCA::Thyra::Group::setParam(int paramID, double val)
{
  this->resetIsValidFlags();
  params.setValue(paramID, val);
}

void
LOCA::Thyra::Group::setParam(std::string paramID, double val)
{
  this->resetIsValidFlags();
  params.setValue(paramID, val);
}

const LOCA::ParameterVector& 
LOCA::Thyra::Group::getParams() const 
{
  return params;
}

double
LOCA::Thyra::Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

double
LOCA::Thyra::Group::getParam(std::string paramID) const
{
  return params.getValue(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Thyra::Group::computeDfDpMulti(const std::vector<int>& paramIDs, 
				     NOX::Abstract::MultiVector& fdfdp, 
				     bool isValidF)
{
  // Currently this does not work because the thyra modelevaluator is not
  // setting the parameter names correctly in the epetraext modelevalator, 
  // so we are disabling this for now
  implement_dfdp = false;

  // Use default implementation if we don't want to use model evaluator, or
  // it doesn't support it
  if (!implement_dfdp || 
      !out_args_.supports(::Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,
			  param_index).supports(::Thyra::ModelEvaluatorBase::DERIV_MV_BY_COL)) {
    NOX::Abstract::Group::ReturnType res = 
      LOCA::Abstract::Group::computeDfDpMulti(paramIDs, fdfdp, isValidF);
    return res;
  }

  // Split fdfdp into f and df/dp
  int num_vecs = fdfdp.numVectors()-1;
  std::vector<int> index_dfdp(num_vecs);
  for (int i=0; i<num_vecs; i++)
    index_dfdp[i] = i+1;
  NOX::Thyra::Vector& f = dynamic_cast<NOX::Thyra::Vector&>(fdfdp[0]);
  Teuchos::RCP<NOX::Abstract::MultiVector> dfdp =
    fdfdp.subView(index_dfdp);

  // Right now this isn't very efficient because we have to compute
  // derivatives with respect to all of the parameters, not just
  // paramIDs.  Will have to work out with Ross how to selectively get
  // parameter derivatives
  int np = params.length();
  Teuchos::RCP<NOX::Thyra::MultiVector> dfdp_full =
    Teuchos::rcp_dynamic_cast<NOX::Thyra::MultiVector>(dfdp->clone(np));

  ::Thyra::ModelEvaluatorBase::DerivativeMultiVector<double> dmv(dfdp_full->getThyraMultiVector(), ::Thyra::ModelEvaluatorBase::DERIV_MV_BY_COL);
  ::Thyra::ModelEvaluatorBase::Derivative<double> deriv(dmv);

  in_args_.set_x(x_vec_->getThyraRCPVector().assert_not_null());
  if (in_args_.supports(::Thyra::ModelEvaluatorBase::IN_ARG_x_dot))
    in_args_.set_x_dot(x_dot_vec);
  in_args_.set_p(param_index, param_thyra_vec);
  if (!isValidF)
    out_args_.set_f(f.getThyraRCPVector().assert_not_null());
  out_args_.set_DfDp(param_index, deriv);

  // Evaluate model
  model_->evalModel(in_args_, out_args_);

  // Copy back dfdp
  for (int i=0; i<num_vecs; i++)
    (*dfdp)[i] = (*dfdp_full)[paramIDs[i]];

  // Reset inargs/outargs
  in_args_.set_x(Teuchos::null);
  in_args_.set_p(param_index, Teuchos::null);
  out_args_.set_f(Teuchos::null);
  out_args_.set_DfDp(param_index, 
		     ::Thyra::ModelEvaluatorBase::Derivative<double>());

  if (out_args_.isFailed())
    return NOX::Abstract::Group::Failed;
  
  return NOX::Abstract::Group::Ok;
}

void
LOCA::Thyra::Group::preProcessContinuationStep(
			     LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (saveDataStrategy != Teuchos::null)
    saveDataStrategy->preProcessContinuationStep(stepStatus);
}

void
LOCA::Thyra::Group::postProcessContinuationStep(
			     LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (saveDataStrategy != Teuchos::null)
    saveDataStrategy->postProcessContinuationStep(stepStatus);
}

void
LOCA::Thyra::Group::projectToDraw(const NOX::Abstract::Vector& x,
				  double *px) const
{
  if (saveDataStrategy != Teuchos::null)
    saveDataStrategy->projectToDraw(x, px);
}

int
LOCA::Thyra::Group::projectToDrawDimension() const
{
  if (saveDataStrategy != Teuchos::null)
    return saveDataStrategy->projectToDrawDimension();
  return 0;
}

double
LOCA::Thyra::Group::computeScaledDotProduct(
				       const NOX::Abstract::Vector& a,
				       const NOX::Abstract::Vector& b) const
{
  return a.innerProduct(b) / a.length();
}

void
LOCA::Thyra::Group::printSolution(const double conParam) const
{
  printSolution(*x_vec_, conParam);
}

void
LOCA::Thyra::Group::printSolution(const NOX::Abstract::Vector& x_,
				  const double conParam) const
{
  if (saveDataStrategy != Teuchos::null)
    saveDataStrategy->saveSolution(x_, conParam);
}

void
LOCA::Thyra::Group::scaleVector(NOX::Abstract::Vector& x) const
{
  x.scale(1.0 / sqrt(static_cast<double>(x.length())));
}

NOX::Abstract::Group::ReturnType
LOCA::Thyra::Group::computeShiftedMatrix(double alpha, double beta)
{
  shared_jacobian_->getObject(this);
  
  in_args_.set_x(x_vec_->getThyraRCPVector());
  if (in_args_.supports(::Thyra::ModelEvaluatorBase::IN_ARG_x_dot))
    in_args_.set_x_dot(x_dot_vec);
  in_args_.set_p(param_index, param_thyra_vec);
  in_args_.set_alpha(-beta);
  in_args_.set_beta(alpha);
  out_args_.set_W_op(lop_);

  model_->evalModel(in_args_, out_args_);

  in_args_.set_x(Teuchos::null);
  in_args_.set_p(param_index, Teuchos::null);
  in_args_.set_alpha(0.0);
  in_args_.set_beta(1.0);
  out_args_.set_W_op(Teuchos::null);

  is_valid_jacobian_ = false;
  is_valid_lows_ = false;

  if (out_args_.isFailed())
    return NOX::Abstract::Group::Failed;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::Thyra::Group::applyShiftedMatrix(const NOX::Abstract::Vector& input,
                                        NOX::Abstract::Vector& result) const
{
  const NOX::Thyra::Vector& thyra_input = 
    dynamic_cast<const NOX::Thyra::Vector&>(input);
  NOX::Thyra::Vector& thyra_result = 
    dynamic_cast<NOX::Thyra::Vector&>(result);

  ::Thyra::apply(*lop_, ::Thyra::NOTRANS,
		 thyra_input.getThyraVector(), thyra_result.getThyraRCPVector().ptr());

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::Thyra::Group::applyShiftedMatrixMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const
{
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
LOCA::Thyra::Group::applyShiftedMatrixInverseMultiVector(
			        Teuchos::ParameterList& lsParams, 
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const
{
  return this->applyJacobianInverseMultiVector(lsParams, input, result);
}

void
LOCA::Thyra::Group::setSaveDataStrategy(
			 const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy>& s)
{
  saveDataStrategy = s;
}
