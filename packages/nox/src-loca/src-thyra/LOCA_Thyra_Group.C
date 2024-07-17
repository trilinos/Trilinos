// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Thyra_Group.H"              // class definition
#include "NOX_Thyra_MultiVector.H"
#include "Teuchos_Assert.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "RTOpPack_Types.hpp"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Thyra::Group::Group(
        const Teuchos::RCP<LOCA::GlobalData>& global_data,
        const NOX::Thyra::Vector& initial_guess,
        const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model,
        const LOCA::ParameterVector& p,
        int p_index,
        bool impl_dfdp,
        const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weight_vector,
        const bool set_transient_in_args) :
  NOX::Thyra::Group(initial_guess, model, weight_vector),
  LOCA::Abstract::Group(global_data),
  globalData(global_data),
  params(p),
  param_index(1,p_index),
  saveDataStrategy(),
  implement_dfdp(impl_dfdp),
  weight_vec_(weight_vector),
  paramsInSeparatePVecs(false),
  set_transient_in_args_(set_transient_in_args)
{
  updateThyraParamView();
  updateThyraXDot();
}

LOCA::Thyra::Group::Group(
        const Teuchos::RCP<LOCA::GlobalData>& global_data,
        const NOX::Thyra::Group& nox_group,
        const LOCA::ParameterVector& p,
        const std::vector<int>& p_index,
        bool impl_dfdp,
        const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weight_vector,
        const bool set_transient_in_args) :
  NOX::Thyra::Group(nox_group, NOX::DeepCopy),
  LOCA::Abstract::Group(global_data),
  globalData(global_data),
  params(p),
  param_index(p_index),
  saveDataStrategy(),
  implement_dfdp(impl_dfdp),
  weight_vec_(weight_vector),
  paramsInSeparatePVecs(true),
  set_transient_in_args_(set_transient_in_args)
{
  updateThyraParamView();
  updateThyraXDot();
}

LOCA::Thyra::Group::Group(const LOCA::Thyra::Group& source,
               NOX::CopyType type) :
  NOX::Thyra::Group(source, type),
  LOCA::Abstract::Group(source, type),
  globalData(source.globalData),
  params(source.params),
  param_index(source.param_index),
  saveDataStrategy(source.saveDataStrategy),
  implement_dfdp(source.implement_dfdp),
  paramsInSeparatePVecs(source.paramsInSeparatePVecs),
  set_transient_in_args_(source.set_transient_in_args_)
{
  updateThyraParamView();
  updateThyraXDot();
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
    paramsInSeparatePVecs = source.paramsInSeparatePVecs;
    set_transient_in_args_ = source.set_transient_in_args_;
    updateThyraParamView();
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

  auto in_args = model_->createInArgs();
  auto out_args = model_->createOutArgs();

  if (this->usingBasePoint())
    in_args = this->base_point_;

  in_args.set_x(x_vec_->getThyraRCPVector().assert_not_null());
  if (in_args.supports(::Thyra::ModelEvaluatorBase::IN_ARG_x_dot) && set_transient_in_args_)
    in_args.set_x_dot(x_dot_vec);
  for (size_t i=0; i < param_index.size(); ++i)
    in_args.set_p(param_index[i], param_thyra_vec[i]);
  out_args.set_f(f_vec_->getThyraRCPVector().assert_not_null());

  model_->evalModel(in_args, out_args);

  is_valid_f_ = true;

  if (out_args.isFailed())
    return NOX::Abstract::Group::Failed;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::Thyra::Group::computeJacobian()
{
  if (this->isJacobian())
    return NOX::Abstract::Group::Ok;

  shared_jacobian_->getObject(this);

  auto in_args = model_->createInArgs();
  auto out_args = model_->createOutArgs();

  if (this->usingBasePoint())
    in_args = this->base_point_;

  in_args.set_x(x_vec_->getThyraRCPVector());
  if (set_transient_in_args_) {
    if (in_args.supports(::Thyra::ModelEvaluatorBase::IN_ARG_x_dot))
      in_args.set_x_dot(x_dot_vec);
    if (in_args.supports(::Thyra::ModelEvaluatorBase::IN_ARG_alpha))
      in_args.set_alpha(0.0);
    if (in_args.supports(::Thyra::ModelEvaluatorBase::IN_ARG_beta))
      in_args.set_beta(1.0);
  }
  for (size_t i=0; i < param_index.size(); ++i)
    in_args.set_p(param_index[i], param_thyra_vec[i]);
  out_args.set_W_op(lop_);

  model_->evalModel(in_args, out_args);

  is_valid_jacobian_ = true;

  if (out_args.isFailed())
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
  updateThyraParamView();
}

void
LOCA::Thyra::Group::setParam(int paramID, double val)
{
  this->resetIsValidFlags();
  params.setValue(paramID, val);
  updateThyraParamView();
}

void
LOCA::Thyra::Group::setParam(std::string paramID, double val)
{
  this->resetIsValidFlags();
  params.setValue(paramID, val);
  updateThyraParamView();
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

  auto in_args = model_->createInArgs();
  auto out_args = model_->createOutArgs();

  if (this->usingBasePoint())
    in_args = this->base_point_;

  // Use default implementation if we don't want to use model evaluator, or
  // it doesn't support it
  if (!implement_dfdp ||
      !out_args.supports(::Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,
                         param_index[0]).supports(::Thyra::ModelEvaluatorBase::DERIV_MV_BY_COL)) {
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

  in_args.set_x(x_vec_->getThyraRCPVector().assert_not_null());
  if (in_args.supports(::Thyra::ModelEvaluatorBase::IN_ARG_x_dot) && set_transient_in_args_)
    in_args.set_x_dot(x_dot_vec);
  for (size_t i=0; i < param_index.size(); ++i)
    in_args.set_p(param_index[i], param_thyra_vec[i]);
  if (!isValidF)
    out_args.set_f(f.getThyraRCPVector().assert_not_null());
  out_args.set_DfDp(param_index[0], deriv);

  // Evaluate model
  model_->evalModel(in_args, out_args);

  // Copy back dfdp
  for (int i=0; i<num_vecs; i++)
    (*dfdp)[i] = (*dfdp_full)[paramIDs[i]];

  if (out_args.isFailed())
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
  // Mimic
  //     LOCA::Epetra::ModelEvaluatorInterface::postProcessContinuationStep.
  // Previously, response functions were not evaluated after each LOCA step. I
  // had initially tried to fix this problem using
  //     Piro::LOCASolver::evalConvergedModel.
  // However, that evaluates the tangent and gradient as well as the
  // response. Additionally, sequencing the response function and observer calls
  // is tricky using that approach.
  //   If there are no responses, then we don't have to call evalModel.
  if (model_->Ng() > 0) {
    auto in_args = model_->createInArgs();
    auto out_args = model_->createOutArgs();

    if (this->usingBasePoint())
      in_args = this->base_point_;

    in_args.set_x(x_vec_->getThyraRCPVector().assert_not_null());
    if (in_args.supports(::Thyra::ModelEvaluatorBase::IN_ARG_x_dot) && set_transient_in_args_)
      in_args.set_x_dot(x_dot_vec);
    for (size_t i=0; i < param_index.size(); ++i)
      in_args.set_p(param_index[i], param_thyra_vec[i]);
    out_args.set_f(f_vec_->getThyraRCPVector().assert_not_null());
    // This is the key part. It makes the model evaluator call the response
    // functions.
    const Teuchos::RCP< ::Thyra::VectorBase<double> >
      g0 = ::Thyra::createMember(model_->get_g_space(0));
    out_args.set_g(0, g0);

    model_->evalModel(in_args, out_args);
  }

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

  auto in_args = model_->createInArgs();
  auto out_args = model_->createOutArgs();

  if (this->usingBasePoint())
    in_args = this->base_point_;

  in_args.set_x(x_vec_->getThyraRCPVector());
  // Don't use the set_transient_in_args_ flag here. Shifted Matrix
  // need special flags.
  if (in_args.supports(::Thyra::ModelEvaluatorBase::IN_ARG_x_dot))
    in_args.set_x_dot(x_dot_vec);
  for (size_t i=0; i < param_index.size(); ++i)
    in_args.set_p(param_index[i], param_thyra_vec[i]);
  if (in_args.supports(::Thyra::ModelEvaluatorBase::IN_ARG_alpha))
    in_args.set_alpha(-beta);
  if (in_args.supports(::Thyra::ModelEvaluatorBase::IN_ARG_beta))
    in_args.set_beta(alpha);
  out_args.set_W_op(lop_);

  model_->evalModel(in_args, out_args);

  is_valid_jacobian_ = false;
  is_valid_lows_ = false;

  if (out_args.isFailed())
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

void
LOCA::Thyra::Group::updateThyraParamView()
{
  if (paramsInSeparatePVecs) {
    param_thyra_vec.resize(param_index.size());
    for (size_t i=0; i < param_thyra_vec.size(); ++i) {
      const RTOpPack::ConstSubVectorView<double> pv(Teuchos::ArrayRCP<const double>(params.getDoubleArrayPointer()+i,0,1,false));
      Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > ps =
        model_->get_p_space(param_index[i]);
      param_thyra_vec[i] = ::Thyra::createMemberView(ps,pv);
    }
  }
  else {
    // Create thyra vector to store parameters that is a view of LOCA
    // parameter vector
    param_thyra_vec.resize(1);
    const RTOpPack::ConstSubVectorView<double> pv(Teuchos::ArrayRCP<const double>(params.getDoubleArrayPointer(),0,params.length(),false));
    Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > ps =
      model_->get_p_space(param_index[0]);
    param_thyra_vec[0] = ::Thyra::createMemberView(ps,pv);
  }
}

void
LOCA::Thyra::Group::updateThyraXDot()
{
  // Create x_dot vector of zeros
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > xs =
    model_->get_x_space();
  const Teuchos::RCP< ::Thyra::VectorBase<double> > x_dot_vec_setup =
    ::Thyra::createMember(xs);
  ::Thyra::put_scalar(0.0, x_dot_vec_setup.ptr());
  x_dot_vec = x_dot_vec_setup;
}
