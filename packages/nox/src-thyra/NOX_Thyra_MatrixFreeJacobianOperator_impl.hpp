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

#include "NOX_Abstract_Group.H"
#include "NOX_Utils.H"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

template<typename Scalar>
MatrixFreeJacobianOperator<Scalar>::
MatrixFreeJacobianOperator(Teuchos::ParameterList& printParams) :
  setup_called_(false),
  utils_(printParams)
{ }

template<typename Scalar>
void MatrixFreeJacobianOperator<Scalar>::
setBaseEvaluationToRawThyra(const Teuchos::RCP<const ::Thyra::VectorBase<Scalar> >& x_base,
			    const Teuchos::RCP<const ::Thyra::VectorBase<Scalar> >& f_base,
			    const Teuchos::RCP< ::Thyra::ModelEvaluator<Scalar> > model)
{
  TEUCHOS_ASSERT(setup_called_);
  TEUCHOS_ASSERT(nonnull(x_base));
  TEUCHOS_ASSERT(nonnull(f_base));
  TEUCHOS_ASSERT(nonnull(model));
  x_base_ = x_base;
  f_base_ = f_base;
  x_perturb_ = x_base_->space()->createMember();
  f_perturb_ = f_base_->space()->createMember();

  if (difference_type_ == Centered)
    f2_perturb_ = f_base_->space()->createMember();
    
  model_ = model;
  base_evaluation_type_ = RawThyra;
}

template<typename Scalar>
void MatrixFreeJacobianOperator<Scalar>::
setBaseEvaluationToNOXGroup(const Teuchos::RCP<const NOX::Abstract::Group>& base_group)
{
  TEUCHOS_ASSERT(setup_called_);
  TEUCHOS_ASSERT(nonnull(base_group));
  base_group_ = Teuchos::rcp_dynamic_cast<const NOX::Thyra::Group>(base_group,true);
  x_base_ = Teuchos::rcp_dynamic_cast<const NOX::Thyra::Vector>(base_group_->getXPtr(),true)->getThyraRCPVector();
  f_base_ = Teuchos::rcp_dynamic_cast<const NOX::Thyra::Vector>(base_group_->getFPtr(),true)->getThyraRCPVector();
  x_perturb_ = ::Thyra::createMember(x_base_->space());
  f_perturb_ = ::Thyra::createMember(f_base_->space());

  if (difference_type_ == Centered)
    f2_perturb_ = ::Thyra::createMember(f_base_->space());
    
  model_ = base_group_->getModel();
  base_evaluation_type_ = NoxGroup;
}

template<typename Scalar>
void MatrixFreeJacobianOperator<Scalar>::setBaseInArgs(const Teuchos::RCP< ::Thyra::ModelEvaluatorBase::InArgs<Scalar> >& base_in_args)
{
  in_args_ = base_in_args;
}

template<typename Scalar>
void MatrixFreeJacobianOperator<Scalar>::
setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& p)
{
  p->validateParametersAndSetDefaults(*this->getValidParameters());

  this->setMyParamList(p);

  difference_type_ = Teuchos::getIntegralValue<E_DifferenceType>(*p,"Difference Type");

  perturbation_type_ = Teuchos::getIntegralValue<E_PerturbationType>(*p,"Perturbation Algorithm");
  
  lambda_ = Teuchos::as<double>(p->get<double>("lambda"));
  
  user_defined_delta_ = Teuchos::as<double>(p->get<double>("User Defined delta Value"));

  setup_called_ = true;
}

template<typename Scalar>
Teuchos::RCP<const Teuchos::ParameterList> 
MatrixFreeJacobianOperator<Scalar>::getValidParameters() const
{
  if (is_null(valid_params_)) {
    valid_params_ = Teuchos::parameterList();

    Teuchos::setStringToIntegralParameter<E_DifferenceType>(
      "Difference Type",
      "Forward",
      "Difference method to be used to calculate the Jacobian-vector product.",
      Teuchos::tuple<std::string>("Forward","Central"),
      valid_params_.get()
      );

    Teuchos::setStringToIntegralParameter<E_PerturbationType>(
      "Perturbation Algorithm",
      "Salinger LOCA v1.0",
      "Determines the algorithm to use to calculation the perturbation delta.  ",
      Teuchos::tuple<std::string>("Salinger LOCA v1.0","KSP NOX 2001", "Knoll Keyes JCP 2004","User Defined"),
      valid_params_.get()
      );
    
    valid_params_->set<double>("lambda",1.0e-6);

    valid_params_->set<double>("User Defined delta Value",1.0e-6);

  }

  return valid_params_;
}

template<typename Scalar>
Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > 
MatrixFreeJacobianOperator<Scalar>::range() const
{
  return model_->get_f_space();
}

template<typename Scalar>
Teuchos::RCP<const ::Thyra::VectorSpaceBase< Scalar > > 
MatrixFreeJacobianOperator<Scalar>::domain() const
{
  return model_->get_x_space();
}

template<typename Scalar>
Teuchos::RCP<const ::Thyra::LinearOpBase< Scalar > > 
MatrixFreeJacobianOperator<Scalar>::clone() const
{
  TEUCHOS_ASSERT(true);
  return Teuchos::null;
}

template<typename Scalar>
bool 
MatrixFreeJacobianOperator<Scalar>::opSupportedImpl (::Thyra::EOpTransp M_trans) const
{
  if ( (M_trans == ::Thyra::TRANS) || (M_trans == ::Thyra::CONJTRANS) )
    return false;

  return true;
}

template<typename Scalar>
void 
MatrixFreeJacobianOperator<Scalar>::applyImpl(const ::Thyra::EOpTransp M_trans,
					      const ::Thyra::MultiVectorBase< Scalar > &thyra_mv_y,
					      const Teuchos::Ptr< ::Thyra::MultiVectorBase< Scalar > > &thyra_mv_u,
					      const Scalar alpha,
					      const Scalar beta) const
{
  TEUCHOS_ASSERT(setup_called_);

  // Use a directional derivative to approximate u = Jy
  
  /*
   * delta = scalar perturbation
   * x     = solution vector used to evaluate f
   * f     = function evaluation (RHS)
   * y     = vector that J is applied to
   *
   *            f(x + delta * y) - f(x)
   * u = Jy =   -----------------------
   *                     delta
   */

  const ::Thyra::VectorBase<Scalar>& y = dynamic_cast<const ::Thyra::VectorBase<Scalar>&>(thyra_mv_y);
  ::Thyra::VectorBase<Scalar>& u = dynamic_cast< ::Thyra::VectorBase<Scalar>&>(*thyra_mv_u);
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm_2_x = ::Thyra::norm(*x_base_);
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm_2_y = ::Thyra::norm(y);

  // Make sure the y-norm is not zero, otherwise we can get an inf perturbation from divide by zero
  if (norm_2_y == 0.0) {
    norm_2_y = 1.0;
    ::Thyra::assign(Teuchos::ptrFromRef(u),0.0);
    return;
  }
 
  if (perturbation_type_ == SalingerLOCA) {
    delta_ = lambda_ * (lambda_ + norm_2_x / norm_2_y);
  }
  else if (perturbation_type_ == KelleySalingerPawlowski) {
    Scalar inner_prod_x_y = ::Thyra::inner(*x_base_,y); 
    
    if (inner_prod_x_y==Teuchos::ScalarTraits<Scalar>::zero()) 
      inner_prod_x_y = 1.0e-12;
    
    delta_ = lambda_ * (1.0e-12 / lambda_ + std::fabs(inner_prod_x_y) / (norm_2_y * norm_2_y) ) * inner_prod_x_y / std::fabs(inner_prod_x_y);
  }
  else if (perturbation_type_ == KnollKeyes) {
    delta_ = lambda_ * ::Thyra::norm_1(*x_base_) / (Teuchos::as<double>(x_base_->space()->dim()) * norm_2_y) + lambda_;
  }
  else {
    delta_ = user_defined_delta_;
  }

  // perturbed solution: x_p = delta * y + x_0
  ::Thyra::V_StVpV(x_perturb_.ptr(),delta_,y,*x_base_);

  if (is_null(in_args_)) {
    in_args_ = Teuchos::rcp(new ::Thyra::ModelEvaluatorBase::InArgs<Scalar>(model_->createInArgs()));
    in_args_->setArgs(model_->getNominalValues());
  }

  in_args_->set_x(x_perturb_);
  
  if (is_null(out_args_))
    out_args_ = Teuchos::rcp(new ::Thyra::ModelEvaluatorBase::OutArgs<Scalar>(model_->createOutArgs()));

  out_args_->set_f(f_perturb_);

  // f_p = f(delta * y + x)
  model_->evalModel(*in_args_,*out_args_);

  // to be safe, remove arguments
  in_args_->set_x(Teuchos::null);
  out_args_->set_f(Teuchos::null);

  Scalar inv_delta = Teuchos::ScalarTraits<Scalar>::one() / delta_;
  
  if ( difference_type_ == Centered ) {
    // perturbed solution: x_p = -delta * y + x_0
    ::Thyra::V_StVpV(x_perturb_.ptr(),-delta_,y,*x_base_);

    out_args_->set_f(f2_perturb_);

    // f_p2 = f(-delta * y + x)
    model_->evalModel(*in_args_,*out_args_);

    // to be safe, remove arguments
    in_args_->set_x(Teuchos::null);
    out_args_->set_f(Teuchos::null);

    ::Thyra::V_StVpStV(ptrFromRef(u),inv_delta,*f_perturb_,-inv_delta,*f2_perturb_);
  }
  else {
    ::Thyra::V_StVpStV(ptrFromRef(u),inv_delta,*f_perturb_,-inv_delta,*f_base_);
  }

}

template<typename Scalar>
void MatrixFreeJacobianOperator<Scalar>::setLambda(double lambda)
{
  lambda_ = lambda;
}

template<typename Scalar>
Scalar MatrixFreeJacobianOperator<Scalar>::getLambda() const
{
  return lambda_;
}

template<typename Scalar>
void MatrixFreeJacobianOperator<Scalar>::setUserDefinedDelta(double delta)
{
  TEUCHOS_ASSERT(perturbation_type_ == UserDefined);
  user_defined_delta_ = delta;
}

template<typename Scalar>
Scalar MatrixFreeJacobianOperator<Scalar>::getUserDefinedDelta() const
{
  return user_defined_delta_;
}

template<typename Scalar>
Scalar MatrixFreeJacobianOperator<Scalar>::getDelta() const
{
  return delta_;
}
