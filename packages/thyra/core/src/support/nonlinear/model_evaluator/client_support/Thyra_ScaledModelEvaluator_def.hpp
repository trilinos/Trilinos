// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DEF_HPP
#define THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DEF_HPP

#include "Thyra_DefaultFiniteDifferenceModelEvaluator_decl.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Thyra {


// Constructors/initializers/accessors/utilities


template<class Scalar>
ScaledModelEvaluator<Scalar>::ScaledModelEvaluator()
{}


// Public functions overridden from Teuchos::Describable


template<class Scalar>
std::string ScaledModelEvaluator<Scalar>::description() const
{
  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::ScaledModelEvaluator{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}


template<class Scalar>
void ScaledModelEvaluator<Scalar>::
set_f_scaling(const RCP<const Thyra::VectorBase<Scalar> >& f_scaling)
{
  f_scaling_ = f_scaling;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
void ScaledModelEvaluator<Scalar>::evalModelImpl(
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Thyra::ScaledModelEvaluator",inArgs,outArgs
    );

  thyraModel->evalModel(inArgs, outArgs);

  if (nonnull(f_scaling_)) {

    const RCP<VectorBase<Scalar> > f = outArgs.get_f();
    if (nonnull(f)) {
      ele_wise_scale(*f_scaling_, f.ptr());
    }
    
    const RCP<LinearOpBase<Scalar> > W_op = outArgs.get_W_op();
    if (nonnull(W_op)) {
      const RCP<ScaledLinearOpBase<Scalar> > W_scaled =
        rcp_dynamic_cast<ScaledLinearOpBase<Scalar> >(W_op, true);
      W_scaled->scaleLeft(*f_scaling_);
    }

  }

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
 
}


} // namespace Thyra


#endif // THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DEF_HPP
