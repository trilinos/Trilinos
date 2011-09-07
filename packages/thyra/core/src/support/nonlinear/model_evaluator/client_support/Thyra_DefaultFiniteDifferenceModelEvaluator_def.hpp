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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DEF_HPP
#define THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DEF_HPP

#include "Thyra_DefaultFiniteDifferenceModelEvaluator_decl.hpp"


namespace Thyra {


// Constructors/initializers/accessors/utilities


template<class Scalar>
DefaultFiniteDifferenceModelEvaluator<Scalar>::DefaultFiniteDifferenceModelEvaluator()
{}


template<class Scalar>
void DefaultFiniteDifferenceModelEvaluator<Scalar>::initialize(
  const RCP<ModelEvaluator<Scalar> > &thyraModel,
  const RCP<DirectionalFiniteDiffCalculator<Scalar> > &direcFiniteDiffCalculator_in
  )
{
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
  direcFiniteDiffCalculator_ = direcFiniteDiffCalculator_in;
}


// Public functions overridden from Teuchos::Describable


template<class Scalar>
std::string DefaultFiniteDifferenceModelEvaluator<Scalar>::description() const
{
  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::DefaultFiniteDifferenceModelEvaluator{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DefaultFiniteDifferenceModelEvaluator<Scalar>::createOutArgsImpl() const
{
  typedef ModelEvaluatorBase MEB;
  typedef MEB::DerivativeSupport DS;
  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  const MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  const int l_Np = wrappedOutArgs.Np(), l_Ng = wrappedOutArgs.Ng();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(l_Np,l_Ng);
  outArgs.setSupports(wrappedOutArgs);
  if (wrappedOutArgs.supports(MEB::OUT_ARG_f)) {
    for( int l = 0; l < l_Np; ++l ) {
      outArgs.setSupports(MEB::OUT_ARG_DfDp,l,MEB::DERIV_MV_BY_COL);
    }
  }
  for( int j = 0; j < l_Ng; ++j ) {
    for( int l = 0; l < l_Np; ++l ) {
      outArgs.setSupports( MEB::OUT_ARG_DgDp , j, l, MEB::DERIV_MV_BY_COL);
    }
  }
  // ToDo: Add support for more derivatives as needed!
  return outArgs;
}


template<class Scalar>
void DefaultFiniteDifferenceModelEvaluator<Scalar>::evalModelImpl(
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  //typedef typename ST::magnitudeType ScalarMag;
  typedef ModelEvaluatorBase MEB;
  namespace DFDCT = DirectionalFiniteDiffCalculatorTypes;

  typedef RCP<VectorBase<Scalar> > V_ptr;
  typedef RCP<const VectorBase<Scalar> > CV_ptr;
  typedef RCP<MultiVectorBase<Scalar> > MV_ptr;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Thyra::DefaultFiniteDifferenceModelEvaluator",inArgs,outArgs
    );

  //
  // Note: Just do derivatives DfDp(l) and DgDp(j,l) for now!
  //

  const RCP<const VectorSpaceBase<Scalar> >
    p_space = thyraModel->get_p_space(0),
    g_space = thyraModel->get_g_space(0);

  //
  // A) Compute the base point
  //

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
    *out << "\nComputing the base point ...\n";

  const int l_Np = outArgs.Np();
  const int l_Ng = outArgs.Ng();
  MEB::InArgs<Scalar> wrappedInArgs = inArgs;
  MEB::OutArgs<Scalar> baseFunc = thyraModel->createOutArgs();
  if( outArgs.supports(MEB::OUT_ARG_f) && outArgs.get_f().get() )
    baseFunc.set_f(outArgs.get_f());
  for( int j = 0; j < l_Ng; ++j ) {
    V_ptr g_j;
    if( (g_j=outArgs.get_g(j)).get() )
      baseFunc.set_g(j,g_j);
  }
  // 2007/08/27: We really should really try to allow some derivatives to pass
  // through and some derivatives to be computed by finite differences. Right
  // now, if you use this class, all derivatives w.r.t. parameters are finite
  // differenced and that is not given the user enough control!

  thyraModel->evalModel(wrappedInArgs,baseFunc);

  bool failed = baseFunc.isFailed();

  //
  // B) Compute the derivatives
  //
 
  if(!failed) {

    // a) Determine what derivatives you need to support first

    Array<int> compute_DfDp;
    Array<Array<int> > compute_DgDp(l_Ng);
    DFDCT::SelectedDerivatives selectedDerivs;

    for ( int l = 0; l < l_Np; ++l ) {

      // DfDp(l)
      if(
        outArgs.supports(MEB::OUT_ARG_DfDp,l).none()==false
        &&
        outArgs.get_DfDp(l).isEmpty()==false
        )
      {
        selectedDerivs.supports(MEB::OUT_ARG_DfDp,l);
        compute_DfDp.push_back(true);
      }
      else
      {
        compute_DfDp.push_back(false);
      }

      // DgDp(j=0...,l)
      for ( int j = 0; j < l_Ng; ++j ) {
        if(
          outArgs.supports(MEB::OUT_ARG_DgDp,j,l).none()==false
          &&
          outArgs.get_DgDp(j,l).isEmpty()==false
          )
        {
          selectedDerivs.supports(MEB::OUT_ARG_DgDp,j,l);
          compute_DgDp[j].push_back(true);
        }
        else
        {
          compute_DgDp[j].push_back(false);
        }
      }
    }

    // b) Create the deriv OutArgs and set the output objects that need to be
    // computed with finite differences
 
    MEB::OutArgs<Scalar>
      deriv = direcFiniteDiffCalculator_->createOutArgs(
        *thyraModel, selectedDerivs );

    for ( int l = 0; l < l_Np; ++l ) {
      if ( compute_DfDp[l] )
        deriv.set_DfDp(l,outArgs.get_DfDp(l));
      for ( int j = 0; j < l_Ng; ++j ) {
        if ( compute_DgDp[j][l] )
          deriv.set_DgDp(j,l,outArgs.get_DgDp(j,l));
      }
    }

    // c) Compute the missing functions with finite differences!

    direcFiniteDiffCalculator_->calcDerivatives(
      *thyraModel,inArgs,baseFunc,deriv
      );

  }

  if(failed) {
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW))
      *out << "\nEvaluation failed, returning NaNs ...\n";
    outArgs.setFailed();
  }

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
 
}


} // namespace Thyra


#endif // THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DEF_HPP
