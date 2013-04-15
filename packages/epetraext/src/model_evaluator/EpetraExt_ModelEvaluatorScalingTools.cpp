//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//@HEADER


#include "EpetraExt_ModelEvaluatorScalingTools.h"
#include "Teuchos_implicit_cast.hpp"
#include "Epetra_RowMatrix.h"

//
// Here in the implementation of scaling we write all scaling routines to
// specifically and individually handle every input and output object and we
// transfer all objects from one [In,Out]Arg container to another to make sure
// that that object has been specifically addressed.  We could just use the
// [In,Out]Args::setArgs(...) function to copy all arguments be default but
// then that would setup subtle errors where a new quality could be silently
// ignored and not be scaled correctly.  Therefore, we feel it is better to
// have to deal with *every* input and output object specifically with respect
// to scaling in order to avoid overlooking scaling.  In this way, if a new
// input or output object is added to [In,Out]Args but the code in this file
// is not updated, then that object will not be passed through and some type
// of error will be generated right away.  We feel this is the best behavior
// and it justifies having every scaling-related function take both an input
// and an output [In,Out]Args object and transferring the objects specifically.
//

namespace {


const std::string fwdScalingVecName = "fwdScalingVec";


// Assert that the input scaling vectors have been set up correctly
void assertModelVarScalings(
  const EpetraExt::ModelEvaluator::InArgs &varScalings
  )
{
  typedef EpetraExt::ModelEvaluator EME;
  TEUCHOS_TEST_FOR_EXCEPTION(
    (varScalings.supports(EME::IN_ARG_x) && varScalings.supports(EME::IN_ARG_x_dot))
    && (varScalings.get_x() != varScalings.get_x_dot()),
    std::logic_error,
    "Error, if scaling for x is given and x_dot is supported, then\n"
    "the scaling for x_dot must also be set and must be the same scaling\n"
    "as is used for x!"
    );
}
  


// Scale a single vector using a templated policy object to take care of what
// vector gets used.
template<class InArgsVectorGetterSetter>
void scaleModelVar(
  InArgsVectorGetterSetter vecGetterSetter, // Templated policy object!
  const EpetraExt::ModelEvaluator::InArgs &origVars,
  const EpetraExt::ModelEvaluator::InArgs &varScalings,
  EpetraExt::ModelEvaluator::InArgs *scaledVars,
  Teuchos::FancyOStream *out,
  Teuchos::EVerbosityLevel verbLevel
  )
{

  using Teuchos::null;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Ptr;
  using Teuchos::rcp_const_cast;


#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!scaledVars);
#endif

  RCP<const Epetra_Vector>
    orig_vec = vecGetterSetter.getVector(origVars);
  if ( !is_null(orig_vec) ) {
    RCP<const Epetra_Vector>
      inv_s_vec = vecGetterSetter.getVector(varScalings);
    if ( !is_null(inv_s_vec) ) {
      RCP<Epetra_Vector>
        scaled_vec = rcp_const_cast<Epetra_Vector>(
          vecGetterSetter.getVector(*scaledVars) );
      if ( is_null(scaled_vec) )
        scaled_vec = rcp(new Epetra_Vector(orig_vec->Map()));
      // See if there is a "hidden" forward scaling vector to use
      Ptr<const RCP<const Epetra_Vector> > fwd_s_vec =
        Teuchos::getOptionalEmbeddedObj<Epetra_Vector, RCP<const Epetra_Vector> >(
          inv_s_vec);
/*
        Teuchos::get_optional_extra_data<const RCP<const Epetra_Vector> >(
          inv_s_vec, fwdScalingVecName );
*/
      if ( !is_null(fwd_s_vec) ) {
        // Use the "hidden" forward scaling vector and multiply
        scaled_vec->Multiply( 1.0, **fwd_s_vec, *orig_vec, 0.0 );
      }
      else {
        // Just use the inverse scaling vector and divide
        EpetraExt::scaleModelVarsGivenInverseScaling(
          *orig_vec, *inv_s_vec, &*scaled_vec );
      }
      vecGetterSetter.setVector( scaled_vec, scaledVars );
    }
    else {
      vecGetterSetter.setVector( orig_vec, scaledVars );
    }
  }
  else {
    vecGetterSetter.setVector( null, scaledVars );
  }

}
  

// Scale variable bounds for a single vector using a templated policy object
// to take care of what vector gets used.
template<class InArgsVectorGetterSetter>
void scaleModelBound(
  InArgsVectorGetterSetter vecGetterSetter, // Templated policy object!
  const EpetraExt::ModelEvaluator::InArgs &origLowerBounds,
  const EpetraExt::ModelEvaluator::InArgs &origUpperBounds,
  const double infBnd,
  const EpetraExt::ModelEvaluator::InArgs &varScalings,
  EpetraExt::ModelEvaluator::InArgs *scaledLowerBounds,
  EpetraExt::ModelEvaluator::InArgs *scaledUpperBounds,
  Teuchos::FancyOStream *out,
  Teuchos::EVerbosityLevel verbLevel
  )
{

  using Teuchos::null;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!scaledLowerBounds);
  TEUCHOS_TEST_FOR_EXCEPT(!scaledUpperBounds);
#endif

  RCP<const Epetra_Vector>
    orig_lower_vec = vecGetterSetter.getVector(origLowerBounds);
  if ( !is_null(orig_lower_vec) ) {
    RCP<const Epetra_Vector>
      inv_s_vec = vecGetterSetter.getVector(varScalings);
    if ( !is_null(inv_s_vec) ) {
      TEUCHOS_TEST_FOR_EXCEPT("Can't handle scaling bounds yet!");
    }
    else {
      vecGetterSetter.setVector( orig_lower_vec, scaledLowerBounds );
    }
  }
  else {
    vecGetterSetter.setVector( null, scaledLowerBounds );
  }

  RCP<const Epetra_Vector>
    orig_upper_vec = vecGetterSetter.getVector(origUpperBounds);
  if ( !is_null(orig_upper_vec) ) {
    RCP<const Epetra_Vector>
      inv_s_vec = vecGetterSetter.getVector(varScalings);
    if ( !is_null(inv_s_vec) ) {
      TEUCHOS_TEST_FOR_EXCEPT("Can't handle scaling bounds yet!");
    }
    else {
      vecGetterSetter.setVector( orig_upper_vec, scaledUpperBounds );
    }
  }
  else {
    vecGetterSetter.setVector( null, scaledUpperBounds );
  }

}


// Unscale a single vector using a templated policy object to take care of
// what vector gets used.
template<class InArgsVectorGetterSetter>
void unscaleModelVar(
  InArgsVectorGetterSetter vecGetterSetter, // Templated policy object!
  const EpetraExt::ModelEvaluator::InArgs &scaledVars,
  const EpetraExt::ModelEvaluator::InArgs &varScalings,
  EpetraExt::ModelEvaluator::InArgs *origVars,
  Teuchos::FancyOStream *out,
  Teuchos::EVerbosityLevel verbLevel
  )
{

  using Teuchos::null;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;
  using Teuchos::includesVerbLevel;


#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!origVars);
#endif

  RCP<const Epetra_Vector>
    scaled_vec = vecGetterSetter.getVector(scaledVars);
  if ( !is_null(scaled_vec) ) {
    RCP<const Epetra_Vector>
      inv_s_vec = vecGetterSetter.getVector(varScalings);
    if ( !is_null(inv_s_vec) ) {
      RCP<Epetra_Vector>
        orig_vec = rcp_const_cast<Epetra_Vector>(
          vecGetterSetter.getVector(*origVars) );
      if ( is_null(orig_vec) )
        orig_vec = rcp(new Epetra_Vector(scaled_vec->Map()));
      EpetraExt::unscaleModelVarsGivenInverseScaling(
        *scaled_vec, *inv_s_vec,  &*orig_vec );
      if (out && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH)) {
        *out << "\nUnscaled vector "<<vecGetterSetter.getName()<<":\n";
        Teuchos::OSTab tab(*out);
        orig_vec->Print(*out);
      }
      vecGetterSetter.setVector( orig_vec, origVars );
    }
    else {
      vecGetterSetter.setVector( scaled_vec, origVars );
    }
  }
  else {
    vecGetterSetter.setVector( null, origVars );
  }

}


// Scale a single vector using a templated policy object to take care of what
// vector gets used.
template<class OutArgsVectorGetterSetter>
void scaleModelFunc(
  OutArgsVectorGetterSetter vecGetterSetter, // Templated policy object!
  const EpetraExt::ModelEvaluator::OutArgs &origFuncs,
  const EpetraExt::ModelEvaluator::OutArgs &funcScalings,
  EpetraExt::ModelEvaluator::OutArgs *scaledFuncs,
  Teuchos::FancyOStream *out,
  Teuchos::EVerbosityLevel verbLevel
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(0==scaledFuncs);
  Teuchos::RCP<Epetra_Vector>
    func = vecGetterSetter.getVector(origFuncs);
  if (!is_null(func) ) {
    Teuchos::RCP<const Epetra_Vector>
      funcScaling = vecGetterSetter.getVector(funcScalings);
    if (!is_null(funcScaling) ) {
      EpetraExt::scaleModelFuncGivenForwardScaling( *funcScaling, &*func );
    }
  }
  vecGetterSetter.setVector( func, scaledFuncs );
}


} // namespace


void EpetraExt::gatherModelNominalValues(
  const ModelEvaluator &model,
  ModelEvaluator::InArgs *nominalValues
  )
{

  using Teuchos::implicit_cast;
  typedef ModelEvaluator EME;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!nominalValues);
#endif

  *nominalValues = model.createInArgs();

  if(nominalValues->supports(EME::IN_ARG_x)) {
    nominalValues->set_x(model.get_x_init());
  }

  if(nominalValues->supports(EME::IN_ARG_x_dot)) {
    nominalValues->set_x_dot(model.get_x_dot_init());
  }
  
  for( int l = 0; l < nominalValues->Np(); ++l ) {
    nominalValues->set_p( l, model.get_p_init(l) );
  }
  
  if(nominalValues->supports(EME::IN_ARG_t)) {
    nominalValues->set_t(model.get_t_init());
  }

}


void EpetraExt::gatherModelBounds(
  const ModelEvaluator &model,
  ModelEvaluator::InArgs *lowerBounds,
  ModelEvaluator::InArgs *upperBounds
  )
{

  using Teuchos::implicit_cast;
  typedef ModelEvaluator EME;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!lowerBounds);
  TEUCHOS_TEST_FOR_EXCEPT(!upperBounds);
#endif

  *lowerBounds = model.createInArgs();
  *upperBounds = model.createInArgs();

  if(lowerBounds->supports(EME::IN_ARG_x)) {
    lowerBounds->set_x(model.get_x_lower_bounds());
    upperBounds->set_x(model.get_x_upper_bounds());
  }
  
  for( int l = 0; l < lowerBounds->Np(); ++l ) {
    lowerBounds->set_p( l, model.get_p_lower_bounds(l) );
    upperBounds->set_p( l, model.get_p_upper_bounds(l) );
  }
  
  if(lowerBounds->supports(EME::IN_ARG_t)) {
    lowerBounds->set_t(model.get_t_lower_bound());
    upperBounds->set_t(model.get_t_upper_bound());
  }

}


void EpetraExt::scaleModelVars(
  const ModelEvaluator::InArgs &origVars,
  const ModelEvaluator::InArgs &varScalings,
  ModelEvaluator::InArgs *scaledVars,
  Teuchos::FancyOStream *out,
  Teuchos::EVerbosityLevel verbLevel
  )
{
  typedef ModelEvaluator EME;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!scaledVars);
  assertModelVarScalings(varScalings);
#endif

  if (origVars.supports(EME::IN_ARG_x)) {
    scaleModelVar( InArgsGetterSetter_x(), origVars, varScalings, scaledVars,
      out, verbLevel );
  }

  if (origVars.supports(EME::IN_ARG_x_dot)) {
    scaleModelVar( InArgsGetterSetter_x_dot(), origVars, varScalings, scaledVars,
      out, verbLevel );
  }

  const int Np = origVars.Np();
  for ( int l = 0; l < Np; ++l ) {
    scaleModelVar( InArgsGetterSetter_p(l), origVars, varScalings, scaledVars,
      out, verbLevel );
  }

  if (origVars.supports(EME::IN_ARG_x_poly)) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      !is_null(varScalings.get_x()), std::logic_error,
      "Error, can't hanlde scaling of x_poly yet!"
      );
    scaledVars->set_x_poly(origVars.get_x_poly());
  }

  if (origVars.supports(EME::IN_ARG_x_dot_poly)) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      !is_null(varScalings.get_x()), std::logic_error,
      "Error, can't hanlde scaling of x_dot_poly yet!"
      );
    scaledVars->set_x_dot_poly(origVars.get_x_dot_poly());
  }

  if (origVars.supports(EME::IN_ARG_t)) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      varScalings.get_t() > 0.0, std::logic_error,
      "Error, can't hanlde scaling of t yet!"
      );
    scaledVars->set_t(origVars.get_t());
  }

  if (origVars.supports(EME::IN_ARG_alpha)) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      varScalings.get_alpha() > 0.0, std::logic_error,
      "Error, can't hanlde scaling of alpha yet!"
      );
    scaledVars->set_alpha(origVars.get_alpha());
  }

  if (origVars.supports(EME::IN_ARG_beta)) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      varScalings.get_beta() > 0.0, std::logic_error,
      "Error, can't hanlde scaling of beta yet!"
      );
    scaledVars->set_beta(origVars.get_beta());
  }

  // ToDo: Support other input arguments?

}


void EpetraExt::scaleModelBounds(
  const ModelEvaluator::InArgs &origLowerBounds,
  const ModelEvaluator::InArgs &origUpperBounds,
  const double infBnd,
  const ModelEvaluator::InArgs &varScalings,
  ModelEvaluator::InArgs *scaledLowerBounds,
  ModelEvaluator::InArgs *scaledUpperBounds,
  Teuchos::FancyOStream *out,
  Teuchos::EVerbosityLevel verbLevel
  )
{

  typedef ModelEvaluator EME;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!scaledLowerBounds);
  TEUCHOS_TEST_FOR_EXCEPT(!scaledUpperBounds);
  assertModelVarScalings(varScalings);
#endif

  if (origLowerBounds.supports(EME::IN_ARG_x)) {
    scaleModelBound(
      InArgsGetterSetter_x(), origLowerBounds, origUpperBounds, infBnd,
      varScalings, scaledLowerBounds, scaledUpperBounds,
      out, verbLevel );
  }

  if (origLowerBounds.supports(EME::IN_ARG_x_dot)) {
    scaleModelBound(
      InArgsGetterSetter_x_dot(), origLowerBounds, origUpperBounds, infBnd,
      varScalings, scaledLowerBounds, scaledUpperBounds,
      out, verbLevel );
  }

  const int Np = origLowerBounds.Np();
  for ( int l = 0; l < Np; ++l ) {
    scaleModelBound(
      InArgsGetterSetter_p(l), origLowerBounds, origUpperBounds, infBnd,
      varScalings, scaledLowerBounds, scaledUpperBounds,
      out, verbLevel );
  }

  // ToDo: Support other input arguments?

}


void EpetraExt::unscaleModelVars(
  const ModelEvaluator::InArgs &scaledVars,
  const ModelEvaluator::InArgs &varScalings,
  ModelEvaluator::InArgs *origVars,
  Teuchos::FancyOStream *out,
  Teuchos::EVerbosityLevel verbLevel
  )
{

  using Teuchos::RCP;
  using Teuchos::includesVerbLevel;
  typedef ModelEvaluator EME;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!origVars);
  assertModelVarScalings(varScalings);
#endif

  // Print scaling vectors

  if (out && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH)) {
    RCP<const Epetra_Vector> inv_s_x;
    if ( scaledVars.supports(EME::IN_ARG_x) &&
      !is_null(inv_s_x=varScalings.get_x()) )
    {
      *out << "\nState inverse scaling vector inv_s_x:\n";
      Teuchos::OSTab tab(*out);
      inv_s_x->Print(*out);
    }
  }

  // Scal the input varaibles
  
  if (scaledVars.supports(EME::IN_ARG_x_dot)) {
    unscaleModelVar( InArgsGetterSetter_x_dot(), scaledVars, varScalings, origVars,
      out, verbLevel );
  }

  if (scaledVars.supports(EME::IN_ARG_x)) {
    unscaleModelVar( InArgsGetterSetter_x(), scaledVars, varScalings, origVars,
      out, verbLevel );
  }

  const int Np = scaledVars.Np();
  for ( int l = 0; l < Np; ++l ) {
    unscaleModelVar( InArgsGetterSetter_p(l), scaledVars, varScalings, origVars,
      out, verbLevel );
  }

  if (scaledVars.supports(EME::IN_ARG_x_dot_poly)) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      !is_null(varScalings.get_x()), std::logic_error,
      "Error, can't hanlde unscaling of x_dot_poly yet!"
      );
    origVars->set_x_dot_poly(scaledVars.get_x_dot_poly());
  }

  if (scaledVars.supports(EME::IN_ARG_x_poly)) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      !is_null(varScalings.get_x()), std::logic_error,
      "Error, can't hanlde unscaling of x_poly yet!"
      );
    origVars->set_x_poly(scaledVars.get_x_poly());
  }

  if (scaledVars.supports(EME::IN_ARG_t)) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      varScalings.get_t() > 0.0, std::logic_error,
      "Error, can't hanlde unscaling of t yet!"
      );
    origVars->set_t(scaledVars.get_t());
  }

  if (scaledVars.supports(EME::IN_ARG_alpha)) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      varScalings.get_alpha() > 0.0, std::logic_error,
      "Error, can't hanlde unscaling of alpha yet!"
      );
    origVars->set_alpha(scaledVars.get_alpha());
  }

  if (scaledVars.supports(EME::IN_ARG_beta)) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      varScalings.get_beta() > 0.0, std::logic_error,
      "Error, can't hanlde unscaling of beta yet!"
      );
    origVars->set_beta(scaledVars.get_beta());
  }

}


void EpetraExt::scaleModelFuncs(
  const ModelEvaluator::OutArgs &origFuncs,
  const ModelEvaluator::InArgs &varScalings,
  const ModelEvaluator::OutArgs &funcScalings,
  ModelEvaluator::OutArgs *scaledFuncs,
  bool *allFuncsWhereScaled,
  Teuchos::FancyOStream *out,
  Teuchos::EVerbosityLevel verbLevel
  )
{

  using Teuchos::RCP;
  typedef ModelEvaluator EME;

  TEUCHOS_TEST_FOR_EXCEPT(0==allFuncsWhereScaled);

  *allFuncsWhereScaled = true;

  const int Np = origFuncs.Np();
  const int Ng = origFuncs.Ng();

  // f
  if ( origFuncs.supports(EME::OUT_ARG_f) && !is_null(origFuncs.get_f()) ) {
    scaleModelFunc( OutArgsGetterSetter_f(), origFuncs, funcScalings, scaledFuncs,
      out, verbLevel );
  }

  // f_poly
  if (
    origFuncs.supports(EME::OUT_ARG_f_poly)
    && !is_null(origFuncs.get_f_poly())
    )
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      !is_null(funcScalings.get_f()), std::logic_error,
      "Error, we can't handle scaling of f_poly yet!"
      );
    scaledFuncs->set_f_poly(origFuncs.get_f_poly());
  }

  // g(j)
  for ( int j = 0; j < Ng; ++j ) {
    scaleModelFunc( OutArgsGetterSetter_g(j), origFuncs, funcScalings, scaledFuncs,
      out, verbLevel );
  }

  // W
  RCP<Epetra_Operator> W;
  if ( origFuncs.supports(EME::OUT_ARG_W) && !is_null(W=origFuncs.get_W()) ) {
    bool didScaling = false;
    scaleModelFuncFirstDerivOp(
      varScalings.get_x().get(), funcScalings.get_f().get(),
      &*W, &didScaling
      );
    if(didScaling)
      scaledFuncs->set_W(W);
    else
      *allFuncsWhereScaled = false;
  }

  // DfDp(l)
  for ( int l = 0; l < Np; ++l ) {
    EME::Derivative orig_DfDp_l;
    if (
      !origFuncs.supports(EME::OUT_ARG_DfDp,l).none()
      && !(orig_DfDp_l=origFuncs.get_DfDp(l)).isEmpty()
      )
    {
      EME::Derivative scaled_DfDp_l;
      bool didScaling = false;
      scaleModelFuncFirstDeriv(
        orig_DfDp_l, varScalings.get_p(l).get(), funcScalings.get_f().get(),
        &scaled_DfDp_l, &didScaling
        );
      if(didScaling)
        scaledFuncs->set_DfDp(l,scaled_DfDp_l);
      else
        *allFuncsWhereScaled = false;
    }

  }

  // DgDx_dot(j), DgDx(j), and DgDp(j,l)
  for ( int j = 0; j < Ng; ++j ) {

    EME::Derivative orig_DgDx_dot_j;
    if (
      !origFuncs.supports(EME::OUT_ARG_DgDx_dot,j).none()
      && !(orig_DgDx_dot_j=origFuncs.get_DgDx_dot(j)).isEmpty()
      )
    {
      EME::Derivative scaled_DgDx_dot_j;
      bool didScaling = false;
      scaleModelFuncFirstDeriv(
        orig_DgDx_dot_j, varScalings.get_x().get(), funcScalings.get_g(j).get(),
        &scaled_DgDx_dot_j, &didScaling
        );
      if(didScaling)
        scaledFuncs->set_DgDx_dot(j,scaled_DgDx_dot_j);
      else
        *allFuncsWhereScaled = false;
    }

    EME::Derivative orig_DgDx_j;
    if (
      !origFuncs.supports(EME::OUT_ARG_DgDx,j).none()
      && !(orig_DgDx_j=origFuncs.get_DgDx(j)).isEmpty()
      )
    {
      EME::Derivative scaled_DgDx_j;
      bool didScaling = false;
      scaleModelFuncFirstDeriv(
        orig_DgDx_j, varScalings.get_x().get(), funcScalings.get_g(j).get(),
        &scaled_DgDx_j, &didScaling
        );
      if(didScaling)
        scaledFuncs->set_DgDx(j,scaled_DgDx_j);
      else
        *allFuncsWhereScaled = false;
    }

    for ( int l = 0; l < Np; ++l ) {
      EME::Derivative orig_DgDp_j_l;
      if (
        !origFuncs.supports(EME::OUT_ARG_DgDp,j,l).none()
        && !(orig_DgDp_j_l=origFuncs.get_DgDp(j,l)).isEmpty()
        )
      {
        EME::Derivative scaled_DgDp_j_l;
        bool didScaling = false;
        scaleModelFuncFirstDeriv(
          orig_DgDp_j_l, varScalings.get_p(l).get(), funcScalings.get_g(j).get(),
          &scaled_DgDp_j_l, &didScaling
          );
        if(didScaling)
          scaledFuncs->set_DgDp(j,l,scaled_DgDp_j_l);
        else
          *allFuncsWhereScaled = false;
      }
    }
  }

}


Teuchos::RCP<const Epetra_Vector>
EpetraExt::createInverseModelScalingVector(
  Teuchos::RCP<const Epetra_Vector> const& scalingVector
  )
{
  Teuchos::RCP<Epetra_Vector> invScalingVector =
    Teuchos::rcpWithEmbeddedObj(
      new Epetra_Vector(scalingVector->Map()),
      scalingVector
      );
  invScalingVector->Reciprocal(*scalingVector);
  return invScalingVector;
  // Above, we embedd the forward scaling vector.  This is done in order to
  // achieve the exact same numerics as before this refactoring and to improve
  // runtime speed and accruacy.
}


void EpetraExt::scaleModelVarsGivenInverseScaling(
  const Epetra_Vector &origVars,
  const Epetra_Vector &invVarScaling,
  Epetra_Vector *scaledVars
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!scaledVars);
  TEUCHOS_TEST_FOR_EXCEPT(!origVars.Map().SameAs(invVarScaling.Map()));
  TEUCHOS_TEST_FOR_EXCEPT(!origVars.Map().SameAs(scaledVars->Map()));
#endif
  const int localDim = origVars.Map().NumMyElements();
  for ( int i = 0; i < localDim; ++i )
    (*scaledVars)[i] = origVars[i] / invVarScaling[i];
}


void EpetraExt::scaleModelVarBoundsGivenInverseScaling(
  const Epetra_Vector &origLowerBounds,
  const Epetra_Vector &origUpperBounds,
  const double infBnd,
  const Epetra_Vector &invVarScaling,
  Epetra_Vector *scaledLowerBounds,
  Epetra_Vector *scaledUpperBounds
  )
{
  TEUCHOS_TEST_FOR_EXCEPT("ToDo: Implement!");
}


void EpetraExt::unscaleModelVarsGivenInverseScaling(
  const Epetra_Vector &origVars,
  const Epetra_Vector &invVarScaling,
  Epetra_Vector *scaledVars
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(0==scaledVars);
  scaledVars->Multiply( 1.0, invVarScaling, origVars, 0.0 );
}


void EpetraExt::scaleModelFuncGivenForwardScaling(
  const Epetra_Vector &fwdFuncScaling,
  Epetra_Vector *funcs
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(0==funcs);
  funcs->Multiply( 1.0,  fwdFuncScaling, *funcs, 0.0 );
  // Note: Above is what Epetra_LinearProblem does to scale the RHS and LHS
  // vectors so this type of argument aliasing must be okay in Epetra!
}


void EpetraExt::scaleModelFuncFirstDerivOp(
  const Epetra_Vector *invVarScaling,
  const Epetra_Vector *fwdFuncScaling,
  Epetra_Operator *funcDerivOp,
  bool *didScaling
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(0==funcDerivOp);
  TEUCHOS_TEST_FOR_EXCEPT(0==didScaling);
  *didScaling = false; // Assume not scaled to start
  Epetra_RowMatrix *funcDerivRowMatrix
    = dynamic_cast<Epetra_RowMatrix*>(funcDerivOp);
  if (funcDerivRowMatrix) {
    if (fwdFuncScaling)
      funcDerivRowMatrix->LeftScale(*fwdFuncScaling);
    if (invVarScaling)
      funcDerivRowMatrix->RightScale(*invVarScaling);
    *didScaling = true;
    // Note that above I do the func scaling before the var scaling since that
    // is the same order it was done for W in Thyra::EpetraModelEvaluator
  }
}


void EpetraExt::scaleModelFuncFirstDeriv(
  const ModelEvaluator::Derivative &origFuncDeriv,
  const Epetra_Vector *invVarScaling,
  const Epetra_Vector *fwdFuncScaling,
  ModelEvaluator::Derivative *scaledFuncDeriv,
  bool *didScaling
  )
{
  using Teuchos::RCP;
  typedef ModelEvaluator EME;
  TEUCHOS_TEST_FOR_EXCEPT(0==scaledFuncDeriv);
  TEUCHOS_TEST_FOR_EXCEPT(0==didScaling);
  *didScaling = false;
  const RCP<Epetra_MultiVector>
    funcDerivMv = origFuncDeriv.getMultiVector();
  const EME::EDerivativeMultiVectorOrientation
    funcDerivMv_orientation = origFuncDeriv.getMultiVectorOrientation();
  if(!is_null(funcDerivMv)) {
    if ( funcDerivMv_orientation == EME::DERIV_MV_BY_COL )
    {
      if (fwdFuncScaling) {
        funcDerivMv->Multiply(1.0, *fwdFuncScaling, *funcDerivMv, 0.0);
      }
      if (invVarScaling) {
        TEUCHOS_TEST_FOR_EXCEPT("ToDo: Scale rows!");
        //funcDerivMv->Multiply(1.0, *funcDerivMv, *invVarScaling, 0.0);
      }
    }
    else if ( funcDerivMv_orientation == EME::DERIV_TRANS_MV_BY_ROW )
    {
      if (invVarScaling) {
        funcDerivMv->Multiply(1.0, *invVarScaling, *funcDerivMv, 0.0);
      }
      if (fwdFuncScaling) {
        TEUCHOS_TEST_FOR_EXCEPT("ToDo: Scale rows!");
        //funcDerivMv->Multiply(1.0, *funcDerivMv, *fwdFuncScaling, 0.0);
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT("Should not get here!");
    }
    *scaledFuncDeriv = EME::Derivative(funcDerivMv,funcDerivMv_orientation);
    *didScaling = true;
  }
  else {
    RCP<Epetra_Operator>
      funcDerivOp = origFuncDeriv.getLinearOp();
    TEUCHOS_TEST_FOR_EXCEPT(is_null(funcDerivOp));
    scaleModelFuncFirstDerivOp( invVarScaling, fwdFuncScaling,
      &*funcDerivOp, didScaling );
    if (didScaling)
      *scaledFuncDeriv = EME::Derivative(funcDerivOp);
  }
}
