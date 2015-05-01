// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_EXPLICIT_MODEL_EVALUATOR_IMPL_HPP
#define PANZER_EXPLICIT_MODEL_EVALUATOR_IMPL_HPP

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

#include "Thyra_get_Epetra_Operator.hpp"
#include "EpetraExt_RowMatrixOut.h"

namespace panzer {

template<typename Scalar>
ExplicitModelEvaluator<Scalar>::
ExplicitModelEvaluator(const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > & model,
                       bool constantMassMatrix,
                       bool useLumpedMass)
   : Thyra::ModelEvaluatorDelegatorBase<Scalar>(model)
   , constantMassMatrix_(constantMassMatrix)
   , massLumping_(useLumpedMass)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // extract a panzer::ModelEvaluator if appropriate
  panzerModel_ = rcp_dynamic_cast<panzer::ModelEvaluator<Scalar> >(model);

  // extract a panzer::ModelEvaluator_Epetra if appropriate
  RCP<Thyra::EpetraModelEvaluator> epME = rcp_dynamic_cast<Thyra::EpetraModelEvaluator>(model);
  if(epME!=Teuchos::null)
    panzerEpetraModel_ = rcp_dynamic_cast<const panzer::ModelEvaluator_Epetra>(epME->getEpetraModel());

  // note at this point its possible that panzerModel_ = panzerEpetraModel_ = Teuchos::null

  buildArgsPrototypes();
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> ExplicitModelEvaluator<Scalar>::
getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar> nomVals = createInArgs();
  nomVals.setArgs(this->getUnderlyingModel()->getNominalValues(),true);

  return nomVals;
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> ExplicitModelEvaluator<Scalar>::
createInArgs() const
{
  return prototypeInArgs_;
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar> ExplicitModelEvaluator<Scalar>::
createOutArgs() const
{
  return prototypeOutArgs_;
}

template<typename Scalar>
Teuchos::RCP<panzer::ModelEvaluator<Scalar> > ExplicitModelEvaluator<Scalar>::
getPanzerUnderlyingModel()
{
  return Teuchos::rcp_dynamic_cast<panzer::ModelEvaluator<Scalar> >(this->getNonconstUnderlyingModel());
}

template<typename Scalar>
void ExplicitModelEvaluator<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP; 

  RCP<const Thyra::ModelEvaluator<Scalar> > under_me = this->getUnderlyingModel();

  // intialize a zero to get rid of the x-dot 
  if(zero_==Teuchos::null) {
    zero_ = Thyra::createMember(*under_me->get_x_space());
    Thyra::assign(zero_.ptr(),0.0);
  }

  MEB::InArgs<Scalar> under_inArgs = under_me->createInArgs();
  under_inArgs.setArgs(inArgs);
  under_inArgs.set_x_dot(zero_); // no time derivative...this may be a problem
                                 // for VMS stabilization in particular
                                 // or more importantly cause performance to suck!

  Teuchos::RCP<Thyra::VectorBase<Scalar> > f = outArgs.get_f();

  MEB::OutArgs<Scalar> under_outArgs = under_me->createOutArgs();
  under_outArgs.setArgs(outArgs);
  if(f!=Teuchos::null) {
    // build a scrap vector that will contain the weak residual
    if(scrap_f_==Teuchos::null)
      scrap_f_ = Thyra::createMember(*under_me->get_f_space());
    
    Thyra::assign(scrap_f_.ptr(),0.0);
    under_outArgs.set_f(scrap_f_);
  }

  under_me->evalModel(under_inArgs,under_outArgs);

  if(f!=Teuchos::null) {
    if(invMassMatrix_==Teuchos::null || constantMassMatrix_==false)
      buildInverseMassMatrix();

    // invert the mass matrix
    Thyra::apply(*invMassMatrix_,Thyra::NOTRANS,*scrap_f_,f.ptr()); 
  }
}

template<typename Scalar>
void ExplicitModelEvaluator<Scalar>::
buildInverseMassMatrix() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;
  using Thyra::createMember;
  
  RCP<const Thyra::ModelEvaluator<Scalar> > me = this->getUnderlyingModel();

  // first allocate space for the mass matrix
  RCP<Thyra::LinearOpBase<Scalar> > mass = me->create_W_op();

  // intialize a zero to get rid of the x-dot 
  if(zero_==Teuchos::null) {
    zero_ = Thyra::createMember(*me->get_x_space());
    Thyra::assign(zero_.ptr(),0.0);
  }
  
  // request only the mass matrix from the physics
  // Model evaluator builds: alpha*u_dot + beta*F(u) = 0
  MEB::InArgs<Scalar>  inArgs  = me->createInArgs();
  inArgs.set_x(zero_);
  inArgs.set_x_dot(zero_);
  inArgs.set_alpha(-1.0);
  inArgs.set_beta(0.0);

  // set the one time beta to ensure dirichlet conditions
  // are correctly included in the mass matrix: do it for
  // both epetra and Tpetra. 
  if(panzerModel_!=Teuchos::null)
    panzerModel_->setOneTimeDirichletBeta(-1.0);
  else if(panzerEpetraModel_!=Teuchos::null)
    panzerEpetraModel_->setOneTimeDirichletBeta(-1.0);
  else {
    // assuming the underlying model is a delegator, walk through
    // the decerator hierarchy until you find a panzer::ME or panzer::EpetraME.
    // If you don't find one, then throw because you are in a load of trouble anyway!
    setOneTimeDirichletBeta(-1.0,*this->getUnderlyingModel());
  }

  // set only the mass matrix
  MEB::OutArgs<Scalar> outArgs = me->createOutArgs();
  outArgs.set_W_op(mass);

  // this will fill the mass matrix operator 
  me->evalModel(inArgs,outArgs);

  // Teuchos::RCP<const Epetra_CrsMatrix> crsMat = Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*mass));
  // EpetraExt::RowMatrixToMatrixMarketFile("expmat.mm",*crsMat);

  if(!massLumping_) {
    invMassMatrix_ = Thyra::inverse<Scalar>(*me->get_W_factory(),mass);
  }
  else {
    // build lumped mass matrix (assumes all positive mass entries, does a simple sum)
    Teuchos::RCP<Thyra::VectorBase<Scalar> > ones = Thyra::createMember(*mass->domain());
    Thyra::assign(ones.ptr(),1.0);

    RCP<Thyra::VectorBase<Scalar> > invLumpMass = Thyra::createMember(*mass->range());
    Thyra::apply(*mass,Thyra::NOTRANS,*ones,invLumpMass.ptr());
    Thyra::reciprocal(*invLumpMass,invLumpMass.ptr());

    invMassMatrix_ = Thyra::diagonal(invLumpMass);
  }
}

template<typename Scalar>
void ExplicitModelEvaluator<Scalar>::
buildArgsPrototypes()
{
  typedef Thyra::ModelEvaluatorBase MEB;
  
  MEB::InArgsSetup<Scalar> inArgs(this->getUnderlyingModel()->createInArgs());
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_alpha,false);
  inArgs.setSupports(MEB::IN_ARG_beta,false);
  inArgs.setSupports(MEB::IN_ARG_x_dot,false);
  prototypeInArgs_ = inArgs;

  MEB::OutArgsSetup<Scalar> outArgs(this->getUnderlyingModel()->createOutArgs());
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_W,false);
  outArgs.setSupports(MEB::OUT_ARG_W_op,false);
  prototypeOutArgs_ = outArgs;
}

template<typename Scalar>
void ExplicitModelEvaluator<Scalar>::
setOneTimeDirichletBeta(double beta,const Thyra::ModelEvaluator<Scalar> & me) const
{
  using Teuchos::Ptr;
  using Teuchos::ptrFromRef;
  using Teuchos::ptr_dynamic_cast;

  // try to extract base classes that support setOneTimeDirichletBeta
  Ptr<const panzer::ModelEvaluator<Scalar> > panzerModel = ptr_dynamic_cast<const panzer::ModelEvaluator<Scalar> >(ptrFromRef(me));
  if(panzerModel!=Teuchos::null) {
    panzerModel->setOneTimeDirichletBeta(beta);
    return;
  }
  else {
    Ptr<const Thyra::EpetraModelEvaluator> epModel = ptr_dynamic_cast<const Thyra::EpetraModelEvaluator>(ptrFromRef(me));
    if(epModel!=Teuchos::null) {
      Ptr<const panzer::ModelEvaluator_Epetra> panzerEpetraModel 
          = ptr_dynamic_cast<const panzer::ModelEvaluator_Epetra>(epModel->getEpetraModel().ptr());

      if(panzerEpetraModel!=Teuchos::null) {
        panzerEpetraModel->setOneTimeDirichletBeta(beta);
        return;
      }
    }
  }

  // if you get here then the ME is not a panzer::ME or panzer::EpetraME, check
  // to see if its a delegator
 
  Ptr<const Thyra::ModelEvaluatorDelegatorBase<Scalar> > delegator
      = ptr_dynamic_cast<const Thyra::ModelEvaluatorDelegatorBase<Scalar> >(ptrFromRef(me));
  if(delegator!=Teuchos::null) {
    setOneTimeDirichletBeta(beta,*delegator->getUnderlyingModel());
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                             "panzer::ExplicitModelEvaluator::setOneTimeDirichletBeta can't find a panzer::ME or a panzer::EpetraME. "
                             "The deepest model is also not a delegator. Thus the recursion failed and an exception was generated.");
}

} // end namespace panzer

#endif
