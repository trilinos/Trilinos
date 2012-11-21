// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER


#include <cmath>

#include "Piro_Epetra_MatrixFreeDecorator.hpp"

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "EpetraExt_ModelEvaluator.h"


Piro::Epetra::MatrixFreeDecorator::MatrixFreeDecorator(
                          Teuchos::RCP<EpetraExt::ModelEvaluator>& model_,
                          double lambda_) :
  model(model_),
  lambda(lambda_)
{
  // Allocate base residual
  fBase = Teuchos::rcp(new Epetra_Vector(*(model->get_f_map())));
}

Piro::Epetra::MatrixFreeDecorator::~MatrixFreeDecorator()
{
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::MatrixFreeDecorator::get_x_map() const
{
  return model->get_x_map();
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::MatrixFreeDecorator::get_f_map() const
{
  return model->get_f_map();
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::MatrixFreeDecorator::get_p_map(int l) const
{
  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::MatrixFreeDecorator::get_g_map(int j) const
{
  return model->get_g_map(j);
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::MatrixFreeDecorator::get_x_init() const
{
  return model->get_x_init();
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::MatrixFreeDecorator::get_p_init(int l) const
{
  return model->get_p_init(l);
}

Teuchos::RCP<Epetra_Operator> Piro::Epetra::MatrixFreeDecorator::create_W() const
{
  return Teuchos::rcp(new Piro::Epetra::MatrixFreeOperator(model, lambda));
}

Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
Piro::Epetra::MatrixFreeDecorator::create_WPrec() const
{
  return model->create_WPrec();
}

EpetraExt::ModelEvaluator::InArgs Piro::Epetra::MatrixFreeDecorator::createInArgs() const
{
  return model->createInArgs();
}

Teuchos::RCP<const Teuchos::Array<std::string> >
  Piro::Epetra::MatrixFreeDecorator::get_p_names(int l) const
{
  return model->get_p_names(l);
}


EpetraExt::ModelEvaluator::OutArgs Piro::Epetra::MatrixFreeDecorator::createOutArgs() const
{
  // Augment outArgs to support W
  // Copy supports from underlying model (can be moved to a gerenal utility)
  typedef EpetraExt::ModelEvaluator EME;

  EME::OutArgsSetup outArgs;
  EME::OutArgs modelOutArgs = model->createOutArgs();

  int np = modelOutArgs.Np();
  int ng = modelOutArgs.Ng();
  outArgs.set_Np_Ng(np, ng);

  for (int i=0; i<EME::NUM_E_OUT_ARGS_MEMBERS; i++) {
   const EME::EOutArgsMembers iarg = static_cast<EME::EOutArgsMembers>(i);
   if (modelOutArgs.supports(iarg)) outArgs.setSupports(iarg);
  }
  
  for (int i=0; i<np; i++) 
   if (!modelOutArgs.supports(OUT_ARG_DfDp, i).none()) {

     outArgs.setSupports(OUT_ARG_DfDp, i, modelOutArgs.supports(OUT_ARG_DfDp, i));
   }
    
  for (int i=0; i<ng; i++) 
   if (!modelOutArgs.supports(OUT_ARG_DgDx, i).none()) {
     outArgs.setSupports(OUT_ARG_DgDx, i, modelOutArgs.supports(OUT_ARG_DgDx, i));
   }
    
  for (int i=0; i<ng; i++) 
   if (!modelOutArgs.supports(OUT_ARG_DgDx_dot, i).none()) {
     outArgs.setSupports(OUT_ARG_DgDx_dot, i, modelOutArgs.supports(OUT_ARG_DgDx_dot, i));
   }

  for (int i=0; i<np; i++) 
    for (int j=0; j<ng; j++) 
     if (!modelOutArgs.supports(OUT_ARG_DgDp, j, i).none()) {
       outArgs.setSupports(OUT_ARG_DgDp, j, i, modelOutArgs.supports(OUT_ARG_DgDp, j, i));
     }
    
  for (int i=0; i<ng; i++) {
    outArgs.setSupports(OUT_ARG_g_sg, i, 
			modelOutArgs.supports(OUT_ARG_g_sg, i));
    outArgs.setSupports(OUT_ARG_DgDx_sg, i, 
			modelOutArgs.supports(OUT_ARG_DgDx_sg, i));
    outArgs.setSupports(OUT_ARG_DgDx_dot_sg, i, 
			modelOutArgs.supports(OUT_ARG_DgDx_dot_sg, i));
    for (int j=0; j<np; j++) 
      outArgs.setSupports(OUT_ARG_DgDp_sg, i, j, 
			  modelOutArgs.supports(OUT_ARG_DgDp_sg, i, j));
  }

  // Add support for W
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties( DerivativeProperties(
      DERIV_LINEARITY_UNKNOWN ,DERIV_RANK_FULL , false));

  return outArgs;
}

void Piro::Epetra::MatrixFreeDecorator::evalModel( const InArgs& inArgs,
                                     const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<Epetra_Operator> W_out = outArgs.get_W();

  if (W_out == Teuchos::null) {
    // Just pass through as is: nothing to Decorate
    model->evalModel(inArgs, outArgs);
  }
  else {

    RCP<Piro::Epetra::MatrixFreeOperator> W_mfo =
      Teuchos::rcp_dynamic_cast<Piro::Epetra::MatrixFreeOperator>(W_out);

    TEUCHOS_TEST_FOR_EXCEPTION(W_mfo==Teuchos::null, std::logic_error, 
      "Epetra_Operator sent as W to Piro::Epetra::MatrixFreeDecorator\n" 
      "be of type Piro::Epetra::MatrixFreeOperator");
   
    // Do base case for MatrixFree: set f instead of W
    OutArgs modelOutArgs(outArgs);
    InArgs modelInArgs(inArgs);

    // Store f_out in case it was also requested
    RCP<Epetra_Vector> f_out = outArgs.get_f();

    modelOutArgs.set_f(fBase);
    modelOutArgs.set_W(Teuchos::null);

    //Evaluate the underlying model
    model->evalModel(modelInArgs, modelOutArgs);

    // If f_out was requested, return it.
    if (f_out != Teuchos::null) *f_out = *fBase;

    // Save unperturbed solution (deep copy inArgs, shallow f)
    InArgs clonedInArgs = inArgs;
    for (int l = 0; l < inArgs.Np(); ++l) {
      const RCP<const Epetra_Vector> p_l = inArgs.get_p(l);
      if (nonnull(p_l))
        clonedInArgs.set_p(l, Teuchos::rcp(new Epetra_Vector(*p_l)));
    }
    clonedInArgs.set_x(Teuchos::rcp(new Epetra_Vector(*inArgs.get_x())));

    bool haveXdot = false;
    if (inArgs.supports(IN_ARG_x_dot)) {
      RCP<const Epetra_Vector> xdot = inArgs.get_x_dot();
      if (nonnull(xdot)) {
        clonedInArgs.set_x_dot(Teuchos::rcp(new Epetra_Vector(*inArgs.get_x_dot())));
        haveXdot = true;
      }
    }
    W_mfo->setBase(clonedInArgs, fBase, haveXdot);
  }
}

/***  Epetra Operator Overloads follow ***/

Piro::Epetra::MatrixFreeOperator::MatrixFreeOperator(
             const Teuchos::RCP<EpetraExt::ModelEvaluator>& model_,
             double lambda_) :
  model(model_),
  modelInArgs(model_->createInArgs()),
  solutionNorm(0.0),
  baseIsSet(false),
  haveXdot(false),
  lambda(lambda_)
{
  using Teuchos::rcp;

  // Save an RCP to a map, preventing it from going out of scope 
  // when references to the map anc comm are later pulled out
  map = model->get_x_map();
 

  // Allocate space for perturbed solution and residuals
  xPert = rcp(new Epetra_Vector(*map));
  fPert = rcp(new Epetra_Vector(*(model->get_f_map())));
  xdotPert = rcp(new Epetra_Vector(*map));
}

void Piro::Epetra::MatrixFreeOperator::setBase(
             const EpetraExt::ModelEvaluator::InArgs modelInArgs_,
             Teuchos::RCP<Epetra_Vector> fBase_,
             const bool haveXdot_)
{
  // Deep copy base solution, shallow copy residual base
  modelInArgs = modelInArgs_;
  haveXdot = haveXdot_;

  fBase = fBase_;
  modelInArgs.get_x()->Norm2(&solutionNorm);
  baseIsSet = true;

}

int  Piro::Epetra::MatrixFreeOperator::Apply
    (const Epetra_MultiVector& V, Epetra_MultiVector& Y) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!baseIsSet, std::logic_error, 
     " Piro::Epetra::MatrixFreeOperator must have Base values set before Apply");

  // Compute  Wv=y  by perturbing x
  //const Epetra_Vector* v = V(0); ??THIS DOESN'T WORK!!! Why???
  Teuchos::RCP<Epetra_Vector> v = Teuchos::rcp(new Epetra_Vector(View, V, 0));
  Teuchos::RCP<Epetra_Vector> y = Teuchos::rcp(new Epetra_Vector(Copy, Y, 0));
  Teuchos::RCP<const Epetra_Vector> xBase = modelInArgs.get_x();
  Teuchos::RCP<const Epetra_Vector> xdotBase;
  if (haveXdot) xdotBase = modelInArgs.get_x_dot();
  double vectorNorm;
  v->Norm2(&vectorNorm);

  // Any operator time zero vector is zero vector
  if (vectorNorm == 0.0) {
    Y.PutScalar(0.0);
    return 0;
  }

  const double eta = lambda * (lambda + solutionNorm/vectorNorm);

  if (haveXdot) {
    const double alpha_W = modelInArgs.get_alpha();
    const double beta_W = modelInArgs.get_beta();

    xPert->Update(1.0, *xBase, beta_W * eta, *v, 0.0);
    xdotPert->Update(1.0, *xdotBase, alpha_W * eta, *v, 0.0);
  } else {
    xPert->Update(1.0, *xBase, eta, *v, 0.0);
  }

  EpetraExt::ModelEvaluator::OutArgs modelOutArgs =
    model->createOutArgs();

  modelInArgs.set_x(xPert);
  if (haveXdot) modelInArgs.set_x_dot(xdotPert);

  // Alert model that this is a perturbed calculation, in case it does something different.
  EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> fPertEval(fPert, EpetraExt::ModelEvaluator::EVAL_TYPE_APPROX_DERIV);
  modelOutArgs.set_f(fPertEval);

  model->evalModel(modelInArgs, modelOutArgs);

  modelInArgs.set_x(xBase);
  if (haveXdot) modelInArgs.set_x_dot(xdotBase);
  modelOutArgs.set_f(fBase);

  y->Update(1.0, *fPert, -1.0, *fBase, 0.0);
  y->Scale(1.0/eta);

  Teuchos::RCP<Epetra_Vector> Y0 = Teuchos::rcp(new Epetra_Vector(View, Y, 0));
  *Y0 = *y; // copy in

  return 0;
}

Piro::Epetra::MatrixFreeOperator::~MatrixFreeOperator() { 
   // relinquish pointers
   modelInArgs.set_x(Teuchos::null); 
   if (haveXdot) modelInArgs.set_x_dot(Teuchos::null); 
   for (int l = 0; l < modelInArgs.Np(); ++l) 
     modelInArgs.set_p(l,Teuchos::null);
}

// The following are trivial implementations
int  Piro::Epetra::MatrixFreeOperator::SetUseTranspose(bool UseTranspose)
{
   TEUCHOS_TEST_FOR_EXCEPTION(UseTranspose, std::logic_error, 
     " Piro::Epetra::MatrixFreeOperator cannot support Transpose operators");
   return 0;
}
int  Piro::Epetra::MatrixFreeOperator::ApplyInverse
    (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
     " Piro::Epetra::MatrixFreeOperator does not support applyInverse");
   return 0;
}
double  Piro::Epetra::MatrixFreeOperator::NormInf() const
{
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
     " Piro::Epetra::MatrixFreeOperator does not support NormInf");
   return 0.0;
}
const char*  Piro::Epetra::MatrixFreeOperator::Label() const
{ return "Piro::Epetra::MatrixFreeOperator"; }
bool  Piro::Epetra::MatrixFreeOperator::HasNormInf() const
{ return false; }
bool  Piro::Epetra::MatrixFreeOperator::UseTranspose() const
{ return false; }
const Epetra_Comm &  Piro::Epetra::MatrixFreeOperator::Comm() const
{ return map->Comm(); }
const Epetra_Map&  Piro::Epetra::MatrixFreeOperator::OperatorDomainMap() const
{ return *map; }
const Epetra_Map&  Piro::Epetra::MatrixFreeOperator::OperatorRangeMap() const
{ return *map; }

