//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluator_impl_hpp
#define Tempus_PhiEvaluator_impl_hpp

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluator_decl.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

#include "Thyra_DefaultInverseLinearOp_decl.hpp"
#include "Thyra_DefaultMultipliedLinearOp_decl.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

template <class Scalar>
PhiEvaluator<Scalar>::PhiEvaluator()
  : name_("Phi Evaluator"), isInitialized_(false)
{
}

template <class Scalar>
PhiEvaluator<Scalar>::PhiEvaluator(std::string name)
  : isInitialized_(false)
{
  setName(name);
}

template <class Scalar>
void PhiEvaluator<Scalar>::copy(
    Teuchos::RCP<const PhiEvaluator<Scalar> > phi)
{
  this->setName(phi->getName());
}

template <class Scalar>
std::string PhiEvaluator<Scalar>::description() const
{
  return ("Tempus::PhiEvaluator - '" + name_ + "'");
}

template <class Scalar>
void PhiEvaluator<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  //TODO
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << "\n--- " << this->description() << " ---" << std::endl;

  if ((Teuchos::as<int>(verbLevel) ==
       Teuchos::as<int>(Teuchos::VERB_DEFAULT)) ||
      (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW))) {
    //*l_out << "  abc     = " << ... <<
    // std::endl;
  }

  //if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_MEDIUM)) {
  //
  //}
  *l_out << std::string(this->description().length() + 8, '-') << std::endl;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
PhiEvaluator<Scalar>::getValidParameters() const
{
  return this->getValidParametersBasic();
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
PhiEvaluator<Scalar>::getValidParametersBasic() const
{
  //TODO
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Phi Evaluator");

  pl->setName(getName());

  //pl->set(
  //    "var", default,
  //    "'var' sets the var.  "
  //    "'opt1' - will do this.  "
  //    "'opt2' - will do that!");

  pl->set(
      "PhiEvaluator Type", "PFD",
      "Method to approximate the phi-function evaluation.");

  //pl->set("?", *member_->getNonconstParameterList());

  return pl;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
PhiEvaluator<Scalar>::getNonconstParameterList()
{
  return Teuchos::rcp_const_cast<Teuchos::ParameterList>(getValidParameters());
}

template <class Scalar>
void PhiEvaluator<Scalar>::initialize() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      appModel_ == Teuchos::null, std::logic_error,
      "Error - PhiEvaluator::initialize() Model not set!\n");

  isInitialized_ = true;  // Only place where this is set to true!
}

template <class Scalar>
void PhiEvaluator<Scalar>::checkInitialized()
{
  if (!this->isInitialized()) {
    this->describe(*(this->getOStream()), Teuchos::VERB_MEDIUM);
    TEUCHOS_TEST_FOR_EXCEPTION(
        !this->isInitialized(), std::logic_error,
        "Error - " << this->description() << " is not initialized!");
  }
}

template <class Scalar>
void PhiEvaluator<Scalar>::setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<double> > appModel)
{
  appModel_ = appModel;

  // TODO: read this in from the xml file
  const bool lumpmass = true;
  phiLinSolv_ = Teuchos::rcp(new PhiLinearSolver<Scalar>(appModel_, lumpmass));
}

template <class Scalar>
void PhiEvaluator<Scalar>::setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs)
{

}


/*
 * PhiLinearSolver methods
 */

template <class Scalar>
void PhiLinearSolver<Scalar>::computeMassMatrix(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs)
{

  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;
  using Thyra::createMember;

  // first allocate space for the mass matrix
  fullMassMatrix_ = appModel_->create_W_op();

  // request only the mass matrix from the physics
  // Model evaluator builds: alpha*u_dot + beta*F(u) = 0
  MEB::InArgs<Scalar>  inArgs_new  = appModel_->createInArgs();
  inArgs_new.setArgs(inArgs);
  inArgs_new.set_x_dot(inArgs.get_x_dot());
  inArgs_new.set_alpha(1.0);
  inArgs_new.set_beta(0.0);

  // TODO:
  // set the one time beta to ensure dirichlet conditions
  // are correctly included in the mass matrix: do it for
  // both epetra and Tpetra.
  //if(panzerModel_!=Teuchos::null)
  //  panzerModel_->setOneTimeDirichletBeta(1.0);
  //else {
    // assuming the underlying model is a delegator, walk through
    // the decerator hierarchy until you find a panzer::ME or panzer::EpetraME.
    // If you don't find one, then throw because you are in a load of trouble anyway!
  //  setOneTimeDirichletBeta(1.0,*this->getUnderlyingModel());
  //}

  // set only the mass matrix
  MEB::OutArgs<Scalar> outArgs = appModel_->createOutArgs();
  outArgs.set_W_op(fullMassMatrix_);

  // this will fill the mass matrix operator
  appModel_->evalModel(inArgs_new, outArgs);

  //TODO:
  //Teuchos::RCP<const Epetra_CrsMatrix> crsMat = Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*fullMassMatrix));
  //EpetraExt::RowMatrixToMatrixMarketFile("fullMassMatrix_mat.mm",*crsMat);

  if(!lumpMass_) {
    invMassMatrix_ = Thyra::inverse<Scalar>(appModel_->get_W_factory(),fullMassMatrix_);
  }
  else {
    // build lumped mass matrix (assumes all positive mass entries, does a simple sum)
    Teuchos::RCP<Thyra::VectorBase<Scalar> > ones = Thyra::createMember(*fullMassMatrix_->domain());
    Thyra::assign(ones.ptr(),1.0);

    RCP<Thyra::VectorBase<Scalar> > lumpMass = Thyra::createMember(*fullMassMatrix_->range());
    RCP<Thyra::VectorBase<Scalar> > invLumpMass = Thyra::createMember(*fullMassMatrix_->range());
    Thyra::apply(*fullMassMatrix_, Thyra::NOTRANS, *ones, lumpMass.ptr());

    //TODO:
    //Teuchos::RCP<const Epetra_Vector> mv = Teuchos::rcp_dynamic_cast<const Epetra_Vector>(Thyra::get_Epetra_Vector(crsMat->RangeMap(), invLumpMass));
    //EpetraExt::VectorToMatrixMarketFile("fullMassMatrix_v.mm", *mv);

    Thyra::reciprocal(*lumpMass, invLumpMass.ptr());

    lumpMassMatrix_ = Thyra::diagonal(lumpMass);
    invMassMatrix_ = Thyra::diagonal(invLumpMass);
  }
}

template <class Scalar>
void PhiLinearSolver<Scalar>::applyMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Mf, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f)
{
  // apply the mass matrix
  if (f != Teuchos::null && Mf != Teuchos::null) {
    if (!lumpMass_)
      Thyra::apply(*fullMassMatrix_, Thyra::NOTRANS, *f, Mf.ptr());
    else
      Thyra::apply(*lumpMassMatrix_, Thyra::NOTRANS, *f, Mf.ptr());
  }
}

template <class Scalar>
void PhiLinearSolver<Scalar>::solveMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> f, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf)
{
  // invert the mass matrix
  if (Mf != Teuchos::null && f != Teuchos::null) {
    Thyra::apply(*invMassMatrix_, Thyra::NOTRANS, *Mf, f.ptr());
  }
}

template <class Scalar>
void PhiLinearSolver<Scalar>::computeJacobian(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs)
{

}

}  // namespace Tempus
#endif  // Tempus_PhiEvaluator_impl_hpp


/**

#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

#include "Panzer_ExplicitModelEvaluator.hpp"

#include "PanzerDiscFE_config.hpp"

namespace panzer {

template<typename Scalar>
ExplicitModelEvaluator<Scalar>::
ExplicitModelEvaluator(const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > & model,
                       bool constantMassMatrix,
                       bool useLumpedMass,
                       bool applyMassInverse)
   : Thyra::ModelEvaluatorDelegatorBase<Scalar>(model)
   , constantMassMatrix_(constantMassMatrix)
   , massLumping_(useLumpedMass)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  this->applyMassInverse_ = applyMassInverse;

  // extract a panzer::ModelEvaluator if appropriate
  panzerModel_ = rcp_dynamic_cast<panzer::ModelEvaluator<Scalar> >(model);

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

  MEB::InArgs<Scalar> under_inArgs = under_me->createInArgs();
  under_inArgs.setArgs(inArgs);

  // read in the supplied time derivative in case it is needed by the explicit model (e.g. for stabilization) and make sure alpha is set to zero
  under_inArgs.set_x_dot(inArgs.get_x_dot());
  under_inArgs.set_alpha(0.0);

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

  // build the mass matrix
  if(invMassMatrix_==Teuchos::null || constantMassMatrix_==false)
    buildInverseMassMatrix(inArgs);

  if(f!=Teuchos::null)
      Thyra::Vt_S(scrap_f_.ptr(),-1.0);

  // invert the mass matrix
  if(f!=Teuchos::null && this->applyMassInverse_) {
    Thyra::apply(*invMassMatrix_,Thyra::NOTRANS,*scrap_f_,f.ptr());
  }
  else if(f!=Teuchos::null){
    Thyra::V_V(f.ptr(),*scrap_f_);
  }
}

template<typename Scalar>
void ExplicitModelEvaluator<Scalar>::
buildInverseMassMatrix(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs) const
{
}

template<typename Scalar>
void ExplicitModelEvaluator<Scalar>::
buildArgsPrototypes()
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgsSetup<Scalar> inArgs(this->getUnderlyingModel()->createInArgs());
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_alpha,true);
  inArgs.setSupports(MEB::IN_ARG_beta,true);
  inArgs.setSupports(MEB::IN_ARG_x_dot,true);
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
**/
