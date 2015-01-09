//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_FORWARD_SENSITIVITY_STEPPER_TESTER_DEF_H
#define Rythmos_FORWARD_SENSITIVITY_STEPPER_TESTER_DEF_H


#include "Rythmos_ForwardSensitivityStepperTester_decl.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_StepperAsModelEvaluator.hpp"
#include "Thyra_DirectionalFiniteDiffCalculator.hpp"
#include "Thyra_TestingTools.hpp"


template<class Scalar>
Teuchos::RCP<Rythmos::ForwardSensitivityStepperTester<Scalar> >
Rythmos::forwardSensitivityStepperTester()
{
  return Teuchos::rcp(new ForwardSensitivityStepperTester<Scalar>);
}


template<class Scalar>
Teuchos::RCP<Rythmos::ForwardSensitivityStepperTester<Scalar> >
Rythmos::forwardSensitivityStepperTester(const RCP<ParameterList> &paramList)
{
  const RCP<Rythmos::ForwardSensitivityStepperTester<Scalar> > fsst =
    forwardSensitivityStepperTester<Scalar>();
  fsst->setParameterList(paramList);
  return fsst;
}


namespace Rythmos {


// Overridden from ParameterListAcceptor


template<class Scalar>
void ForwardSensitivityStepperTester<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  namespace FSSTU = ForwardSensitivityStepperTesterUtils;
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  this->setMyParamList(paramList);
  errorTol_ = Teuchos::getParameter<double>(*paramList, FSSTU::ErrorTol_name);
}


template<class Scalar>
RCP<const ParameterList>
ForwardSensitivityStepperTester<Scalar>::getValidParameters() const
{
  namespace FSSTU = ForwardSensitivityStepperTesterUtils;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    const RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(FSSTU::ErrorTol_name, FSSTU::ErrorTol_default);
    pl->sublist(FSSTU::FdCalc_name).disableRecursiveValidation().setParameters(
      *Thyra::DirectionalFiniteDiffCalculator<Scalar>().getValidParameters()
      );
    validPL = pl;
  }
  return validPL;
}


template<class Scalar>
bool ForwardSensitivityStepperTester<Scalar>::testForwardSens(
  const Ptr<IntegratorBase<Scalar> > &fwdSensIntegrator
  )
{

  using Teuchos::outArg;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  using Teuchos::sublist;
  typedef Thyra::ModelEvaluatorBase MEB;
  namespace FSSTU = ForwardSensitivityStepperTesterUtils;

  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const RCP<Teuchos::FancyOStream> out = this->getOStream();

  const RCP<ParameterList> paramList = this->getMyNonconstParamList();

  bool success = true;

  // A) Extract out and dynamic cast the forward sens stepper
  RCP<ForwardSensitivityStepper<Scalar> > fwdSensStepper =
    rcp_dynamic_cast<ForwardSensitivityStepper<Scalar> >(
      fwdSensIntegrator->getNonconstStepper()
      );

  // B) Extract out the fwd state stepper
  RCP<StepperBase<Scalar> > stateStepper = fwdSensStepper->getNonconstStateStepper();

  // C) Get the initial condition for the state

  MEB::InArgs<Scalar> state_and_sens_ic = fwdSensStepper->getInitialCondition();
  MEB::InArgs<Scalar> state_ic = 
    extractStateInitialCondition(*fwdSensStepper, state_and_sens_ic);

  // D) Create a StepperAsModelEvalutor for the fwd state calculation

  RCP<StepperAsModelEvaluator<Scalar> >
    stateIntegratorAsModel = stepperAsModelEvaluator(
      stateStepper, fwdSensIntegrator->cloneIntegrator(), state_ic
      );
  //stateIntegratorAsModel->setVerbLevel(verbLevel);

  // E) Compute discrete forward sensitivities using fwdSensIntegrator

  const Scalar finalTime = fwdSensIntegrator->getFwdTimeRange().upper();

  *out << "\nCompute discrete forward sensitivities using fwdSensIntegrator ...\n";

  RCP<const VectorBase<Scalar> > x_bar_final, x_bar_dot_final;
  {
    OSTab tab(out);
    get_fwd_x_and_x_dot( *fwdSensIntegrator, finalTime,
      outArg(x_bar_final), outArg(x_bar_dot_final) );
    *out
      << "\nx_bar_final = x_bar(p,finalTime) evaluated using stateAndSensIntegratorAsModel:\n"
      << describe(*x_bar_final, verbLevel);
  }
  
  // F) Compute discrete forward sensitivities by finite differences

  *out << "\nCompute discrete forward sensitivities by finite differences ...\n";

  RCP<Thyra::MultiVectorBase<Scalar> > DxDp_fd_final;
  {
    OSTab tab(out);

    Thyra::DirectionalFiniteDiffCalculator<Scalar> fdCalc;
    if (nonnull(paramList))
      fdCalc.setParameterList(sublist(paramList, FSSTU::FdCalc_name));
    fdCalc.setOStream(out);
    fdCalc.setVerbLevel(Teuchos::VERB_NONE);
    
    MEB::InArgs<Scalar>
      fdBasePoint = stateIntegratorAsModel->createInArgs();
    
    fdBasePoint.setArgs(state_ic, true); // Set the parameters
    fdBasePoint.set_t(finalTime);
    
    const int g_index = 0;
    const int p_index = fwdSensStepper->getFwdSensModel()->get_p_index();
    
    DxDp_fd_final =
      createMembers(
        stateIntegratorAsModel->get_g_space(g_index),
        stateIntegratorAsModel->get_p_space(p_index)->dim()
        );
    
    typedef Thyra::DirectionalFiniteDiffCalculatorTypes::SelectedDerivatives
      SelectedDerivatives;
    
    MEB::OutArgs<Scalar> fdOutArgs =
      fdCalc.createOutArgs(
        *stateIntegratorAsModel,
        SelectedDerivatives().supports(MEB::OUT_ARG_DgDp, g_index, p_index)
        );
    fdOutArgs.set_DgDp(g_index, p_index, DxDp_fd_final);
    
    // Silence the model evaluators that are called.  The fdCal object
    // will show all of the inputs and outputs for each call.
    stateStepper->setVerbLevel(Teuchos::VERB_NONE);
    stateIntegratorAsModel->setVerbLevel(Teuchos::VERB_NONE);
    
    fdCalc.calcDerivatives(
      *stateIntegratorAsModel, fdBasePoint,
      stateIntegratorAsModel->createOutArgs(), // Don't bother with function value
      fdOutArgs
      );
    
    *out
      << "\nFinite difference DxDp_fd_final = DxDp(p,finalTime): "
      << describe(*DxDp_fd_final, verbLevel);

  }

  // G) Compare the integrated and finite difference sensitivities

  *out << "\nChecking that integrated DxDp(p,finalTime) and finite-diff DxDp(p,finalTime) are similar ...\n";
  
  {
    
    Teuchos::OSTab tab(out);
    
    RCP<const Thyra::VectorBase<Scalar> >
      DxDp_vec_final = Thyra::productVectorBase<Scalar>(x_bar_final)->getVectorBlock(1);
    
    RCP<const Thyra::VectorBase<Scalar> >
      DxDp_fd_vec_final = Thyra::multiVectorProductVector(
        rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> >(
          DxDp_vec_final->range()
          ),
        DxDp_fd_final
        );

    const bool result = Thyra::testRelNormDiffErr(
      "DxDp_vec_final", *DxDp_vec_final,
      "DxDp_fd_vec_final", *DxDp_fd_vec_final,
      "maxSensError", errorTol_,
      "warningTol", 1.0, // Don't warn
      &*out, verbLevel
      );
    if (!result)
      success = false;
    
  }

  return success;

}


template<class Scalar>
ForwardSensitivityStepperTester<Scalar>::ForwardSensitivityStepperTester()
  : errorTol_(ForwardSensitivityStepperTesterUtils::ErrorTol_default)
{}


} // namespace Rythmos


// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define FORWARD_SENSITIVITY_STEPPER_INSTANT(SCALAR) \
  \
  template class ForwardSensitivityStepperTester<SCALAR >; \
  \
  template RCP<ForwardSensitivityStepperTester<SCALAR > > \
  forwardSensitivityStepperTester(); \
  \
  template RCP<ForwardSensitivityStepperTester<SCALAR> > \
  forwardSensitivityStepperTester(const RCP<ParameterList> &paramList);


#endif //Rythmos_FORWARD_SENSITIVITY_STEPPER_TESTER_DEF_H
