// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_Stepper.hpp"
#include "Tempus_Stepper_impl.hpp"

namespace Tempus {

  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(Stepper)

  // Validate that the model supports explicit ODE evaluation, f(x,t) [=xdot]
  template void validExplicitODE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model);

  // Validate that the model supports explicit second order ODE evaluation, f(x,xdot,t) [=xdotdot]
  template void validSecondOrderExplicitODE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model);

  // Validate ME supports implicit ODE/DAE evaluation, f(xdot,x,t) [= 0]
  template void validImplicitODE_DAE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model);

  // Validate ME supports 2nd order implicit ODE/DAE evaluation, f(xdotdot,xdot,x,t) [= 0]
  template void validSecondOrderODE_DAE(
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model);

  // Returns the default solver ParameterList for implicit Steppers.
  Teuchos::RCP<Teuchos::ParameterList> defaultSolverParameters();

}

#endif
