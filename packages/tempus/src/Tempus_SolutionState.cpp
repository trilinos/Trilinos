//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionState_impl.hpp"

namespace Tempus {

TEMPUS_INSTANTIATE_TEMPLATE_CLASS(SolutionState)

// Nonmember constructor from non-const solution vectors, x.
template Teuchos::RCP<SolutionState<double> > createSolutionStateX(
    const Teuchos::RCP<Thyra::VectorBase<double> >& x,
    const Teuchos::RCP<Thyra::VectorBase<double> >& xdot,
    const Teuchos::RCP<Thyra::VectorBase<double> >& xdotdot);

// Nonmember constructor from const solution vectors, x.
template Teuchos::RCP<SolutionState<double> > createSolutionStateX(
    const Teuchos::RCP<const Thyra::VectorBase<double> >& x,
    const Teuchos::RCP<const Thyra::VectorBase<double> >& xdot,
    const Teuchos::RCP<const Thyra::VectorBase<double> >& xdotdot);

// Nonmember constructor from const solution vectors, x.
template Teuchos::RCP<SolutionState<double> > createSolutionStateME(
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model,
    const Teuchos::RCP<StepperState<double> >& stepperState,
    const Teuchos::RCP<PhysicsState<double> >& physicsState);

}  // namespace Tempus

#endif
