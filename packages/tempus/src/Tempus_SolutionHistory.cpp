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
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_SolutionHistory_impl.hpp"

namespace Tempus {

TEMPUS_INSTANTIATE_TEMPLATE_CLASS(SolutionHistory)

// Nonmember constructor
template Teuchos::RCP<SolutionHistory<double> > createSolutionHistory();

// Nonmember constructor from a ParameterList
template Teuchos::RCP<SolutionHistory<double> > createSolutionHistoryPL(
    Teuchos::RCP<Teuchos::ParameterList> pList);

// Nonmember contructor from a SolutionState.
template Teuchos::RCP<SolutionHistory<double> > createSolutionHistoryState(
    const Teuchos::RCP<SolutionState<double> >& state);

// Nonmember contructor from a Thyra ModelEvaluator.
template Teuchos::RCP<SolutionHistory<double> > createSolutionHistoryME(
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model);

}  // namespace Tempus

#endif
