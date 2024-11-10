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
#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeStepControl_impl.hpp"

namespace Tempus {

TEMPUS_INSTANTIATE_TEMPLATE_CLASS(TimeStepControl)

// Nonmember constructor from ParameterList.
template Teuchos::RCP<TimeStepControl<double> > createTimeStepControl(
    Teuchos::RCP<Teuchos::ParameterList> const& pList, bool runInitialize);

}  // namespace Tempus

#endif
