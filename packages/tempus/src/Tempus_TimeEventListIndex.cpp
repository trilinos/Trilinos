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
#include "Tempus_TimeEventListIndex.hpp"
#include "Tempus_TimeEventListIndex_impl.hpp"

namespace Tempus {

TEMPUS_INSTANTIATE_TEMPLATE_CLASS(TimeEventListIndex)

// Nonmember constructor
template Teuchos::RCP<TimeEventListIndex<double> > createTimeEventListIndex(
    Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus

#endif
