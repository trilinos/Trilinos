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
#include "Tempus_PhiEvaluatorTaylor.hpp"
#include "Tempus_PhiEvaluatorTaylor_impl.hpp"

namespace Tempus {

TEMPUS_INSTANTIATE_TEMPLATE_CLASS(PhiEvaluatorTaylor)

// Nonmember constructor from a ParameterList
template Teuchos::RCP<PhiEvaluatorTaylor<double> > createPhiEvaluatorTaylor(
	    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

}  // namespace Tempus

#endif
