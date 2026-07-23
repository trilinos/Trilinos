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
#include "Tempus_PhiEvaluatorPFD.hpp"
#include "Tempus_PhiEvaluatorPFD_impl.hpp"

namespace Tempus {

TEMPUS_INSTANTIATE_TEMPLATE_CLASS(PhiEvaluatorPFD)

// Nonmember constructor from a ParameterList
template Teuchos::RCP<PhiEvaluatorPFD<double> > createPhiEvaluatorPFD(
	    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

}  // namespace Tempus

#endif
