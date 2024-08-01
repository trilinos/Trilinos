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
#ifdef TEMPUS_ENABLE_EPETRA_STACK
#include "CDR_Model.hpp"
#include "CDR_Model_impl.hpp"
#endif
#ifdef TEMPUS_ENABLE_TPETRA_STACK
#include "CDR_Model_Tpetra.hpp"
#include "CDR_Model_Tpetra_impl.hpp"
#endif

namespace Tempus_Test {
#ifdef TEMPUS_ENABLE_EPETRA_STACK
TEMPUS_INSTANTIATE_TEMPLATE_CLASS(CDR_Model)
#endif

#ifdef TEMPUS_ENABLE_TPETRA_STACK
// Get default Tpetra template types
using SC   = Tpetra::Vector<>::scalar_type;
using LO   = Tpetra::Vector<>::local_ordinal_type;
using GO   = Tpetra::Vector<>::global_ordinal_type;
using Node = Tpetra::Vector<>::node_type;

TEMPUS_INSTANTIATE_TEMPLATE_CLASS_TPETRA(CDR_Model_Tpetra, SC, LO, GO, Node)
#endif
}  // namespace Tempus_Test

#endif  // HAVE_TEMPUS_EXPLICIT_INSTANTIATION
