// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tpetra_Details_printOnce.hpp
/// \brief Declaration of Tpetra::Details::printOnce.

#ifndef TPETRA_DETAILS_PRINTONCE_HPP
#define TPETRA_DETAILS_PRINTONCE_HPP

#include "TpetraCore_config.h"
#include <ostream>
#include <string>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
// Forward declaration of Comm.
template <class OrdinalType> class Comm;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {
namespace Details {
  
/// \brief Print on one process of the given communicator, or at least
///   try to do so (if MPI is not initialized).
///
/// \param out [out] Output stream to which to print.  If MPI is
///   initialized, then it need only be valid on Process 0 of the
///   given communicator.  Otherwise, it must be valid on all
///   processes of the given communicator.
///
/// \param s [in] String to print.
///
/// \param comm [in] Communicator; if nullptr, print on all processes,
///   else, print based on above rule.
void
printOnce (std::ostream& out,
	   const std::string& s,
	   const Teuchos::Comm<int>* comm);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_PRINTONCE_HPP
