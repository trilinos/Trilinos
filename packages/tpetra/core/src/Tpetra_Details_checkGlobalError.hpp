// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tpetra_Details_globalError.hpp
/// \brief Declaration of a function that reports a global error.
///
/// \warning This is an implementation detail of Tpetra.  Users may
///   not rely on this header file or any declarations or definitions
///   in it.  They may disappear or change at any time.

#ifndef TPETRA_DETAILS_GLOBALERROR_HPP
#define TPETRA_DETAILS_GLOBALERROR_HPP

#include "TpetraCore_config.h"

namespace Teuchos {
  template<class OrdinalType>
  class Comm;
} // namespace Teuchos

#include <ostream>

namespace Tpetra {
namespace Details {

void
checkGlobalError(std::ostream& globalOutputStream,
                 const bool localSuccess,
                 const char localErrorMessage[],
                 const char globalErrorMessageHeader[],
                 const Teuchos::Comm<int>& comm);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_GLOBALERROR_HPP
