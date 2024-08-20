// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_MAKEVALIDVERBOSESTREAM_HPP
#define TPETRA_DETAILS_MAKEVALIDVERBOSESTREAM_HPP

#include "TpetraCore_config.h"
#include "Teuchos_FancyOStream.hpp"

namespace Tpetra {
namespace Details {

Teuchos::RCP<Teuchos::FancyOStream>
makeValidVerboseStream (const Teuchos::RCP<Teuchos::FancyOStream>& out =
                        Teuchos::null);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_MAKEVALIDVERBOSESTREAM_HPP
