// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_makeValidVerboseStream.hpp"
#include <iostream>

namespace Tpetra {
namespace Details {

Teuchos::RCP<Teuchos::FancyOStream>
makeValidVerboseStream (const Teuchos::RCP<Teuchos::FancyOStream>& out)
{
  if (out.is_null ()) {
    return Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr));
  }
  else {
    return out;
  }
}

} // namespace Details
} // namespace Tpetra
