// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_Details_throwBecauseDeprecated.hpp"
#include "Teuchos_TestForException.hpp"
#include <stdexcept>

namespace Ifpack2 {
namespace Details {

void
throwBecauseDeprecated (const char functionName[])
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "\"" << functionName
     << "\" is DEPRECATED and will be removed soon!");
}

} // namespace Details
} // namespace Ifpack2
