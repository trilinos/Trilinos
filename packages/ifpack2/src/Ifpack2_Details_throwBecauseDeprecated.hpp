// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_THROWBECAUSEDEPRECATED_HPP
#define IFPACK2_DETAILS_THROWBECAUSEDEPRECATED_HPP

#include "Ifpack2_config.h"

namespace Ifpack2 {
namespace Details {

/// \brief Call this in deprecated methods that will be removed soon;
///   it will throw std::logic_error with a helpful message.
void
throwBecauseDeprecated (const char functionName[]);

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_THROWBECAUSEDEPRECATED_HPP
