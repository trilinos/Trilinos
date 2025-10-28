// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_ENVVARIABLES_HPP
#define TEUCHOS_ENVVARIABLES_HPP

#include <string_view>

namespace Teuchos {

template <typename T>
T getEnvironmentVariable(const std::string_view environmentVariableName,
                         const T defaultValue);

bool idempotentlyGetNamedEnvironmentVariableAsBool(
    const char name[], bool &initialized, const char environmentVariableName[],
    const bool defaultValue);

/** \brief Read a variable from the environment.
    Example usage:

    \code
    constexpr bool defaultValue = true;
    static bool value = defaultValue;
    static bool initialized = false;
    idempotentlyGetEnvironmentVariable(value, initialized, "TEUCHOS_VARIABLE", defaultValue);

    \endcode
 */
template <typename T>
T idempotentlyGetEnvironmentVariable(
    T &value, bool &initialized, const std::string_view environmentVariableName,
    const T defaultValue);

} // namespace Teuchos

#endif
