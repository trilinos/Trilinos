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

#include <string>
#include <string_view>

namespace Teuchos {

/** \brief Read a process environment variable by name.
 *
 * \param[in] name Variable name (must not contain '=').
 * \return Pointer to the value, or nullptr if \c name is invalid or the variable
 *         is unset. On MSVC the value is held in per-thread storage inside this
 *         function; on other platforms the pointer refers to the C runtime
 *         environment table (see \c std::getenv).
 *
 * \warning <b>Lifetime.</b> Do not retain the returned pointer across later calls
 *          to \c getEnvironmentVariableValue, \c setEnvironmentVariable,
 *          \c unsetEnvironmentVariable, or other APIs that modify the process
 *          environment. Earlier pointers may be invalidated (see
 *          <a href="https://stackoverflow.com/questions/48568707/getenv-function-may-be-unsafe-really">discussion of getenv safety</a>).
 *          Copy into an owning \c std::string if the value must outlive the call.
 *
 * \warning <b>Concurrency.</b> Since C++11, \c std::getenv is thread-safe only while
 *          the environment is not modified. Concurrent \c setenv, \c unsetenv, or
 *          \c putenv (and the Teuchos wrappers that call them) introduce data races
 *          with reads unless you provide external synchronization.
 *
 * \throws std::runtime_error if the underlying platform read fails (for example
 *         \c ENOMEM from \c _dupenv_s on MSVC).
 *
 * \sa <a href="https://en.cppreference.com/w/cpp/utility/program/getenv">std::getenv</a>,
 *     <a href="https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/dupenv-s-wdupenv-s?view=msvc-170">_dupenv_s</a>
 */
const char* getEnvironmentVariableValue(const char name[]);

/** \brief Overload for \c std::string - avoids requiring <tt>.c_str()</tt> at call sites. */
inline const char* getEnvironmentVariableValue(const std::string& name) {
  return getEnvironmentVariableValue(name.c_str());
}

/** \brief Set a process environment variable.
 *
 * \param[in] name Variable name (must not contain '=').
 * \param[in] value Value to assign.
 * \param[in] overwrite If zero and the variable already exists, the call is a
 *            no-op; otherwise the value is replaced (POSIX \c setenv semantics).
 *
 * \note No-op if \c name or \c value is nullptr or \c name contains '=' (no
 *       exception; POSIX would report \c EINVAL for an invalid name).
 *
 * \warning Modifies the process environment. Not safe to call concurrently with
 *          \c getEnvironmentVariableValue, \c std::getenv, or other readers without
 *          synchronization. May invalidate pointers returned by prior \c std::getenv
 *          calls on platforms that use the C environment table directly.
 *
 * \throws std::runtime_error if the underlying platform call fails (for example
 *         \c ENOMEM on POSIX or a non-zero return from MSVC \c _putenv_s).
 *
 * \sa <a href="https://man7.org/linux/man-pages/man3/setenv.3.html">setenv(3)</a>,
 *     <a href="https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/putenv-s-wputenv-s?view=msvc-170">_putenv_s</a>,
 *     <a href="https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/dupenv-s-wdupenv-s?view=msvc-170">_dupenv_s</a>
 */
void setEnvironmentVariable(const char name[], const char value[], int overwrite);

/** \brief Overload for \c std::string - avoids requiring <tt>.c_str()</tt> at call sites. */
inline void setEnvironmentVariable(const std::string& name, const std::string& value, int overwrite) {
  setEnvironmentVariable(name.c_str(), value.c_str(), overwrite);
}

/** \brief Remove a process environment variable.
 *
 * \param[in] name Variable name (must not contain '=').
 *
 * \note No-op if \c name is nullptr or contains '=' (no exception).
 *
 * \warning Same concurrency and pointer-invalidation considerations as
 *          \c setEnvironmentVariable.
 *
 * \throws std::runtime_error if the underlying platform call fails.
 *
 * \sa <a href="https://pubs.opengroup.org/onlinepubs/9690949399/functions/unsetenv.html">unsetenv</a>,
 *     <a href="https://man7.org/linux/man-pages/man3/setenv.3.html">unsetenv(3)</a>,
 *     <a href="https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/putenv-s-wputenv-s?view=msvc-170">_putenv_s</a> (empty value removes the variable)
 */
void unsetEnvironmentVariable(const char name[]);

/** \brief Overload for \c std::string - avoids requiring <tt>.c_str()</tt> at call sites. */
inline void unsetEnvironmentVariable(const std::string& name) {
  unsetEnvironmentVariable(name.c_str());
}

template <typename T>
T getEnvironmentVariable(std::string_view environmentVariableName,
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
    T &value, bool &initialized, std::string_view environmentVariableName,
    const T defaultValue);

} // namespace Teuchos

#endif
