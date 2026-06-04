// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_EnvVariables.hpp"

#include <stddef.h>
#include <string>

#if !defined(_MSC_VER)
#include <cerrno>  // errno
#endif

#include <cstdlib> // std::free
#include <cstring> // std::strchr, std::strerror
#include <memory>  // std::addressof
#include <sstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "Teuchos_StrUtils.hpp"
#include "Teuchos_TestForException.hpp"

namespace Teuchos {

namespace {

bool environmentVariableNameIsValid(const char name[]) {
  return name != nullptr && std::strchr(name, '=') == nullptr;
}

void throwEnvironmentVariablePlatformError(const char operation[],
                                         const char name[],
                                         const std::string& detail) {
  std::ostringstream oss;
  oss << "Teuchos environment variable " << operation << " failed";
  if(name != nullptr) {
    oss << " for \"" << name << "\"";
  }
  oss << ": " << detail;
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, oss.str());
}

// setenv/unsetenv are POSIX; MSVC CRT provides _putenv_s instead.
// See setenv(3), unsetenv, _putenv_s, and _dupenv_s documentation linked from
// Teuchos_EnvVariables.hpp.
// Other Windows toolchains (e.g. MinGW) may still expose setenv/unsetenv.
#if defined(_MSC_VER)
void setEnvironmentVariableImpl(const char* name, const char* value, int overwrite) {
  if(!environmentVariableNameIsValid(name) || value == nullptr) {
    return;
  }
  if(overwrite == 0) {
    char* buf{};
    size_t bufSize{};
    const errno_t dupStatus = _dupenv_s(std::addressof(buf), std::addressof(bufSize), name);
    if(dupStatus != 0) {
      std::ostringstream oss;
      oss << "_dupenv_s returned " << dupStatus;
      throwEnvironmentVariablePlatformError("set", name, oss.str());
    }
    if(buf != nullptr) {
      std::free(buf);
      return;
    }
  }
  const errno_t putStatus = _putenv_s(name, value);
  if(putStatus != 0) {
    std::ostringstream oss;
    oss << "_putenv_s returned " << putStatus;
    throwEnvironmentVariablePlatformError("set", name, oss.str());
  }
}

void unsetEnvironmentVariableImpl(const char* name) {
  if(!environmentVariableNameIsValid(name)) {
    return;
  }
  const errno_t putStatus = _putenv_s(name, "");
  if(putStatus != 0) {
    std::ostringstream oss;
    oss << "_putenv_s returned " << putStatus;
    throwEnvironmentVariablePlatformError("unset", name, oss.str());
  }
}
#else  // !defined(_MSC_VER)
void setEnvironmentVariableImpl(const char* name, const char* value, int overwrite) {
  if(!environmentVariableNameIsValid(name) || value == nullptr) {
    return;
  }
  if(setenv(name, value, overwrite) != 0) {
    const int err = errno;
    std::ostringstream oss;
    oss << "setenv returned -1, errno=" << err;
    if(err != 0) {
      oss << " (" << std::strerror(err) << ")";
    }
    throwEnvironmentVariablePlatformError("set", name, oss.str());
  }
}

void unsetEnvironmentVariableImpl(const char* name) {
  if(!environmentVariableNameIsValid(name)) {
    return;
  }
  if(unsetenv(name) != 0) {
    const int err = errno;
    std::ostringstream oss;
    oss << "unsetenv returned -1, errno=" << err;
    if(err != 0) {
      oss << " (" << std::strerror(err) << ")";
    }
    throwEnvironmentVariablePlatformError("unset", name, oss.str());
  }
}
#endif

#if defined(_MSC_VER)
const char* getEnvironmentVariableValueImpl(const char* name) {
  char* buf{};
  size_t bufSize{};
  const errno_t dupStatus = _dupenv_s(std::addressof(buf), std::addressof(bufSize), name);
  if(dupStatus != 0) {
    std::ostringstream oss;
    oss << "_dupenv_s returned " << dupStatus;
    throwEnvironmentVariablePlatformError("get", name, oss.str());
  }
  if(buf == nullptr) {
    return nullptr;
  }
  thread_local std::string valueStorage;
  valueStorage.assign(buf);
  std::free(buf);
  return valueStorage.c_str();
}
#else  // !defined(_MSC_VER)
const char* getEnvironmentVariableValueImpl(const char* name) {
  // C++11 and later: concurrent std::getenv calls do not data-race only while
  // the process environment is unchanged. setenv/unsetenv/putenv (including
  // Teuchos set/unset helpers) require external synchronization with reads.
  // The returned pointer is owned by the C runtime; copy to std::string if the
  // value must outlive this call or any later environment access.
  return std::getenv(name);
}
#endif

} // namespace

const char* getEnvironmentVariableValue(const char name[]) {
  if(!environmentVariableNameIsValid(name)) {
    return nullptr;
  }
  return getEnvironmentVariableValueImpl(name);
}

namespace {

const char* getEnvironmentVariableValueFromView(std::string_view environmentVariableName) {
  const std::string nameString(environmentVariableName);
  return getEnvironmentVariableValue(nameString.c_str());
}

} // namespace

void setEnvironmentVariable(const char name[], const char value[], int overwrite) {
  setEnvironmentVariableImpl(name, value, overwrite);
}

void unsetEnvironmentVariable(const char name[]) { unsetEnvironmentVariableImpl(name); }

namespace {

enum EnvironmentVariableState {
  EnvironmentVariableIsSet_ON,
  EnvironmentVariableIsSet_OFF,
  EnvironmentVariableIsSet,
  EnvironmentVariableIsNotSet
};

EnvironmentVariableState
environmentVariableState(const std::string &environmentVariableValue) {
  std::string v = StrUtils::allCaps(environmentVariableValue);
  if (v == "1" || v == "YES" || v == "TRUE" || v == "ON")
    // Environment variable is "ON"
    return EnvironmentVariableIsSet_ON;
  else if (v == "0" || v == "NO" || v == "FALSE" || v == "OFF")
    // Environment variable is "OFF"
    return EnvironmentVariableIsSet_OFF;
  // Environment has some other non-boolean value
  return EnvironmentVariableIsSet;
}

void setEnvironmentVariableMap(
    const char environmentVariableName[],
    std::map<std::string, std::map<std::string, bool>> &valsMap,
    const bool defaultValue) {
  using std::map;
  using std::string;
  using std::vector;

  // Set the default value for this variable
  valsMap[environmentVariableName] =
      map<string, bool>{{"DEFAULT", defaultValue}};

  const char *varVal = getEnvironmentVariableValue(environmentVariableName);
  if (varVal == nullptr) {
    // Environment variable is not set, use the default value for any named
    // variants
    return;
  }

  // Variable is not empty.
  const string varStr(varVal);
  // vector<string> names;
  // split(varStr, [&](const string &x) { names.push_back(x); });
  auto names = StrUtils::splitString(varStr, ',');
  for (auto const &name : names) {
    auto state = environmentVariableState(name);
    if (state == EnvironmentVariableIsSet_ON) {
      // Environment variable was set as ENVAR_NAME=[1,YES,TRUE,ON]
      // Global value takes precedence
      valsMap[environmentVariableName]["DEFAULT"] = true;
    } else if (state == EnvironmentVariableIsSet_OFF) {
      // Environment variable was set as ENVAR_NAME=[0,NO,FALSE,OFF]
      // Global value takes precedence
      valsMap[environmentVariableName]["DEFAULT"] = false;
    } else {
      // Environment variable was set as ENVAR_NAME=...:name:...
      // So we set the mapping true for this named variant
      valsMap[environmentVariableName][name] = true;
    }
  }
  return;
}
} // namespace

template <typename T>
T getEnvironmentVariable(std::string_view environmentVariableName,
                         const T defaultValue) {
  const char* varVal = getEnvironmentVariableValueFromView(environmentVariableName);
  if (varVal == nullptr) {
    return defaultValue;
  } else {
    std::stringstream ss(varVal);
    T parsed;
    ss >> parsed;

    TEUCHOS_TEST_FOR_EXCEPTION(
        !ss, std::out_of_range,
        "Environment variable \""
            << environmentVariableName << "\" has a value " << varVal
            << " that cannot be parsed as a " << typeid(T).name() << ".");

    return parsed;
  }
}

// full specialization of bool to preserve historical parsing behavior
template <>
bool getEnvironmentVariable<bool>(
    std::string_view environmentVariableName, const bool defaultValue) {
  const char* varVal = getEnvironmentVariableValueFromView(environmentVariableName);
  bool retVal = defaultValue;
  if (varVal != nullptr) {
    auto state = environmentVariableState(std::string(varVal));
    if (state == EnvironmentVariableIsSet_ON)
      retVal = true;
    else if (state == EnvironmentVariableIsSet_OFF)
      retVal = false;
  }
  return retVal;
}

/*! full specialization of size_t to preserve historical parsing behavior

Parse it as a long long. If negative, return max size_t.
Else, return cast to size_t
*/
template <>
size_t
getEnvironmentVariable<size_t>(std::string_view environmentVariableName,
                               const size_t defaultValue) {
  const char* varVal = getEnvironmentVariableValueFromView(environmentVariableName);
  if (varVal == nullptr) {
    return defaultValue;
  } else {
    long long val = std::stoll(StrUtils::allCaps(varVal));
    if (val < static_cast<long long>(0)) {
      // If negative - user has requested threshold be lifted
      return std::numeric_limits<size_t>::max();
    }
    if (sizeof(long long) > sizeof(size_t)) {
      // It's hard to test this code, but I want to try writing it
      // at least, in case we ever have to run on 32-bit machines or
      // machines with sizeof(long long)=16 and sizeof(size_t)=8.
      constexpr long long maxSizeT =
          static_cast<long long>(std::numeric_limits<size_t>::max());
      TEUCHOS_TEST_FOR_EXCEPTION(
          val > maxSizeT, std::out_of_range,
          "Environment variable \""
              << environmentVariableName << "\" has a value " << val
              << " larger than the largest size_t value " << maxSizeT << ".");
    }
    return static_cast<size_t>(val);
  }
}

bool idempotentlyGetNamedEnvironmentVariableAsBool(
    const char name[], bool &initialized, const char environmentVariableName[],
    const bool defaultValue) {
  static std::map<std::string, std::map<std::string, bool>> namedVariableMap_;
  if (!initialized) {
    setEnvironmentVariableMap(environmentVariableName, namedVariableMap_,
                              defaultValue);
    initialized = true;
  }
  auto thisEnvironmentVariableMap = namedVariableMap_[environmentVariableName];
  auto thisEnvironmentVariable = thisEnvironmentVariableMap.find(name);
  if (thisEnvironmentVariable != thisEnvironmentVariableMap.end())
    return thisEnvironmentVariable->second;
  return thisEnvironmentVariableMap["DEFAULT"];
}

template <typename T>
T idempotentlyGetEnvironmentVariable(
    T &value, bool &initialized, std::string_view environmentVariableName,
    const T defaultValue) {
  if (!initialized) {
    value = getEnvironmentVariable<T>(environmentVariableName, defaultValue);
    initialized = true;
  }
  return value;
}


template std::string getEnvironmentVariable<std::string>(std::string_view, const std::string);

template std::string idempotentlyGetEnvironmentVariable<std::string>(std::string&, bool&, std::string_view, const std::string);
template int idempotentlyGetEnvironmentVariable<int>(int&, bool&, std::string_view, const int);
template unsigned long idempotentlyGetEnvironmentVariable<unsigned long>(unsigned long&, bool&, std::string_view, const unsigned long);
template bool idempotentlyGetEnvironmentVariable<bool>(bool&, bool&, std::string_view, const bool);
template unsigned long long idempotentlyGetEnvironmentVariable<unsigned long long>(unsigned long long&, bool&, std::string_view, const unsigned long long);

} // namespace Teuchos
