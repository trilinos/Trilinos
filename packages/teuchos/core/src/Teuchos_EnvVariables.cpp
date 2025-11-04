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

#include <cstdlib> // std::getenv
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Teuchos_StrUtils.hpp"
#include "Teuchos_TestForException.hpp"

namespace Teuchos {

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
  using std::getenv;
  using std::map;
  using std::string;
  using std::vector;

  // Set the default value for this variable
  valsMap[environmentVariableName] =
      map<string, bool>{{"DEFAULT", defaultValue}};

  const char *varVal = getenv(environmentVariableName);
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
T getEnvironmentVariable(const std::string_view environmentVariableName,
                         const T defaultValue) {
  const char *varVal = std::getenv(environmentVariableName.data());
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
    const std::string_view environmentVariableName, const bool defaultValue) {
  const char *varVal = std::getenv(environmentVariableName.data());
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
getEnvironmentVariable<size_t>(const std::string_view environmentVariableName,
                               const size_t defaultValue) {
  const char *varVal = std::getenv(environmentVariableName.data());
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
    T &value, bool &initialized, const std::string_view environmentVariableName,
    const T defaultValue) {
  if (!initialized) {
    value = getEnvironmentVariable<T>(environmentVariableName, defaultValue);
    initialized = true;
  }
  return value;
}


template std::string getEnvironmentVariable<std::string>(std::string_view, const std::string);

template std::string idempotentlyGetEnvironmentVariable<std::string>(std::string&, bool&, const std::string_view, const std::string);
template int idempotentlyGetEnvironmentVariable<int>(int&, bool&, const std::string_view, const int);
template unsigned long idempotentlyGetEnvironmentVariable<unsigned long>(unsigned long&, bool&, const std::string_view, const unsigned long);
template bool idempotentlyGetEnvironmentVariable<bool>(bool&, bool&, const std::string_view, const bool);
template unsigned long long idempotentlyGetEnvironmentVariable<unsigned long long>(unsigned long long&, bool&, const std::string_view, const unsigned long long);

} // namespace Teuchos
