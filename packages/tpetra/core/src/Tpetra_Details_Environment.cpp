#include <cstdlib>
#include <string>
#include <iostream>
#include <algorithm>

#include "Tpetra_Details_Environment.hpp"

namespace Tpetra {
namespace Details {

std::list<std::string>
DefaultEnvironmentVariables::getDefaults ()
{
  return {"TPETRA_DEBUG", "TPETRA_USE_BLAS"};
}

Environment&
Environment::getInstance ()
{
  // The following construction guarantees that theInstance_ will be
  // destroyed and is instantiated only on first use.
  static Environment theInstance_;
  return theInstance_;
}

Environment::Environment ()
{
  // Initialize the instance
  std::string variableValue;
  std::list<std::string>::iterator variableName;
  std::list<std::string> std_vars(DefaultEnvironmentVariables::getDefaults());
  for (variableName = std_vars.begin();
       variableName != std_vars.end();
       ++variableName) {
    // By getting the value, it will be cached
    variableValue = getValue(*variableName);
  }
}

void Environment::cacheVariable(const std::string& variableName,
                                const char* variableValue) {
  // Cache the environment variable
  environCache_[variableName] = variableValue;
}

bool Environment::variableExists(const std::string& variableName) {
  // Check if variable is in the cache
  if (environCache_.find(variableName) != environCache_.end()) {
    // variableName is cached
    return environCache_[variableName] != NULL;
  } else {
    // Variable has not been cached
    const char* variableValue = std::getenv(variableName.c_str());
    cacheVariable(variableName, variableValue);
    return variableValue != NULL;
  }
}

bool Environment::variableIsCached(const std::string& variableName) {
  // Check if variable is in the cache
  if (environCache_.find(variableName) != environCache_.end()) {
    return true;
  } else {
    // Variable has not been cached
    return false;
  }
}

bool Environment::getBooleanValue(const std::string& variableName,
                                  const bool defaultValue) {


  // Get the value of the environment variable variableName.
  std::string variableValue (getValue (variableName, defaultValue ? "TRUE" : "FALSE"));
  std::transform(variableValue.begin(), variableValue.end(),
                 variableValue.begin(), ::toupper);

  if (variableValue.empty() ||
      variableValue == "0"  ||
      variableValue == "NO" ||
      variableValue == "FALSE") {
    return false;
  } else {
    return true;
  }
}

std::string Environment::getValue(const std::string& variableName,
                                  const std::string& defaultValue) {
  // Get the value of the environment variable variableName.
  const char* tmpValue;
  if (variableExists(variableName)) {
    tmpValue = environCache_[variableName];
  } else {
    // First time encountering this variable, get it from the system and cache
    // it
    tmpValue = std::getenv(variableName.c_str());
    cacheVariable(variableName, tmpValue);
  }

  std::string variableValue;
  if (tmpValue != NULL) {
    variableValue = std::string(tmpValue);
  } else {
    // Initialize variableValue to the default
    variableValue = defaultValue;
  }
  return variableValue;
}

// void Environment::setValue(const std::string& variableName,
//                            const std::string& variableValue,
//                            const int overwrite) {
//   // Set the environment variable
//   bool exists = variableExists(variableName);
//   if ((exists && overwrite) || (!exists)) {
//     const char* tmpValue(variableValue.c_str());
//     environCache_[variableName] = tmpValue;
//     setenv(variableName.c_str(), variableValue.c_str(), 1);
//   }
// }

void Environment::clearCache() {
  environCache_.clear();
}

} // namespace Details
} // namespace Tpetra

