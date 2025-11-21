// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_SystemInformation.hpp"
#include "Teuchos_EnvVariables.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_TestForException.hpp"

#if !defined(__GNUC__) || (__GNUC__ >= 10)
#include <filesystem>
#endif
#include <iostream>
#include <stdexcept>
#include <set>

// environ should be available on posix platforms
#if not(defined(WIN) && (_MSC_VER >= 1900))
// needs to be in the global namespace
extern char **environ;
#endif
#if defined(_WIN32)
#include <stdio.h>
FILE *popen(const char *command, const char *mode)
{ return _popen(command, mode); }

int pclose(FILE *stream)
{ return _pclose(stream);}
#endif

namespace Teuchos::SystemInformation {

namespace {

std::map<std::string, std::pair<std::string, std::string>>& getCommandsMap() {
  static std::map<std::string, std::pair<std::string, std::string>> commands;
  return commands;
}

std::set<std::string>& getEnvironmentVariablesSet() {
  static std::set<std::string> environmentVariables;
  return environmentVariables;
}

std::set<std::string>& getEnvironmentVariablePrefixSet() {
  static std::set<std::string> environmentVariablePrefixes;
  return environmentVariablePrefixes;
}


}

bool commandIsAvailable(const std::string &command) {
  return !std::system(
      std::string("command -v " + command + " > /dev/null 2>&1").c_str());
}

std::string runCommandAndCaptureOutput(const std::string &command) {
  char buffer[128];
  std::string result = "";
  FILE *pipe = popen(command.c_str(), "r");
  if (!pipe)
    return "command \"" + command + "\"failed";
  try {
    while (fgets(buffer, sizeof buffer, pipe) != NULL) {
      result += buffer;
    }
  } catch (...) {
    pclose(pipe);
    return "command \"" + command + "\"failed";
  }
  pclose(pipe);
  return result;
}

RegistrationResult
registerEnvironmentVariable(const std::string &variableName) {
  auto &environmentVariables = getEnvironmentVariablesSet();

  if (auto search = environmentVariables.find(variableName);
      search != environmentVariables.end()) {
    // variable is already present
    return ALREADY_PRESENT;
  } else {
    // variable not found
    environmentVariables.insert(variableName);
    return REGISTERED;
  }
}

void registerAllPrefixedVariables(const std::string &prefix) {
  char **env;
#if defined(WIN) && (_MSC_VER >= 1900)
  env = *__p__environ();
#else
  env = environ; // defined at the top of this file as extern char **environ;
#endif
  for (; *env; ++env) {
    // const std::string_view ev(*env);

    // split name=value on the first =, everything before = is name
    auto substrings = StrUtils::splitString(*env, '=');
    if (substrings.size() > 0) {
      std::string name = substrings[0];

      if (name.size() >= prefix.size() &&
          name.substr(0, prefix.size()) == prefix) {
        Teuchos::SystemInformation::registerEnvironmentVariable(name);
      }
    }
  }
}

RegistrationResult
registerEnvironmentVariablePrefix(const std::string &prefix) {
  auto &environmentVariablePrefixes = getEnvironmentVariablePrefixSet();
  if (auto search = environmentVariablePrefixes.find(prefix);
      search != environmentVariablePrefixes.end()) {
    // variable is already present
    return ALREADY_PRESENT;
  } else {
    // variable not found
    environmentVariablePrefixes.insert(prefix);
    return REGISTERED;
  }

}

RegistrationResult
registerCommand(const std::string &commandLabel,
                const std::string &commandToRunAndCapture,
                const std::string &commandToCheckForExistence) {
  auto &commands = getCommandsMap();

  std::string myCommandToRunAndCapture = commandToRunAndCapture;
  if (myCommandToRunAndCapture.empty())
    myCommandToRunAndCapture = commandLabel;

  std::string myCommandToCheckForExistence = commandToCheckForExistence;
  if (myCommandToCheckForExistence.empty())
    myCommandToCheckForExistence = myCommandToRunAndCapture;

  if (auto search = commands.find(commandLabel); search != commands.end()) {
    if ((commands[commandLabel].first == myCommandToCheckForExistence) &&
        (commands[commandLabel].second == myCommandToRunAndCapture))
      return ALREADY_PRESENT;
    else {
      std::cerr
          << "Teuchos::SystemInformation: Attempted to register a command (\""
          << commandLabel << "\", \"" << myCommandToRunAndCapture << "\", \""
          << myCommandToCheckForExistence << "\") "
          << "that clashes with already registered command: (" << commandLabel
          << "\", \"" << commands[commandLabel].first << "\", \""
          << commands[commandLabel].second << "\")." << std::endl;
      return FAILURE;
    }
  } else {
    // command not found
    commands[commandLabel] = {myCommandToCheckForExistence, myCommandToRunAndCapture};
    return REGISTERED;
  }
}

void initializeCollection() {

  std::string userProvidedEnvVariables =
      Teuchos::getEnvironmentVariable<std::string>(
          "TEUCHOS_USER_ENVIRONMENT_VARIABLES", "");
  if (!userProvidedEnvVariables.empty()) {
    auto strings = StrUtils::splitString(userProvidedEnvVariables, ';');
    bool isValid = (strings.size() >= 1);
    if (!isValid)
      std::cerr
          << "Teuchos::SystemInformation: The value of the environment "
             "variable "
             "TEUCHOS_USER_ENVIRONMENT_VARIABLES needs to be a semi-colon "
             "seperated string. Value: "
          << userProvidedEnvVariables << std::endl;
    for (int itemNo = 0; itemNo < strings.size(); ++itemNo) {
      std::string &variableName = strings[itemNo];
      if (registerEnvironmentVariable(variableName) == FAILURE)
        isValid = false;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
        !isValid, std::runtime_error,
        "Teuchos::SystemInformation: Invalid environment variable "
        "TEUCHOS_USER_ENVIRONMENT_VARIABLES");
  }

  std::string userProvidedCommands =
      Teuchos::getEnvironmentVariable<std::string>("TEUCHOS_USER_COMMANDS", "");
  if (!userProvidedCommands.empty()) {
    auto strings = StrUtils::splitString(userProvidedEnvVariables, ';');
    bool isValid = (strings.size() % 3 == 0) && (strings.size() >= 3);
    if (!isValid)
      std::cerr << "Teuchos::SystemInformation: The value of the environment "
                   "variable TEUCHOS_USER_COMMANDS "
                   "needs to be a semi-colon seperated string with a number of "
                   "elements that is a multiple of 3. Value: "
                << userProvidedCommands << std::endl;
    int tupleNo = 0;
    while (isValid && (3 * tupleNo + 2 < strings.size())) {
      std::string &commandLabel = strings[3 * tupleNo];
      std::string &commandToRunAndCapture = strings[3 * tupleNo + 1];
      std::string &commandToCheckForExistence = strings[3 * tupleNo + 2];
      if (registerCommand(commandLabel, commandToRunAndCapture,
                          commandToCheckForExistence) == FAILURE) {
        isValid = false;
      }
      ++tupleNo;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!isValid, std::runtime_error,
                               "Teuchos::SystemInformation: Invalid "
                               "environment variable TEUCHOS_USER_COMMANDS");
  }

#if !defined(__GNUC__) || (__GNUC__ >= 10)
  try {
    const std::string executable = std::filesystem::canonical("/proc/self/exe").string();
    if (!executable.empty()) {
      registerCommand("ldd", "ldd " + executable, "ldd");
    }
  } catch (std::filesystem::filesystem_error &) {
    // ignore
  }
#endif

  // CPUs
  registerCommand("lscpu");

  // temperature sensore
  registerCommand("sensors");

  // OpenMP
  registerEnvironmentVariablePrefix("OMP");

  // OpenMPI
  registerCommand("ompi_info");
  registerEnvironmentVariablePrefix("OMPI");

  // MPICH
  registerCommand("mpichinfo");
  registerEnvironmentVariablePrefix("MPICH");

  // CRAY
  registerEnvironmentVariablePrefix("CRAY");

  // modules
  registerCommand("module", "module list 2>&1", "module");

  // CUDA
  registerCommand("nvidia-smi", "nvidia-smi --query", "nvidia-smi");
  registerEnvironmentVariablePrefix("CUDA");

  // ROCm
  registerCommand("rocm-smi", "rocm-smi --showallinfo", "rocm-smi");

  // SYCL
  registerCommand("sycl-ls", "sycl-ls --verbose", "sycl-ls");

  // package namespaced environment variables
  for (auto &prefix : {"TEUCHOS", "KOKKOS", "TPETRA", "STK"})
    registerEnvironmentVariablePrefix(prefix);
}

std::map<std::string, std::string> collectSystemInformation() {

  std::map<std::string, std::string> data;

  const std::string DATA_NOT_AVAILABLE = "NOT AVAILABLE";

  auto &commands = getCommandsMap();
  for (auto &command : commands) {
    const bool isAvailable = commandIsAvailable(command.second.first);
    if (isAvailable) {
      data[command.first] = runCommandAndCaptureOutput(command.second.second);
    } else {
      data[command.first] = DATA_NOT_AVAILABLE;
    }
  }

  auto &environmentVariablePrefixes = getEnvironmentVariablePrefixSet();
  for (auto &prefix : environmentVariablePrefixes) {
    registerAllPrefixedVariables(prefix);
  }

  auto &environmentVariables = getEnvironmentVariablesSet();
  for (auto &envVariable : environmentVariables) {
    const char *varVal = std::getenv(envVariable.c_str());
    if (varVal == nullptr)
      data[envVariable] = "NOT SET";
    else {
      data[envVariable] =
          Teuchos::getEnvironmentVariable<std::string>(envVariable, "");
    }
  }

  return data;
}

} // namespace Teuchos::SystemInformation
