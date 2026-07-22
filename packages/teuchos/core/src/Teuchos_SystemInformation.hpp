// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_SYSTEMINFORMATION_HPP
#define TEUCHOS_SYSTEMINFORMATION_HPP

/// \file Teuchos_SystemInformation.hpp
/// \brief Collect information about the runtime environment
///
/// This tool collects the values of environment variables and the output of
/// commands for debugging purposes. Several useful variables and commands are
/// pre-registered. Additional environment variables can be added by setting
/// export
/// TEUCHOS_USER_ENVIRONMENT_VARIABLES=MY_FANCY_VARIABLE;MY_LESS_FANCY_VARIABLE
/// Additional commands can be added by setting
/// export
/// TEUCHOS_USER_COMMANDS=label_for_the_command;command_to_call;executable_that_is_checked_for_availabilty;....
/// The collection of data can be triggered by passing a --print-system-info to
/// any Teuchos::CommandLineProcessor

#include <map>
#include <string>

namespace Teuchos::SystemInformation {

/// Check whether a command is available on the system.
bool commandIsAvailable(const std::string &command);

/// Run a command and capture its output
std::string runCommandAndCaptureOutput(const std::string &command);

enum RegistrationResult { REGISTERED, ALREADY_PRESENT, FAILURE };

/// Register an environment variable that should be tracked.
RegistrationResult registerEnvironmentVariable(const std::string &variableName);

/// Register all variables with a given prefix that can be found in the
/// environment.
void registerAllPrefixedVariables(const std::string &prefix);

/// Register an environment variable prefix that should be tracked.
RegistrationResult registerEnvironmentVariablePrefix(const std::string &prefix);

/// Register a command.
RegistrationResult
registerCommand(const std::string &commandLabel,
                const std::string &commandToRunAndCapture = "",
                const std::string &commandToCheckForExistence = "");

/// Track commonly used environment variables and commands.
void initializeCollection();

/// Collect information about the system.
std::map<std::string, std::string> collectSystemInformation();

} // namespace Teuchos::SystemInformation

#endif
