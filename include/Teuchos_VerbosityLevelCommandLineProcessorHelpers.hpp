// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_VERBOSITY_LEVEL_COMMANDLINE_PROCESSOR_HELPERS_HPP
#define TEUCHOS_VERBOSITY_LEVEL_COMMANDLINE_PROCESSOR_HELPERS_HPP

#include "Teuchos_VerbosityLevel.hpp"


namespace Teuchos {


class CommandLineProcessor;


/** \brief Set a verbosity level parameter on a CommandLineProcessor object..
 *
 * \relates CommandLineProcessor
 */
void setVerbosityLevelOption(
  const std::string &optionName,
  EVerbosityLevel *verbLevel,
  const std::string &docString,
  CommandLineProcessor *clp,
  const bool required = false
  );


} // namespace Teuchos


#endif // TEUCHOS_VERBOSITY_LEVEL_COMMANDLINE_PROCESSOR_HELPERS_HPP
