// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_as.hpp"


Teuchos::ArrayView<const Teuchos::EVerbosityLevel> Teuchos::getValidVerbLevels()
{
  static const Tuple<Teuchos::EVerbosityLevel, EVerbosityLevel_size>
    verbLevelArray = tuple<Teuchos::EVerbosityLevel>(
      Teuchos::VERB_DEFAULT,
      Teuchos::VERB_NONE,
      Teuchos::VERB_LOW,
      Teuchos::VERB_MEDIUM,
      Teuchos::VERB_HIGH,
      Teuchos::VERB_EXTREME
      );
  return verbLevelArray();
}


Teuchos::ArrayView<const std::string> Teuchos::getValidVerbLevelsNames()
{
  static const Tuple<std::string, EVerbosityLevel_size>
    verbLevelNamesArray = tuple<std::string>(
      "VERB_DEFAULT",
      "VERB_NONE",
      "VERB_LOW",
      "VERB_MEDIUM",
      "VERB_HIGH",
      "VERB_EXTREME"
      );
  return verbLevelNamesArray();
}


Teuchos::ArrayView<const char * const> Teuchos::getValidVerbLevelsNamesRawStrings()
{
  ArrayView<const std::string> verbLevelNamesArray = getValidVerbLevelsNames();
  static const Tuple<const char*, EVerbosityLevel_size>
    verbLevelNamesRawStringsArray;
  for (int i = 0; i < EVerbosityLevel_size; ++i) {
    verbLevelNamesRawStringsArray[i] = verbLevelNamesArray[i].c_str();
  }
  return verbLevelNamesRawStringsArray();
}


std::string Teuchos::toString(const EVerbosityLevel verbLevel)
{
  switch (verbLevel) {
  case VERB_DEFAULT:
    return "VERB_DEFAULT";
  case VERB_NONE:
    return "VERB_NONE";
  case VERB_LOW:
    return "VERB_LOW";
  case VERB_MEDIUM:
    return "VERB_MEDIUM";
  case VERB_HIGH:
    return "VERB_HIGH";
  case VERB_EXTREME:
    return "VERB_EXTREME";
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(
       true, std::invalid_argument, "Teuchos::toString(const Teuchos::"
       "EVerbosityLevel): Input argument " << verbLevel << " has an invalid "
       "value.  Valid values are VERB_DEFAULT=" << VERB_DEFAULT << ", VERB_NONE"
       "=" << VERB_NONE << ", VERB_LOW=" << VERB_LOW << ", VERB_MEDIUM="
       << VERB_MEDIUM << ", VERB_HIGH=" << VERB_HIGH << ", AND VERB_EXTREME="
       << VERB_EXTREME << ".");
  }

  // NOTE (mfh 15 Sep 2014): Most compilers have figured out that the
  // return statement below is unreachable.  Some older compilers
  // might not realize this.  That's why the return statement was put
  // there, so that those compilers don't warn that this function
  // doesn't return a value.  If it's a choice between one warning and
  // another, I would prefer the choice that produces less code and
  // doesn't have unreachable code (which never gets tested).

  //return ""; // Never get here!
}


bool Teuchos::includesVerbLevel(
  const EVerbosityLevel verbLevel,
  const EVerbosityLevel requestedVerbLevel,
  const bool isDefaultLevel
  )
{
  return (
    ( as<int>(verbLevel) >= as<int>(requestedVerbLevel) )
    ||
    ( verbLevel == VERB_DEFAULT && isDefaultLevel )
    );
}


Teuchos::EVerbosityLevel
Teuchos::incrVerbLevel(
  const EVerbosityLevel inputVerbLevel,
  const int numLevels
  )
{
  if (inputVerbLevel == VERB_DEFAULT)
    return VERB_DEFAULT;
  if (inputVerbLevel == VERB_EXTREME)
    return VERB_EXTREME;
  const int intVerbLevel = as<int>(inputVerbLevel) + numLevels;
  if (intVerbLevel < as<int>(VERB_NONE))
    return VERB_NONE;
  else if (intVerbLevel > as<int>(VERB_EXTREME))
    return VERB_EXTREME;
  // If we get here, then intVerbLevel is a valid verbosity level.
  return getValidVerbLevels()[intVerbLevel];
}
