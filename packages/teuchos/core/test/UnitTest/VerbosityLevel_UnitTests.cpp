// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_VerbosityLevel.hpp"


namespace {


using Teuchos::ArrayView;
using Teuchos::EVerbosityLevel;


TEUCHOS_UNIT_TEST( VerbosityLevel, EVerbosityLevel_size )
{
  TEST_EQUALITY_CONST(Teuchos::EVerbosityLevel_size, 6);
}


TEUCHOS_UNIT_TEST( VerbosityLevel, getValidVerbLevels )
{
  ECHO(ArrayView<const EVerbosityLevel> validVerbLevels = Teuchos::getValidVerbLevels());
  TEST_EQUALITY(validVerbLevels.size(), Teuchos::EVerbosityLevel_size);
  TEST_EQUALITY_CONST(validVerbLevels[0], Teuchos::VERB_DEFAULT);
  TEST_EQUALITY_CONST(validVerbLevels[1], Teuchos::VERB_NONE);
  TEST_EQUALITY_CONST(validVerbLevels[2], Teuchos::VERB_LOW);
  TEST_EQUALITY_CONST(validVerbLevels[3], Teuchos::VERB_MEDIUM);
  TEST_EQUALITY_CONST(validVerbLevels[4], Teuchos::VERB_HIGH);
  TEST_EQUALITY_CONST(validVerbLevels[5], Teuchos::VERB_EXTREME);
}


TEUCHOS_UNIT_TEST( VerbosityLevel, getValidVerbLevelsNames )
{
  ECHO(ArrayView<const std::string> validVerbLevelsNames =
    Teuchos::getValidVerbLevelsNames());
  TEST_EQUALITY(validVerbLevelsNames.size(), Teuchos::EVerbosityLevel_size);
  TEST_EQUALITY_CONST(validVerbLevelsNames[0], "VERB_DEFAULT");
  TEST_EQUALITY_CONST(validVerbLevelsNames[1], "VERB_NONE");
  TEST_EQUALITY_CONST(validVerbLevelsNames[2], "VERB_LOW");
  TEST_EQUALITY_CONST(validVerbLevelsNames[3], "VERB_MEDIUM");
  TEST_EQUALITY_CONST(validVerbLevelsNames[4], "VERB_HIGH");
  TEST_EQUALITY_CONST(validVerbLevelsNames[5], "VERB_EXTREME");
}


TEUCHOS_UNIT_TEST( VerbosityLevel, getValidVerbLevelsNamesRawStrings )
{
  ECHO(ArrayView<const char * const> validVerbLevelsNamesRawStrings =
    Teuchos::getValidVerbLevelsNamesRawStrings());
  TEST_EQUALITY(validVerbLevelsNamesRawStrings.size(), Teuchos::EVerbosityLevel_size);
  TEST_EQUALITY_CONST(std::strcmp(validVerbLevelsNamesRawStrings[0], "VERB_DEFAULT"), 0);
  TEST_EQUALITY_CONST(std::strcmp(validVerbLevelsNamesRawStrings[1], "VERB_NONE"), 0);
  TEST_EQUALITY_CONST(std::strcmp(validVerbLevelsNamesRawStrings[2], "VERB_LOW"), 0);
  TEST_EQUALITY_CONST(std::strcmp(validVerbLevelsNamesRawStrings[3], "VERB_MEDIUM"), 0);
  TEST_EQUALITY_CONST(std::strcmp(validVerbLevelsNamesRawStrings[4], "VERB_HIGH"), 0);
  TEST_EQUALITY_CONST(std::strcmp(validVerbLevelsNamesRawStrings[5], "VERB_EXTREME"), 0);
}


} // namespace
