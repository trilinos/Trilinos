/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

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
