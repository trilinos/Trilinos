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

#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_as.hpp"


namespace {


const Teuchos::Array<Teuchos::EVerbosityLevel> verbLevelArray =
  Teuchos::tuple<Teuchos::EVerbosityLevel>(
    Teuchos::VERB_NONE,
    Teuchos::VERB_LOW,
    Teuchos::VERB_MEDIUM,
    Teuchos::VERB_HIGH,
    Teuchos::VERB_EXTREME
    );


} // namespace



std::string Teuchos::toString(const EVerbosityLevel verbLevel)
{
  switch(verbLevel) {
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
      TEUCHOS_TEST_FOR_EXCEPT("Should never get here!");
  }
  return ""; // Never get here!
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
  return verbLevelArray[intVerbLevel];
}
