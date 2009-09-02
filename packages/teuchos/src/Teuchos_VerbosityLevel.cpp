// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_as.hpp"


namespace {


const Teuchos::Array<Teuchos::EVerbosityLevel>
verbLevelArray = Teuchos::tuple<Teuchos::EVerbosityLevel>(
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
      TEST_FOR_EXCEPT("Should never get here!");
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
