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

#ifndef TEUCHOS_VERBOSITY_LEVEL_HPP
#define TEUCHOS_VERBOSITY_LEVEL_HPP

/*! \file Teuchos_VerbosityLevel.hpp
    \brief .
*/

#include "Teuchos_TestForException.hpp"


namespace Teuchos {


/** \brief Verbosity level.
 *
 * \ingroup teuchos_outputting_grp
 */
enum EVerbosityLevel {
	VERB_DEFAULT=-1  ///< Generate output as defined by the object
	,VERB_NONE=0     ///< Generate no output
	,VERB_LOW=1      ///< Generate only a minimal amount of output
	,VERB_MEDIUM=2   ///< Generate more output
	,VERB_HIGH=3     ///< Generate a high level of output
	,VERB_EXTREME=4  ///< Generate the most output possible
};


/** \brief Return a string representation of the verbosity level.
 *
 * \ingroup teuchos_outputting_grp
 */
inline
std::string toString(const EVerbosityLevel verbLevel)
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


} // namespace Teuchos

#endif // TEUCHOS_VERBOSITY_LEVEL_HPP
