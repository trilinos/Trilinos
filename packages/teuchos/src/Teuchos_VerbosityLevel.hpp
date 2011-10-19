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

#ifndef TEUCHOS_VERBOSITY_LEVEL_HPP
#define TEUCHOS_VERBOSITY_LEVEL_HPP


/*! \file Teuchos_VerbosityLevel.hpp
    \brief .
*/


#include "Teuchos_Assert.hpp"
#include "Teuchos_iostream_helpers.hpp"


namespace Teuchos {


/** \brief Verbosity level.
 *
 * \ingroup teuchos_outputting_grp
 */
enum EVerbosityLevel {
	VERB_DEFAULT=-1,  ///< Generate output as defined by the object
	VERB_NONE=0,      ///< Generate no output
	VERB_LOW=1,       ///< Generate only a minimal amount of output
	VERB_MEDIUM=2,    ///< Generate more output
	VERB_HIGH=3,      ///< Generate a high level of output
	VERB_EXTREME=4    ///< Generate the most output possible
};


/** Needed for serialization KLN 23/09/2010 */
TEUCHOS_ENUM_INPUT_STREAM_OPERATOR(EVerbosityLevel)


/** \brief Return a std::string representation of the verbosity level.
 *
 * \ingroup teuchos_outputting_grp
 */
TEUCHOS_LIB_DLL_EXPORT std::string toString(const EVerbosityLevel verbLevel);


/** \brief Return true if the verbosity level includes the given level.
 *
 * \param  verbLevel
 *           [in] The verbosity level that is in effect.
 * \param  requestedVerbLevel
 *           [in] The verbosity level the client is asking if
 *           is included in <tt>verbLevel</tt>.
 * \param  isDefaultLevel
 *           [in] Set to <tt>true</tt> if the level in
 *           <tt>requestedVerbLevel</tt> is the default verbosity level.  In
 *           this case, if <tt>verbLevel==VERB_DEFAULT</tt>, then this function
 *           will return <tt>true</tt>.  The default value is <tt>false</tt>.
 */
TEUCHOS_LIB_DLL_EXPORT bool includesVerbLevel(
  const EVerbosityLevel verbLevel,
  const EVerbosityLevel requestedVerbLevel,
  const bool isDefaultLevel = false
  );


/** \brief Return an increased or decreased verbosity level.
 *
 * \param  inputVerbLevel
 *           [in] The base verbosity level.
 * \param  numLevels
 *           [in] The number of levels to increase (>0) or decrease (<0).
 *
 * See the function implementation for details on what it does!
 */
TEUCHOS_LIB_DLL_EXPORT EVerbosityLevel incrVerbLevel(
  const EVerbosityLevel inputVerbLevel,
  const int numLevels
  );


} // namespace Teuchos


#endif // TEUCHOS_VERBOSITY_LEVEL_HPP
