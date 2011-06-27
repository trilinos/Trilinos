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

#ifndef TEUCHOS_VERBOSE_OBJECT_PARAMETER_LIST_HELPERS_HPP
#define TEUCHOS_VERBOSE_OBJECT_PARAMETER_LIST_HELPERS_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"


namespace Teuchos {


/** \brief Return the sublist of valid parameters for the "VerboseObject"
 * sublist.
 *
 * This function need not be directly called by clients since the function
 * <tt>setupVerboseObjectSublist()</tt> sets up the sublist automatically.
 *
 * \relates VerboseObject
 */
TEUCHOS_LIB_DLL_EXPORT RCP<const ParameterList> getValidVerboseObjectSublist();


/** \brief Setup a sublist called "VerboseObject" in the given parameter list.
 *
 * \param paramList
 *          [in/out] The parameter list hat the "VerboseObject" sublist will
 *          be added to
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>paramList!=0</tt>
 * </ul>
 *
 * \relates VerboseObject
 */
TEUCHOS_LIB_DLL_EXPORT void setupVerboseObjectSublist( ParameterList* paramList );

/** \brief Read the parameters in the "VerboseObject" sublist and set them on
 * the given VerboseObject.
 *
 * \param paramList
 *          [in/out] On input, contains the user's parameter list for the
 *          given objet for which "VerboseObject" can be a sublist of.
 * \param oStream
 *          [out] The oStream object to be used.  On output,
 *          <tt>oStream->get()!=0</tt> if an output stream was specified by
 *          the parameter sublist.
 * \param verbLevel
 *          [out] The verbosity level to be used.  On output,
 *          <tt>*verbLevel</tt> gives the verbosity level set in the parameter
 *          list.  If no verbosity level was set, then a value of
 *          <tt>*verbLevel==VERB_DEFAULT</tt> will be set on return.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>oStream!=0</tt>
 * <li><tt>verbLevel!=0</tt>
 * </ul>
 *
 * \relates VerboseObject
 */
TEUCHOS_LIB_DLL_EXPORT void readVerboseObjectSublist(
  ParameterList* paramList,
  RCP<FancyOStream> *oStream, EVerbosityLevel *verbLevel
  );


/** \brief Read the parameters in the "VerboseObject" sublist and set them on
 * the given VerboseObject.
 *
 * \param paramList
 *          [in/out] On input, contains the user's parameter list for the
 *          given objet for which "VerboseObject" can be a sublist of.
 * \param verboseObject
 *          [in/out] The verbose object that will have its verbosity level
 *          and/or output stream set.
 *
 * This function just calls the above nontemplated
 * <tt>readVerboseObjectSublist()</tt> to validate and and read the verbosity
 * and output stream from the "VerboseObject" sublist.
 *
 * \relates VerboseObject
 */
template<class ObjectType>
void readVerboseObjectSublist(
  ParameterList* paramList, VerboseObject<ObjectType> *verboseObject
  );


} // namespace Teuchos


// /////////////////////////////////
// Implementations


template<class ObjectType>
void Teuchos::readVerboseObjectSublist(
  ParameterList* paramList, VerboseObject<ObjectType> *verboseObject
  )
{
  TEST_FOR_EXCEPT(0==paramList);
  TEST_FOR_EXCEPT(0==verboseObject);
  const EVerbosityLevel bogusVerbLevel = static_cast<EVerbosityLevel>(-50);
  RCP<FancyOStream> oStream = null;
  EVerbosityLevel verbLevel = bogusVerbLevel;
  readVerboseObjectSublist(paramList,&oStream,&verbLevel);
  verboseObject->setOverridingOStream(oStream);
  verboseObject->setOverridingVerbLevel(verbLevel);
}


#endif // TEUCHOS_VERBOSE_OBJECT_PARAMETER_LIST_HELPERS_HPP
