// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT RCP<const ParameterList> getValidVerboseObjectSublist();


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
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT void setupVerboseObjectSublist( ParameterList* paramList );

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
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT void readVerboseObjectSublist(
  ParameterList* paramList,
  RCP<FancyOStream> *oStream, EVerbosityLevel *verbLevel
  );


/** \brief Read the parameters in the "VerboseObject" sublist and set them on
 * the given VerboseObject.
 *
 * \param paramList
 *          [in/out] On input, contains the user's parameter list for the
 *          given object of which "VerboseObject" can be a sublist.
 * \param verboseObject
 *          [in/out] The verbose object that will have its verbosity level
 *          and/or output stream set.
 *
 * This function just calls the above nontemplated
 * <tt>readVerboseObjectSublist()</tt> to validate and and read the
 * verbosity and output stream from the "VerboseObject" sublist.
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
  TEUCHOS_TEST_FOR_EXCEPT(0==paramList);
  TEUCHOS_TEST_FOR_EXCEPT(0==verboseObject);
  const EVerbosityLevel bogusVerbLevel = static_cast<EVerbosityLevel>(-50);
  RCP<FancyOStream> oStream = null;
  EVerbosityLevel verbLevel = bogusVerbLevel;
  readVerboseObjectSublist(paramList,&oStream,&verbLevel);
  verboseObject->setOverridingOStream(oStream);
  verboseObject->setOverridingVerbLevel(verbLevel);
}


#endif // TEUCHOS_VERBOSE_OBJECT_PARAMETER_LIST_HELPERS_HPP
