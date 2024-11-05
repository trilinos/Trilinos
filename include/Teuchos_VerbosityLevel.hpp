// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_VERBOSITY_LEVEL_HPP
#define TEUCHOS_VERBOSITY_LEVEL_HPP


/*! \file Teuchos_VerbosityLevel.hpp
    \brief .
*/


#include "Teuchos_Assert.hpp"
#include "Teuchos_ArrayView.hpp"
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


/** Number of valid EVerbosityLevel enum vaules */
constexpr const int EVerbosityLevel_size = 6;


/** Get static array view of verbosity levels enums.
 *
 * \ingroup teuchos_outputting_grp
 */
TEUCHOSCORE_LIB_DLL_EXPORT
ArrayView<const EVerbosityLevel> getValidVerbLevels();


/** Get static array view of verbosity levels string names as array of 
 *  std::string strings.
 *
 * \ingroup teuchos_outputting_grp
 */
TEUCHOSCORE_LIB_DLL_EXPORT
ArrayView<const std::string> getValidVerbLevelsNames();


/** Get static array view of verbosity levels string names as array of 
 *  null-terminated strings.
 *
 * \ingroup teuchos_outputting_grp
 */
TEUCHOSCORE_LIB_DLL_EXPORT
ArrayView<const char * const> getValidVerbLevelsNamesRawStrings();


/** Needed for serialization KLN 23/09/2010 */
TEUCHOS_ENUM_INPUT_STREAM_OPERATOR(EVerbosityLevel)


/** \brief Return a std::string representation of the verbosity level.
 *
 * \ingroup teuchos_outputting_grp
 */
TEUCHOSCORE_LIB_DLL_EXPORT std::string toString(const EVerbosityLevel verbLevel);


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
TEUCHOSCORE_LIB_DLL_EXPORT bool includesVerbLevel(
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
TEUCHOSCORE_LIB_DLL_EXPORT EVerbosityLevel incrVerbLevel(
  const EVerbosityLevel inputVerbLevel,
  const int numLevels
  );


} // namespace Teuchos


#endif // TEUCHOS_VERBOSITY_LEVEL_HPP
