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

#ifndef TEUCHOS_DESCRIBABLE_HPP
#define TEUCHOS_DESCRIBABLE_HPP

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

/** \brief Verbosity level */
enum EVerbosityLevel {
	VERB_DEFAULT=0  ///< Generate output as defined by the object
	,VERB_LOW=1     ///< Generate only a minimal amount of output
	,VERB_MEDIUM=2  ///< Generate more output
	,VERB_HIGH=3    ///< Generate a high level of output
	,VERB_EXTREME=4 ///< Generate the most output possible
};

/** \brief Base class for all objects that can describe themselves and
 * their current state.
 * 
 * This base class is designed to be a minimally invasive approach for
 * allowing subclasses to optionally provide detailed debug-style
 * information about their current state.  This interface has just two
 * virtual member functions, <tt>describe(void)</tt> and
 * <tt>describe()</tt>, which both have default implementations.  The
 * shorter version of <tt>describe(void)</tt> (which takes no arguments
 * and returns an <tt>std::string</tt> object) is meant for very short
 * descriptions while the longer version of <tt>describe()</tt> takes
 * and returns an <tt>std::ostream</tt> argument and is designed for
 * more detailed formated output.
 *
 * Since both of these <tt>describe()</tt> functions have reasonable
 * default implementations, when a subclass inherits from this base
 * class, no virtual functions need to be overridden to start with.
 * However, when debugging time comes, one or both of these functions
 * should be overridden to provide more useful information.
 *
 * ToDo: Include an example/testing function for a few different use
 * cases to demonstrate how to use this interface properly.
 */
class Describable {
public:

	/// Default value for <tt>verLevel</tt> in <tt>describe()</tt>
	static const EVerbosityLevel   verbLevel_default;
	/// Default value for <tt>leadingIndent</tt> in <tt>describe()</tt>
	static const std::string       leadingIndent_default;
	/// Default value for <tt>indentSpacer</tt> in <tt>describe()</tt>
	static const std::string       indentSpacer_default;

	/** \brief . */
	virtual ~Describable() {}

	/** \name Public virtual member functions */
	//@{

	/** \brief Return a simple description (usually just one line) of this object.
	 *
	 * The default implementation just returns
	 * <tt>typeid(*this).name()</tt> but a subclass can modify this if
	 * needed.  Note that some compilers return a mangled name from
	 * <tt>std::type_info::name()</tt> (e.g. g++ version 3.4.x and
	 * before) that is hard for non-g++ developers to read.  Therefore,
	 * it is usually beneficial to override this function an build a more
	 * human-readable name for a subclass, especially if templating is
	 * used.
	 */
	virtual std::string describe() const;

	/** \brief Print the object using given indentation with some
	 * verbosity level to an <tt>std::ostream</tt> object.
	 *
	 * \param  out   
	 *               [in] The <tt>std::ostream</tt> object that output is sent to.
	 * \param  verbLevel
	 *               [in] Determines the level of verbosity for which the
	 *               the object will be printed.  If <tt>verbLevel==VERB_DEFAULT</tt>
	 *               (which is the default value), then the verbosity level will
	 *               be determined by the <tt>*this</tt> object (i.e. perhaps through the
	 *               <tt>ObjectWithVerbosity</tt> interface).  It is up to <tt>*this</tt>
	 *               how to interpret the level represented by <tt>verbLevel</tt>.
	 *               The default value is <tt>VERB_DEFAULT</tt>.
	 * \param  leadingIndent
	 *               [in] This is the leading indentation space that should be inserted
	 *               before each new line is printed to <tt>out</tt>.  Default value is
	 *               <tt>""</tt> (i.e. no initial indentation).
	 * \param  indentSpacer
	 *               [in] This is the character string to use to create internal indentation.
	 *               This allows a consistent level of indentation for all output.  The
	 *               default value is <tt>"  "</tt> (i.e. two spaces).
	 *
	 * In order for this function to work effectively for independently
	 * developed classes, a general consensus needs be reached as to
	 * what the various verbosity levels represented in
	 * <tt>verbLevel</tt> mean in relation to the amount of output
	 * produced.
	 *
	 * A default implementation of this function is provided that simply
	 * performs:
   \verbatim

   return out << leadingIndent << this->describe() << std::endl; \endverbatim
	 *
	 * A subclass should override this function to provide more
	 * interesting and more useful information about the object.
	 */
	virtual std::ostream& describe(
		std::ostream                &out
		,const EVerbosityLevel      verbLevel     = verbLevel_default
		,const std::string          leadingIndent = leadingIndent_default
		,const std::string          indentSpacer  = indentSpacer_default
		) const;

};

/** \defgroup Teuchos_Describable_Stream_Manipulators Stream Manipulators for Describable Objects.
 *
 * This code allows a client to call the function
 * <tt>Describable::describe(...)</tt> inside of a string of output stream
 * calls.
 */

/** \brief Describable stream manipulator state. */
struct DescribableStreamManipulatorState {
  const Describable          &describable;
	const EVerbosityLevel      verbLevel;
	const std::string          leadingIndent;
	const std::string          indentSpacer;
  DescribableStreamManipulatorState(
    const Describable          &_describable
    ,const EVerbosityLevel     _verbLevel
    ,const std::string         _leadingIndent
    ,const std::string         _indentSpacer
    )
    :describable(_describable)
    ,verbLevel(_verbLevel)
    ,leadingIndent(_leadingIndent)
    ,indentSpacer(_indentSpacer)
    {}
};

/** \brief Describable output stream maniuplator */
inline DescribableStreamManipulatorState describe(
  const Describable          &describable
  ,const EVerbosityLevel     verbLevel     = Describable::verbLevel_default
  ,const std::string         leadingIndent = Describable::leadingIndent_default
  ,const std::string         indentSpacer  = Describable::indentSpacer_default
  )
{
  return DescribableStreamManipulatorState(describable,verbLevel,leadingIndent,indentSpacer);
}

/** \brief Output stream operator for Describable manipulator */
inline
std::ostream& operator<<( std::ostream& os, const DescribableStreamManipulatorState& d )
{
  return d.describable.describe(os,d.verbLevel,d.leadingIndent,d.indentSpacer);
}

} // namespace Teuchos

#endif // TEUCHOS_DESCRIBABLE_HPP
