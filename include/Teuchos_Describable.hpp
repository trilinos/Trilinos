// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_DESCRIBABLE_HPP
#define TEUCHOS_DESCRIBABLE_HPP

#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_LabeledObject.hpp"


namespace Teuchos {


/** \brief Base class for all objects that can describe themselves.
 *
 * This base class is designed to be a minimally invasive approach for
 * letting subclasses provide detailed debug-style information about
 * their current state.  This interface has just two virtual member
 * functions, describe() and description(), which both have default
 * implementations.  The shorter version description() generates a
 * brief one-line description, while the longer version describe() is
 * designed for more detailed multiline formatted output.
 *
 * Since both of these functions have reasonable default
 * implementations, when a subclass inherits from this base class, no
 * virtual functions need to be overridden to start with.  However,
 * when debugging time comes, one or both of these functions should be
 * overridden to provide more useful information.
 *
 * This interface derives from the LabeledObject interface.
 * Therefore, a user may set an object-specific label on every
 * Describable object, and the label may be incorporated into the
 * description of the object.
 *
 * \ingroup teuchos_outputting_grp
 */
class TEUCHOSCORE_LIB_DLL_EXPORT Describable : virtual public LabeledObject {
public:
  //! Default value for the \c verbLevel argument of describe().
  static const EVerbosityLevel verbLevel_default;

  //! @name Public virtual member functions
  //@{

  /** \brief Return a simple one-line description of this object.
   *
   * The default implementation just returns <tt>typeName(*this)</tt>,
   * along with the object's label if defined.  The function
   * <tt>typeName(*this)</tt> guarantees that a demangled,
   * human-readable name is returned on most platforms.  Even if
   * subclasses choose to override this function, this default
   * implementation can still be called as
   * <tt>Teuchos::Describable::description()</tt> in order to print
   * the label name along with the class name.
   */
  virtual std::string description() const;

  /** \brief Print the object with some verbosity level to a FancyOStream.
   *
   * \param out [in/out] The output stream to which to write.
   * \param verbLevel [in] Verbosity level for printing the object.
   *   If <tt>verbLevel==VERB_DEFAULT</tt> (which is the default
   *   value), then the verbosity level will be determined by the
   *   object itself (e.g., through the ObjectWithVerbosity
   *   interface).  It is up to each subclass instance how to
   *   interpret the input verbosity level.
   *
   * You may wrap an std::ostream in a Teuchos::FancyOStream by
   * including "Teuchos_FancyOStream.hpp" and calling
   * Teuchos::getFancyOStream().  For example:
   * \code
   * using Teuchos::RCP;
   * using Teuchos::rcpFromRef;
   * using Teuchos::FancyOStream;
   *
   * // Wrap std::cout in a FancyOStream.
   * RCP<FancyOStream> wrappedCout = getFancyOStream (rcpFromRef (std::cout));
   *
   * // Wrap an output file in a FancyOStream.
   * // Use an RCP rather than a reference, to ensure that the file
   * // isn't closed until the method is completely done with it.
   * RCP<std::ofstream> outFile (new std::ofstream ("myFile.txt"));
   * RCP<FancyOStream> wrappedFile = getFancyOStream (outFile);
   * \endcode
   *
   * A subclass' implementation of this method should begin by tabbing
   * one increment using the OSTab class.  That way, users need not
   * include tabs themselves.
   *
   * This class does not mandate what the different verbosity levels
   * mean.  It would be a good idea for Trilinos developers to reach a
   * consensus on this.  For example, if the subclass implements a
   * large data structure like a sparse matrix, VERB_NONE should print
   * nothing at all, VERB_LOW should print \f$ O(1)\f$ data relative to
   * the dimensions and number of entries of the matrix, and
   * VERB_EXTREME may print all the entries of the matrix.  The
   * subclass must decide how to interpret verbosity levels in between
   * these extremes.
   *
   * This class provides a default implementation of this method that
   * simply performs:
   * \code
   * OSTab tab (out);
   * return out << this->description () << std::endl;
   * \endcode
   * Subclasses should override this method to provide more useful
   * information about the object.
   */
  virtual void
  describe (FancyOStream &out,
            const EVerbosityLevel verbLevel = verbLevel_default) const;

  /// \brief Version of describe() that takes an std::ostream instead
  ///   of a FancyOStream.
  ///
  /// This method just wraps \c out in a FancyOStream and calls the
  /// other describe() method.  We provide this as a convenience to
  /// users who can't remember or don't want to make a FancyOStream
  /// out of their std::ostream.
  ///
  /// Subclasses may not override this method.  Instead, they must
  /// override the version of describe() that takes a FancyOStream.
  void
  describe (std::ostream &out,
            const EVerbosityLevel verbLevel = verbLevel_default) const;

  //! Destructor (marked virtual for memory safety of derived classes).
  virtual ~Describable ();
};


// Describable stream manipulator state class
//
// This is not a class that a user needs to see and that is why it is not
// being given doxygen documentation!
struct DescribableStreamManipulatorState {
  const Describable &describable;
  const EVerbosityLevel verbLevel;
  DescribableStreamManipulatorState(
    const Describable &_describable,
    const EVerbosityLevel _verbLevel = VERB_MEDIUM
    )
    :describable(_describable)
    ,verbLevel(_verbLevel)
    {}
};


/** \brief Describable output stream manipulator.
 *
 * This simple function allows you to insert output from
 * <tt>Describable::describe()</tt> right in the middle of a chain of
 * insertion operations.  For example, you can write:

 \code

  void someFunc( const Teuchos::Describable &obj )
  {
    ...
    std::cout
      << "The object is described as "
      << describe(obj,Teuchos::VERB_MEDIUM);
    ...
  }

 \endcode

 * \relates Describable
 */
inline DescribableStreamManipulatorState describe(
  const Describable &describable,
  const EVerbosityLevel verbLevel = Describable::verbLevel_default
  )
{
  return DescribableStreamManipulatorState(describable,verbLevel);
}


/** \brief Output stream operator for Describable manipulator.
 *
 * To call this function use something like:

 \code

  void someFunc( const Teuchos::Describable &obj )
  {
    ...
    std::cout
      << "The object is described as "
      << describe(obj,Teuchos::VERB_MEDIUM);
    ...
  }

 \endcode

 * Note: The input <tt>std::ostream</tt> is casted to a <tt>FancyOStream</tt>
 * object before calling <tt>Describable::describe()</tt> on the underlying
 * <tt>Describable</tt> object.  There is no way around this since this
 * function must be written in terms of <tt>std::ostream</tt> rather than
 * <tt>FancyOStream</tt> if one is to write compound output statements
 * involving primitive data types.
 *
 * \relates Describable
 */
inline
std::ostream& operator<<(
  std::ostream& os, const DescribableStreamManipulatorState& d
  )
{
  d.describable.describe(*getFancyOStream(Teuchos::rcp(&os,false)),d.verbLevel);
  return os;
}

//
// RAB: Note: The above function works with an std::ostream object even
// through Describable::describe(...) requires a FancyOStream object.  We must
// write the stream manipulator in terms of std::ostream, or compound output
// statements like:
//
//  void foo( FancyOStream &out, Describable &d, EVerbLevel verbLevel )
//    {
//      out << "\nThis is the describable object d:" << describe(d,verbLevel);
//    }
//
// will not work correctly.  The problem is that the first output
//
//   out << "\nThis is the describable object d:"
//
//  must return a reference to an std::ostream object.  This should mean that
//  the next statement, which is basically:
//
//    static_cast<std::ostream&>(out) << DescribableStreamManipulatorState
//
// should not even compile.  However, under gcc 3.4.3, the code did compile
// but did not call the above function.  Instead, it set up some type of
// infinite recursion that resulted in a segfault due to the presence of the
// Teuchos::any class!
//


} // namespace Teuchos

#endif // TEUCHOS_DESCRIBABLE_HPP
