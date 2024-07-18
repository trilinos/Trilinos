// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_PARAMETER_LIST_ACCEPTOR_HPP
#define TEUCHOS_PARAMETER_LIST_ACCEPTOR_HPP

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

class ParameterList;
class DependencySheet;
template<class T> class RCP;

/**
 * \brief Interface for objects that can accept a ParameterList.
 *
 * \section Teuchos_ParameterList_Users Summary for users
 *
 * Most users only need to know about two methods:
 *   - setParameterList()
 *   - getValidParameters()
 *
 * setParameterList() lets users set this object's parameters.
 * getValidParameters() returns a default list of parameters,
 * including any documentation and/or validators that the subclass may
 * provide.  If you call setParameterList(), implementations will fill
 * in any missing parameters with their default values, and do
 * validation.  That's all you really need to know!  If you want more
 * details about semantics, though, please read on.
 *
 * \section Teuchos_ParameterList_Semantics Semantics
 *
 * \subsection Teuchos_ParameterList_Semantics_Delta Complete state or delta?
 *
 * This interface does not define the semantics of calling
 * setParametersList() twice with different lists.  For example,
 * suppose that the class \c SomeClass takes two parameters:
 * <ol>
 * <li> An \c int parameter "Integer parameter" </li>
 * <li> A \c bool parameter "Boolean parameter" </li>
 * </ol>
 * The default value of the first is 0, and the default value of the
 * second is false.  In the following code sample, what is the final
 * state of x's parameters?
 * \code
 * SomeClass x;
 * RCP<ParameterList> pl1 = parameterList ();
 * pl1->set ("Integer parameter", 42);
 * // set parameters the first time
 * x.setParameterList (pl1);
 *
 * RCP<ParameterList> pl2 = parameterList ();
 * pl2->set ("Boolean parameter", true);
 * // set parameters the second time
 * x.setParameterList (pl2);
 * \endcode
 * The answer is that we can't tell without knowing more about
 * <tt>SomeClass</tt>.  There are at least two possibilities:
 * <ol>
 * </li> "Integer parameter" is 0, and "Boolean parameter" is true </li>
 * </li> "Integer parameter" is 42, and "Boolean parameter" is true </li>
 * </ol>
 *
 * The first possibility says that the input ParameterList expresses
 * the <i>complete state</i> of the object.  Any missing parameters in
 * subsequent calls get filled in with their default values.  The
 * second possibility says that the input ParameterList expresses a
 * "delta," a difference from its current state.  You must read the
 * subclass' documentation to determine which of these it implements.
 *
 * \section Teuchos_ParameterList_DevNotes Notes for developers
 *
 * Developers who would like a simpler interface from which to inherit
 * may prefer the subclass ParameterListAcceptorDefaultBase.  That
 * class provides default implementations of all but two of this
 * class' methods.
 *
 * It's tempting to begin setParameterList() as follows:
 * \code
 * paramList->validateParametersAndSetDefaults (*getValidParameters ());
 * \endcode
 * That's correct, but be aware that this can only be used to
 * implement "complete state" semantics, not "delta" semantics.
 * This is because validateParametersAndSetDefaults() fills in
 * default values, as its name suggests.
 *
 * Before ParameterList had the validation feature, many
 * implementations of setParameterList() would use the two-argument
 * version of <tt>ParameterList::get()</tt>, and supply the current
 * value of the parameter as the default if that parameter didn't
 * exist in the input list.  This implemented delta semantics.  It is
 * unclear whether implementers knew what semantics they were
 * implementing, but that was the effect of their code.
 *
 * If you want to implement delta semantics, and also want to exploit
 * the validation feature, you have at least two options.  First, you
 * could use the validation method that does not set defaults:
 * \code
 * paramList->validateParameters (*getValidParameters ());
 * \endcode
 * and then use the two-argument version of
 * <tt>ParameterList::get()</tt> in the way discussed above, so that
 * existing parameter values don't get replaced with defaults.
 *
 * The second option is to keep a copy of the ParameterList from the
 * previous call to setParameterList().  This must be a deep copy,
 * because users might have changed parameters since then.  (It is
 * likely that they have just one ParameterList, which they change as
 * necessary.)  You may then use that list -- not the result of
 * getValidParameters() -- as the input argument of
 * validateParametersAndSetDefaults().
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT ParameterListAcceptor {
public:
  //! Destructor.
  virtual ~ParameterListAcceptor ();

  //! @name Pure virtual functions that must be overridden in subclasses
  //@{

  /** \brief Set parameters from a parameter list and return with default values.
   *
   * \param paramList [in/out] On input: contains the parameters set
   *   by the client.  On output: the same list, possibly filled with
   *   default values, depending on the implementation.
   *
   * Implementations of this method generally read parameters out of
   * \c paramList, and use them to modify the state or behavior of
   * this object.  Implementations may validate input parameters, and
   * throw an exception or set an error state if any of them are
   * invalid.  "Validation
   *
   * \pre <tt> ! paramList.is_null () </tt>
   * \post <tt>this->getParameterList().get() == paramList.get()</tt>
   *
   * This object "remembers" \c paramList until it is "unset" using
   * unsetParameterList().  When the input ParameterList is passed in,
   * we assume that the client has finished setting parameters in the
   * ParameterList.  If the client changes \c paramList after calling
   * this method, this object's behavior is undefined.  This is
   * because the object may read the options from \c paramList at any
   * time.  It may either do so in this method, or it may wait to read
   * them at some later time.  Users should not expect that if they
   * change a parameter, that this object will automatically recognize
   * the change.  To change even one parameter, this method must be
   * called again.
   */
  virtual void setParameterList (const RCP<ParameterList>& paramList) = 0;

  /** \brief Get a nonconst version of the parameter list that was set
   *    using setParameterList().
   *
   * The returned ParameterList should be the same object (pointer
   * equality) as the object given to setParameterList().  If
   * setParameterList() has not yet been called on this object, the
   * returned RCP may be null, but need not necessarily be.  If
   * unsetParameterList()
   */
  virtual RCP<ParameterList> getNonconstParameterList() = 0;

  /** \brief Unset the parameter list that was set using <tt>setParameterList()</tt>.
   *
   * This does <i>not</i> undo the effect of setting the parameters
   * via a call to setParameterList().  It merely "forgets" the RCP,
   * so that getParameterList() and getNonconstParameterList() both
   * return null.
   *
   * \post <tt> this->getParameter().is_null () </tt>
   * \post <tt> this->getNonconstParameter().is_null () </tt>
   */
  virtual RCP<ParameterList> unsetParameterList() = 0;

  //@}
  //! @name Virtual functions with default implementation
  //@{

  /** \brief Get const version of the parameter list that was set
   *    using <tt>setParameterList()</tt>.
   *
   * The default implementation returns:
   \code
   return const_cast<ParameterListAcceptor*>(this)->getParameterList();
   \endcode
   */
  virtual RCP<const ParameterList> getParameterList() const;

  /** \brief Return a ParameterList containing all of the valid
   *   parameters that <tt>this->setParameterList(...)</tt> will
   *   accept, along with any validators.
   *
   * Implementations of setParameterList() may use the list returned
   * by getValidParameters() to validate the input ParameterList.
   *
   * The default implementation returns <tt>null</tt>.
   */
  virtual RCP<const ParameterList> getValidParameters() const;

  /**
   * \brief Rreturn a const DependencySheet of all the dependencies
   *   that should be applied to the ParameterList returned by
   *   <tt>this->getValidParameters()</tt>.
   *
   * The default implementation returns <tt>Teuchos::null</tt>.
   */
  virtual RCP<const DependencySheet> getDependencies() const;

  //@}
};

} // end namespace Teuchos

#endif // TEUCHOS_PARAMETER_LIST_ACCEPTOR_HPP
