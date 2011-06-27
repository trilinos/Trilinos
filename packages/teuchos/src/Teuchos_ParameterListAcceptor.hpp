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

#ifndef TEUCHOS_PARAMETER_LIST_ACCEPTOR_HPP
#define TEUCHOS_PARAMETER_LIST_ACCEPTOR_HPP

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

class ParameterList;
class DependencySheet;
template<class T> class RCP;

/** \brief Base class objects that can accept a parameter list.
 *
 * ToDo: Finish Documentation!
 */
class TEUCHOS_LIB_DLL_EXPORT ParameterListAcceptor {
public:

  /** \brief . */
  virtual ~ParameterListAcceptor();

  //! @name Pure virtual functions that must be overridden in subclasses 
  //@{

  /** \brief Set parameters from a parameter list and return with default values.
   *
   * \param  paramList [in] On input contains the parameters set by the client.
   *                   Note that <tt>*paramList</tt> may have parameters set to their
   *                   default values added while the list is being parsed either right
   *                   away or later.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>paramList.get() != NULL</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->getParameterList().get() == paramList.get()</tt>
   * </ul>
   *
   * This is parameter list is "remembered" by <tt>*this</tt> object until it is
   * unset using <tt>unsetParameterList()</tt>.
   *
   * <b>Note:</b> When this parameter list is passed in it is assumed that the
   * client has finished setting all of the values that they want to set so
   * that the list is completely ready to read (and be validated) by
   * <tt>*this</tt> object.  If the client is to change this parameter list by
   * adding new options or changing the value of current options, the behavior
   * of <tt>*this</tt> object is undefined.  This is because, the object may
   * read the options from <tt>*paramList</tt> right away or may wait to read
   * some options until a later time.  There should be no expectation that if
   * an option is changed by the client that this will automatically be
   * recognized by <tt>*this</tt> object.  To change even one parameter, this
   * function must be called again, with the entire sublist.
   */
  virtual void setParameterList(RCP<ParameterList> const& paramList) = 0;

  /** \brief Get the parameter list that was set using <tt>setParameterList()</tt>.
   */
  virtual RCP<ParameterList> getNonconstParameterList() = 0;

  /** \brief Unset the parameter list that was set using <tt>setParameterList()</tt>.
   *
   * This just means that the parameter list that was set using
   * <tt>setParameterList()</tt> is detached from this object.  This does not
   * mean that the effect of the parameters is undone.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->getParameterList().get() == NULL</tt>
   * </ul>
   */
  virtual RCP<ParameterList> unsetParameterList() = 0;

  //@}

  //! @name Virtual functions with default implementation 
  //@{

  /** \brief Get const version of the parameter list that was set using <tt>setParameterList()</tt>.
   *
   * The default implementation returns:
   \code
   return const_cast<ParameterListAcceptor*>(this)->getParameterList();
   \endcode
   */
  virtual RCP<const ParameterList> getParameterList() const;

  /** \brief Return a const parameter list of all of the valid parameters that
   * <tt>this->setParameterList(...)</tt> will accept.
   *
   * The default implementation returns <tt>Teuchos::null</tt>.
   */
  virtual RCP<const ParameterList> getValidParameters() const;

  /**
   * \brief Rreturn a const Dependency Sheet of all the 
   * dependencies that should be applied to the parameter
   * list return by <tt>this->getValidParameters()</tt>.
   *
   * The default implementation returns <tt>Teuchos::null</tt>.
   */
  virtual RCP<const DependencySheet> getDependencies() const;

  //@}

};

} // end namespace Teuchos

#endif // TEUCHOS_PARAMETER_LIST_ACCEPTOR_HPP
