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

#ifndef TEUCHOS_PARAMETER_LIST_ACCEPTOR_DEFAULT_BASE_HPP
#define TEUCHOS_PARAMETER_LIST_ACCEPTOR_DEFAULT_BASE_HPP

#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"


namespace Teuchos {


/** \brief Intermediate node base class for objects that accept parameter lists
 * that implements some of the needed behavior automatically.
 *
 * Subclasses just need to implement <tt>setParameterList()</tt> and
 * <tt>getValidParameters()</tt>.  The underlying parameter list is accessed
 * using the non-virtual protected members <tt>setMyParamList()</tt> and
 * <tt>getMyParamList()</tt>.
 */
class ParameterListAcceptorDefaultBase : virtual public ParameterListAcceptor {
public:

  /** \name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const ParameterList> getParameterList() const;

  //@}

protected:

  /** \name Protected accessors to actual parameter list object. */
  //@{

  /** \brief . */
  void setMyParamList( const RCP<ParameterList> &paramList );

  /** \brief . */
  RCP<ParameterList> getMyNonconstParamList();

  /** \brief . */
  RCP<const ParameterList> getMyParamList() const;

  //@}

private:

  RCP<ParameterList> paramList_;

};


//
// Inline definitions
//


inline
void ParameterListAcceptorDefaultBase::setMyParamList(
  const RCP<ParameterList> &paramList
  )
{
  paramList_ = paramList;
}


inline
RCP<ParameterList>
ParameterListAcceptorDefaultBase::getMyNonconstParamList()
{
  return paramList_;
}


inline
RCP<const ParameterList>
ParameterListAcceptorDefaultBase::getMyParamList() const
{
  return paramList_;
}


} // end namespace Teuchos


#endif // TEUCHOS_PARAMETER_LIST_ACCEPTOR_DEFAULT_BASE_HPP
