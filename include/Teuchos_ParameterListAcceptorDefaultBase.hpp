// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT ParameterListAcceptorDefaultBase : virtual public ParameterListAcceptor {
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
