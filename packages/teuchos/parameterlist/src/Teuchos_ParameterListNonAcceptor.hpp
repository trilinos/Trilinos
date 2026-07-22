// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_PARAMETER_LIST_NON_ACCEPTOR_HPP
#define TEUCHOS_PARAMETER_LIST_NON_ACCEPTOR_HPP

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

namespace Teuchos {


/** \brief Mix-in implementation subclass to be inherited by concrete
 * subclasses who's interface says that they take a parameter list but do not
 * have any parameters yet.
 *
 * ToDo: Finish documention.
 */
class ParameterListNonAcceptor
  : virtual public ParameterListAcceptorDefaultBase
{
public:


  /** \name Overridden from ParameterListAcceptor */
  //@{

  /** \brief Accepts a parameter list but asserts that it is empty. */
  void setParameterList(RCP<ParameterList> const& paramList);

  /** \brief Returns a non-null but empty parameter list.  */
  RCP<const ParameterList> getValidParameters() const;

  //@}

};


} // end namespace Teuchos


#endif // TEUCHOS_PARAMETER_LIST_NON_ACCEPTOR_HPP
