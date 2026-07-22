// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SG_OPERATOR_HPP
#define STOKHOS_SG_OPERATOR_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_EpetraOperatorOrthogPoly.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

namespace Stokhos {

  /*! 
   * \brief An abstract class to represent a generic stochastic Galerkin 
   * operator as an Epetra_Operator.
   */
  class SGOperator : public virtual Epetra_Operator {
  public:

    //! Constructor
    SGOperator() {}

    //! Destructor
    virtual ~SGOperator() {}

    //! Setup operator
    virtual void setupOperator(
      const Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly >& poly) = 0;

    //! Get SG polynomial
    virtual Teuchos::RCP< Stokhos::EpetraOperatorOrthogPoly > 
    getSGPolynomial() = 0;

    //! Get SG polynomial
    virtual Teuchos::RCP<const Stokhos::EpetraOperatorOrthogPoly > 
    getSGPolynomial() const = 0;

  }; // class SGOperator

} // namespace Stokhos

#endif // STOKHOS_SG_OPERATOR_HPP
