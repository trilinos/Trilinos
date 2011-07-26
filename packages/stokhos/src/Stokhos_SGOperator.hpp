// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
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
