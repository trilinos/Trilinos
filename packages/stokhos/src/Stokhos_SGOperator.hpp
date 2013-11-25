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
