// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_DIAGPRECONDITIONER_HPP
#define STOKHOS_DIAGPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {
  
  template <typename ordinal_type, typename value_type>
  class DiagPreconditioner : 
    public Stokhos::Operator<ordinal_type,value_type> {
  public:

    //! Constructor 
    DiagPreconditioner(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A_) : A(A_) {}
    
    //! Destructor
    virtual ~DiagPreconditioner() {} 
  
    virtual ordinal_type ApplyInverse(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Input, 
      Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Result, 
      ordinal_type m) const {
      ordinal_type n=Input.numRows();
      for (ordinal_type i=0; i<n; i++){
	Result(i,0)=Input(i,0)/A(i,i);
      }
      return 0;
    }
   
  protected:
    const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A;
  }; // class DiagPreconditioner

} // namespace Stokhos

#endif // STOKHOS_DIAGPRECONDITIONER_HPP

