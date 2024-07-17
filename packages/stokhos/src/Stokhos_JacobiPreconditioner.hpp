// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_JACOBIPRECONDITIONER_HPP
#define STOKHOS_JACOBIPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {

   template <typename ordinal_type, typename value_type>
  class JacobiPreconditioner : 
    public Stokhos::Operator<ordinal_type,value_type> {

  public:

    //! Constructor 
    JacobiPreconditioner(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A_) : A(A_) {}
  
    //! Destructor
    virtual ~JacobiPreconditioner() {}
    
    virtual ordinal_type ApplyInverse(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Input, 
      Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Result, 
      ordinal_type m) const {
      ordinal_type n=Input.numRows();
      Teuchos::SerialDenseMatrix<ordinal_type, value_type> G(Teuchos::Copy,A);
      Teuchos::SerialDenseMatrix<ordinal_type, value_type> z(n,1);
      for (ordinal_type j=0; j<m; j++){
	if (j==0){  // Compute z=D-1r
	  for (ordinal_type i=0; i<n; i++)
	    z(i,0)=Input(i,0)/A(i,i);
	}
	else {
	  //Compute G=invD(-L-U)=I-inv(D)A 
	  for (ordinal_type i=0; i<n; i++){
	    for (ordinal_type j=0; j<n; j++){
	      if (j==i)
		G(i,j)=0;
	      else 
		G(i,j)=-A(i,j)/A(i,i);
	    }
	  }
	  
	  Result.assign(z);
	  //z=Gz+inv(D)r
	  Result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0, G, z, 1.0);
	  
	}
      }

      return 0;
    }
   
  protected:
     const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& A;
  }; // class JacobiPreconditioner

} // namespace Stokhos

#endif // STOKHOS_JACOBIPRECONDITIONER_HPP

