// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_GSPRECONDITIONER_HPP
#define STOKHOS_GSPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class GSPreconditioner : 
    public Stokhos::Operator<ordinal_type,value_type> {
  public:

    //! Constructor 
    GSPreconditioner(
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type> & A_, 
      const ordinal_type sym_) : A(A_), sym(sym_) {}
    
    //! Destructor
    virtual ~GSPreconditioner() {}
    
    // m is the number of GS iterations to solve Mz=r
    // If sym=0 then do symmetric Gauss Seidel
    virtual ordinal_type ApplyInverse(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Input, 
      Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Result, 
      ordinal_type m) const {
      Result.assign(Input);
      ordinal_type n=A.numRows();
      ordinal_type info;
      Teuchos::LAPACK<ordinal_type, value_type> lapack;
      
      //Get lower triangular part of A  
      Teuchos::SerialDenseMatrix<ordinal_type, value_type> L(A);
      for (ordinal_type i=0; i<n; i++){
	for (ordinal_type j=0; j<n; j++){
	  if (j>i) 
	    L(i,j)=0;
	}
      }
      
      if (sym==0){
  	//Get inv(diag(A))=D
  	Teuchos::SerialDenseMatrix<ordinal_type, value_type> D(n,n);
  	for (ordinal_type i=0; i<n; i++){
	  for (ordinal_type j=0; j<n; j++){
	    D(i,i)=1/A(i,i);
	  }
  	}
	
  	//Get upper triangular part of A  
  	Teuchos::SerialDenseMatrix<ordinal_type, value_type> U(A);
  	for (ordinal_type i=0; i<n; i++){
	  for (ordinal_type j=0; j<n; j++){
	    if (i>j)
	      U(i,j)=0;
	  }
	}
  
 	Result.assign(Input);
	
   	//compute M=(L+D)inv(diagA)
   	Teuchos::SerialDenseMatrix<ordinal_type, value_type> M(n,n);
   	M.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0, L, D, 0.0);
        
   	//Forward solve Lz=r
    	lapack.TRTRS('L', 'N', 'N', M.numRows(), 1, M.values(), M.stride(), Result.values(), Result.stride(),&info);
	
   	//Backward solve Uw=z
    	lapack.TRTRS('U', 'N', 'N', U.numRows(), 1, U.values(), U.stride(), Result.values(), Result.stride(),&info);
   	
      }
      else{
	lapack.TRTRS('L', 'N', 'N', L.numRows(), 1, L.values(), L.stride(), Result.values(), Result.stride(),&info);
      }
      
      return 0;
    }
   
  protected:
    const Teuchos::SerialDenseMatrix<ordinal_type,value_type> & A;

    const ordinal_type sym;
  }; // class GSPreconditioner

} // namespace Stokhos

#endif // STOKHOS_GSPRECONDITIONER_HPP

