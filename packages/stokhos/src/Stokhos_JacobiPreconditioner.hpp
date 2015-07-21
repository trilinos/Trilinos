/*
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
*/


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

