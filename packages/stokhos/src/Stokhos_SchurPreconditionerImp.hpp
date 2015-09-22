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

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_RCP.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::SchurPreconditioner<ordinal_type,value_type>::
SchurPreconditioner(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& K_, 
  const ordinal_type p_, const ordinal_type m_, const ordinal_type diag_) :
  K(K_),
  p(p_),
  m(m_),
  diag(diag_)
{
}

template <typename ordinal_type, typename value_type>
Stokhos::SchurPreconditioner<ordinal_type,value_type>::
~SchurPreconditioner()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::SchurPreconditioner<ordinal_type,value_type>::
fact (ordinal_type n) const
{
  if (n > 1)
    n = n*fact (n-1);
  else
    n = 1;
  return n;
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::SchurPreconditioner<ordinal_type,value_type>::
size (ordinal_type n, ordinal_type m) const
{
  //n is the polynomial order and m is the number of random variables
  // return (fact(n+m)/(fact(n)*fact(m)));
  ordinal_type min;
  if (n == 0 ) 
    return 1;
  else {
    if (n<=m){
      min = n;
    }
    else {
      min = m;
    }
    
    ordinal_type num = n+m;
    for (ordinal_type i=1; i<=min-1; i++)
      num = num*(n+m-i);
    return num/fact(min); 
  }
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::SchurPreconditioner<ordinal_type,value_type>::
ApplyInverse(const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Input,
	     Teuchos::SerialDenseMatrix<ordinal_type, value_type>& U, 
	     const ordinal_type n) const
{
  //p: polynomial order; m: number of random variables; n: level to stop and form actual Schur Complement; diag=0: Use only diagonal of block D
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Resid(Teuchos::Copy, Input);
  ordinal_type ret;
  ordinal_type lmin;
  if (n<=1)
    lmin=1;  
  else 
    lmin=n;
  
  Teuchos::SerialDenseSolver<ordinal_type, value_type> solver;
  for (ordinal_type l=p; l>=lmin; l--){ 
    ordinal_type c=size(l,m);
    ordinal_type s = size(l-1,m);   
    
    //split residual
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> rminus(Teuchos::View, Resid, s, 1);
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> r(Teuchos::Copy, Resid, c-s, 1, s, 0);
    
    //Compute pre-correction
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> B(Teuchos::View, K, s, c-s, 0, s);   
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> D(Teuchos::View, K, c-s, c-s, s,s); 
    
    //Computing inv(D)r
    if (diag == 0){
      //For D diagonal
      for (ordinal_type i=0; i<c-s; i++)
	r(i,0)=r(i,0)/D(i,i);
      
    }
    else{
      //For full D block
      Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > DD, RR;
      DD = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (Teuchos::Copy, D));
      
      RR = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (r));
      
      
      // Setup solver
      solver.setMatrix(DD);
      solver.setVectors(RR, RR);
      //Solve D*r=r
      if (solver.shouldEquilibrate()) {
	solver.factorWithEquilibration(true);
	solver.equilibrateMatrix();
      }
      solver.solve();
      
      
      for (ordinal_type i=0; i<c-s; i++){
	r(i,0)=(*RR)(i,0);
      }
      
    }
    
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> rm(Teuchos::Copy, rminus);  
    ret = rm.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS, -1.0, B, r, 1.0);
    TEUCHOS_ASSERT(ret == 0);
    
    //set r(l-1)=g(l-1)
    if (l>lmin){
      for (ordinal_type i=0; i<s; i++)
        Resid(i,0)=rm(i,0); 
    }
    else {
      //For l=1, solve A0u0=g0
      if (n<1){
	U(0,0)=rm(0,0)/K(0,0);
      }
      //Compute the Schur completement
      else {
	for (ordinal_type i=0; i<s; i++)
          Resid(i,0)=rm(i,0);
	
	Teuchos::SerialDenseMatrix<ordinal_type, value_type> S(Teuchos::Copy, K, s, s);    
	Teuchos::SerialDenseMatrix<ordinal_type, value_type> BinvD(s,c-s);   
	for (ordinal_type i=0; i<c-s; i++)
          for (ordinal_type j=0; j<s; j++)
	    BinvD(j,i)=B(j,i)/D(i,i);
	
	S.multiply(Teuchos::NO_TRANS,Teuchos::TRANS, -1.0, BinvD, B, 1.0);
	Teuchos::SerialDenseMatrix<ordinal_type, value_type> Ul(Teuchos::View, U, s, 1);
	Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > SS, UU, RR; 
	SS = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (S));
	UU = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (Teuchos::Copy, Ul));
	RR = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (rm));
	//Setup solver
	Teuchos::SerialDenseSolver<ordinal_type, value_type> solver2;
	solver2.setMatrix(SS);
	solver2.setVectors(UU, RR);
	//Solve S*u=rm
	if (solver2.shouldEquilibrate()) {
	  solver2.factorWithEquilibration(true);
	  solver2.equilibrateMatrix();
	}
	solver2.solve();
	
	for (ordinal_type i=0; i<s; i++)
	  U(i,0)=(*UU)(i,0);
      }
    }
  }
  
  for (ordinal_type l=lmin; l<=p; l++){
    //compute post-correction
    ordinal_type c = size(l,m);
    ordinal_type s = size(l-1,m);
    
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> B(Teuchos::View, K, s, c-s, 0, s);
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> D(Teuchos::Copy, K, c-s, c-s, s,s);
    //split residual
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> r(Teuchos::Copy, Resid, c-s, 1, s, 0);
    //Get first s values of U
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> Ul(Teuchos::View, U, s, 1);
    ret = r.multiply(Teuchos::TRANS,Teuchos::NO_TRANS, -1.0, B, Ul, 1.0);
    if (diag == 0) {
      //For diag D
      //Concatenate U
      for (ordinal_type i=s; i<c; i++)
        U(i,0)=r(-s+i,0)/D(-s+i,-s+i);
      
    }
    else {
      //For full block D
      
      Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > DD, RR;
      DD = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (D));
      RR = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (r));
      // Setup solver
      solver.setMatrix(DD);
      solver.setVectors(RR, RR);
      //Solve D*r=r
      if (solver.shouldEquilibrate()) {
	solver.factorWithEquilibration(true);
	solver.equilibrateMatrix();
      }
      solver.solve();
      for (ordinal_type i=s; i<c; i++)
	U(i,0)=(*RR)(-s+i,0); 
    } 
    
    
  }
  
  return 0;
}

 
