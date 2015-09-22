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

//Computes the exact Schur complement block LU decomposition

template <typename ordinal_type, typename value_type>
Stokhos::BlockPreconditioner<ordinal_type, value_type>::
BlockPreconditioner(
  const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& K_,const ordinal_type p_, const ordinal_type m_) : 
       K(K_),
       p(p_),
       m(m_)
{
}

template <typename ordinal_type, typename value_type>
Stokhos::BlockPreconditioner<ordinal_type, value_type>::
~BlockPreconditioner()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::BlockPreconditioner<ordinal_type, value_type>::
facto(ordinal_type n) const
{
  if (n > 1)
    return (n * facto(n-1));
  else
    return (1);
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::BlockPreconditioner<ordinal_type, value_type>::
siz (ordinal_type n, ordinal_type m) const
{
  //n is the polynomial order and m is the number of random variables
  return (facto(n+m)/(facto(n)*facto(m)));
 }

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::BlockPreconditioner<ordinal_type, value_type>::
ApplyInverse(const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Input,
	     Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Result, 
	     ordinal_type n) const
{ //Solve M*Result=Input
  ordinal_type c=siz(p,m);
  ordinal_type s = siz(p-1,m);
  
  //Split residual
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> r1(Teuchos::Copy, Input, s, 1);
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> r2(Teuchos::Copy, Input, c-s, 1, s, 0);
  
  //Split Result     
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> u1(Teuchos::Copy, Result, s, 1);
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> u2(Teuchos::Copy, Result, c-s, 1, s, 0);
  
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> B(Teuchos::View, K, s, c-s, 0, s);
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> D(Teuchos::View, K, c-s, c-s, s,s);
  
  //rD=inv(D)r2
  
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> Dr(c-s,1);
  
  for (ordinal_type i=0; i<c-s; i++)
    Dr(i,0)=r2(i,0)/D(i,i);
  
  ordinal_type ret = r1.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS, -1.0, B, Dr, 1.0);
  TEUCHOS_ASSERT(ret == 0);
  
  //Compute S=A-B*inv(D)*Bt
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> S(Teuchos::Copy, K, s, s);
  //Compute B*inv(D)
  Teuchos::SerialDenseMatrix<ordinal_type, value_type> BinvD(s,c-s);
  for (ordinal_type i=0; i<c-s; i++) //col
    for (ordinal_type j=0; j<s; j++) //row
      BinvD(j,i)=B(j,i)/D(i,i);
  
  S.multiply(Teuchos::NO_TRANS,Teuchos::TRANS, -1.0, BinvD, B, 1.0);
  
  Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > SS, w, rr;
  SS = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (S));
  w = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (s,1));
  rr = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (r1));
  
  
  // Setup solver
  Teuchos::SerialDenseSolver<ordinal_type, value_type> solver;
  solver.setMatrix(SS);
  solver.setVectors(w, rr);
  //Solve S*w=r1
  if (solver.shouldEquilibrate()) {
    solver.factorWithEquilibration(true);
    solver.equilibrateMatrix();
  }
  solver.solve();
  
  for (ordinal_type i=0; i<s; i++)
    Result(i,0)=(*w)(i,0);
  
  ret = r2.multiply(Teuchos::TRANS,Teuchos::NO_TRANS, -1.0, B, *w, 1.0);
  TEUCHOS_ASSERT(ret == 0);
  
  for (ordinal_type i=s; i<c; i++)
    Result(i,0)=r2(-s+i,0)/D(-s+i, -s+i);

  return 0;
}
