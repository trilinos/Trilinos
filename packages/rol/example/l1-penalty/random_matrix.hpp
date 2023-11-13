// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#pragma once
#ifndef RANDOM_MATRIX_HPP
#define RANDOM_MATRIX_HPP

#include <random>
#include "Teuchos_SerialDenseMatrix.hpp"
#include "ROL_Teuchos_LinearOperator.hpp"

template<typename Ordinal,typename Real, typename Distribution, typename Generator>
ROL::Ptr<Teuchos::SerialDenseMatrix<Ordinal,Real>>
random_matrix( Ordinal rows, Ordinal cols, Distribution dist, Generator gen ) {
  auto M = ROL::makePtr<Teuchos::SerialDenseMatrix<Ordinal,Real>>(rows,cols);
  for(Ordinal i=0; i<rows; ++i) 
    for(Ordinal j=0; j<cols; ++j) 
      (*M)(i,j) = dist(gen);
  return M; 
}

template<typename Ordinal,typename Real>
ROL::Ptr<Teuchos::SerialDenseMatrix<Ordinal,Real>>
randn_matrix( Ordinal rows, Ordinal cols, Real mu=0.0, Real sigma=1.0 ) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<Real> dist{mu,sigma};
  return random_matrix<Ordinal,Real>(rows,cols,dist,gen); 
}

template<typename Ordinal, typename Real>
ROL::Ptr<ROL::LinearOperator<Real>>
create_random_operator( Ordinal rows, Ordinal cols ) {
  auto A = randn_matrix<Ordinal,Real>(rows,cols);
  return ROL::makePtr<ROL::TeuchosLinearOperator<Ordinal,Real>>(A);
}

template<typename Ordinal, typename Real>
ROL::Ptr<ROL::Vector<Real>>
create_random_vector( Ordinal n ) {
  ROL::Ptr<ROL::Vector<Real>> v = ROL::makePtr<ROL::TeuchosVector<Ordinal,Real>>(n);
  v->randomize(-1.0,1.0);
  return v;
}


#endif // RANDOM_MATRIX_HPP

