// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

