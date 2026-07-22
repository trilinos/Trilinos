// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef HERMITIANMATRIX_HPP
#define HERMITIANMATRIX_HPP

#include <initializer_list>

#include "ROL_LinearOperator.hpp"
#include "ROL_ComplexStdVector.hpp"

/** \class HermitianMatrix
    \brief Implementation of a Hermitian Matrix
*/

template<typename Real>
class HermitianMatrix : public ROL::LinearOperator<Real> {
public:

  using value_type = std::complex<Real>;
  using size_type  = typename std::vector<value_type>::size_type;

  HermitianMatrix( size_type N ) {
    std::vector<value_type> temp(N);
    for( size_type i=0u; i<N; ++i ) mat_.emplace_back(temp);   
  }

  HermitianMatrix( std::initializer_list<std::vector<value_type>> m ) :
    mat_(m) {
    bool is_hermitian = true;
    for( size_type row=0u; row<m.size(); ++row ) 
      for( size_type col=row; col<m.size(); ++col ) 
        is_hermitian = is_hermitian && ( mat_[row][col] == std::conj( mat_[col][row] ) );
    if(!is_hermitian) throw std::logic_error("Matrix is not Hermitian");
  }

  const value_type& operator() ( size_type i, size_type j ) const {
    return mat_[i][j];
  }

  value_type& operator() ( size_type i, size_type j ) {
    return mat_[i][j];
  }

  size_type size() const { return mat_.size(); } 

  void apply( ROL::Vector<Real>& Hv, 
              const ROL::Vector<Real>& v, 
              Real& tol ) const override {
    auto& Hvd = ROL::complex_cast(Hv);
    const auto& vd = ROL::complex_cast(v);
    auto N = size();
    for( size_type row = 0; row<N; ++row  ) {
      Hvd[row] = value_type(0,0);
      for( size_type col = 0; col<N; ++col ) 
        Hvd[row] += mat_[row][col]*vd[col];
    }
  }

  void applyAdjoint( ROL::Vector<Real>& Hv, 
                     const ROL::Vector<Real>& v,
                     Real& tol ) const override {
    return apply(Hv,v,tol);
  }

  static HermitianMatrix example_Matrix() {
    value_type a(-0.5,0.0), b(-0.5,0.5), c(-0.5,-0.5);

    HermitianMatrix<Real> A = {{a,b,a,c},
                               {c,a,b,a},
                               {a,c,a,b},
                               {b,a,c,a}};
    return A;
  }

  static std::vector<Real> example_eigenvalues() {
    return std::vector<Real>({-2,-1,0,1});
  }

  static std::vector<ROL::ComplexStdVector<Real>> 
  example_eigenvectors() {

    ROL::ComplexStdVector<Real> v1 = {{ 0.5, 0.0},
                                      { 0.5, 0.0},
                                      { 0.5, 0.0},
                                      { 0.5, 0.0}};     

    ROL::ComplexStdVector<Real> v2 = {{ 0.5, 0.0},
                                      { 0.0,-0.5},
                                      {-0.5, 0.0},
                                      { 0.0, 0.5}};
 
    ROL::ComplexStdVector<Real> v3 = {{ 0.5, 0.0},
                                      {-0.5, 0.0},
                                      { 0.5, 0.0},
                                      {-0.5, 0.0}};

    ROL::ComplexStdVector<Real> v4 = {{ 0.5, 0.0},
                                      { 0,   0.5},
                                      {-0.5, 0.0},
                                      { 0.0,-0.5}};
    return std::vector<ROL::ComplexStdVector<Real>>({v1,v2,v3,v4});
  }

private:
  std::vector<std::vector<value_type>> mat_;
};





template<typename Real>
std::ostream& operator << ( std::ostream& os, const HermitianMatrix<Real>& H ) {
  using size_type = typename HermitianMatrix<Real>::size_type;
  auto N = H.size();
  for( size_type i=0; i<N; ++i ) {
    for( size_type j=0; j<N; ++j ) 
      os << std::setprecision(4) << std::left << std::setw(16) << H(i,j);
    os << std::endl;
  }
  return os;
}

#endif // HERMITIANMATRIX_HPP

