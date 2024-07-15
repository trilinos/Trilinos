// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_LINEARALGEBRA_HPP
#define ROL_LINEARALGEBRA_HPP

#include <Eigen/Dense>
#include <ostream>

#include "ROL_Ptr.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Stream.hpp"
#include "ROL_Types.hpp"

/** \file  ROL_LinearAlgebra.hpp
    \brief Provides basic capabilities in solving dense
           linear systems and eigenvalue problems using
           Eigen to provide the implementation */


namespace ROL {

namespace LA {

enum ETransp { NO_TRANS, TRANS, CONJ_TRANS };
enum DataAccess { Copy, View };


template<typename Real>
using EVector = Eigen::Matrix<Real,1,Eigen::Dynamic>;

template<typename Real>
using EMatrix = Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic>;

template<typename Real>
class Vector{
  private:
    EVector<Real> v_;
  public:
    Vector(int size) : v_(size) {}
    Vector() : v_() {}
    Real* values(){ return const_cast<Real*>(v_.data()); }
    void resize(int n) { v_.resize(n); }
    void size(int n) { v_.resize(n); v_.setZero();}
    Real& operator()(int i) {
      return v_(i);
    }
    Real& operator[](int i) { return v_[i];}
    Real dot(Vector<Real> x) {
      return v_.dot(x.v_);
    }
    void operator -=(Vector<Real> x) { v_-= x.v_;}
    void operator +=(Vector<Real> x) { v_+= x.v_;}
    void scale(Real alpha) { v_ *= alpha; }
    int numRows() { return v_.size(); }
    int stride() { return v_.outerStride(); }
};

template<typename Real>
class Matrix{
  private:
    EMatrix<Real> M_;
  public:
    Matrix() : M_() {}
    Matrix(int rows, int columns) : M_(rows, columns) {}
    Matrix(DataAccess access, Matrix<Real> A, int rows, int cols, int rowstart = 0, int colstart=0)
    {
      if(access == Copy)
        M_ = A.M_.block(rowstart, colstart, rows, cols).eval();
      else
        M_ = A.M_.block(rowstart, colstart, rows, cols); // Does this DTRT?

    }
    Real* values(){ return const_cast<Real*>(M_.data()); }
    int stride() { return M_.outerStride(); }
    void reshape(int m, int n) { M_.resize(m, n); }
    Real normOne() { return M_.template lpNorm<1>(); }
    Eigen::PartialPivLU<EMatrix<Real>> partialPivLu(bool inplace)
    {
      if(inplace)
        return Eigen::PartialPivLU<Eigen::Ref<EMatrix<Real>>>(M_);
      else
        return M_.partialPivLu();
    }
    Real& operator()(int i, int j) {
      return M_(i, j);
    }

    void multiply	(ETransp transa, ETransp transb, Real alpha, const Matrix&	A,
                   const Matrix& B, Real 	beta) {

      EMatrix<Real> AA;
      if(transa == NO_TRANS)
        AA = A.M_;
      else if(transa == TRANS)
        AA = A.M_.transpose();
      else
        AA = A.M_.conjugate();
      EMatrix<Real> BB;
      if(transa == NO_TRANS)
        BB = B.M_;
      else if(transa == TRANS)
        BB = B.M_.transpose();
      else
        BB = B.M_.conjugate();
      if(beta != Real(1))
        M_ *= beta;
      M_.noalias() += alpha * AA * BB;
    }

    int numRows() { return M_.rows(); }
    int numCols() { return M_.cols(); }

};

//template<typename Real>
//class MatrixFactorization {
//  public:
//    MatrixFactorization(Ptr<LA::Matrix<Real>> A, bool inplace = false);
//    Ptr<LA::Vector<Real>> solve(Ptr<LA::Vector<Real>> b);
//};

//template<typename Real>
//class LuFactorization : public MatrixFactorization<Real>{
//  private:
//    Eigen::PartialPivLU<LA::Matrix<Real>> _lu;
//    Ptr<LA::Matrix<Real>> _A;

//  public:
//    LuFactorization(Ptr<LA::Matrix<Real>> A, bool inplace = false)
//    {
//      _lu = A.partialPivLu(inplace);
//      if(inplace)
//        _A = A; // we need to keep the matrix A alive if we do the decomposition in place.
//    }

//    Ptr<LA::Vector<Real>> solve(Ptr<LA::Vector<Real>> b)
//    {
//      return makePtr<LA::Vector<Real>>(_lu.solve(*b));
//    }

//};

//template<typename Real>
//class LinearSolver {
//private:
  
// Ptr<LA::Matrix<Real>>                A_;
// Ptr<LA::MatrixFactorization<Real>>   P_;
// ParameterList                        options_;

//public:
//  LinearSolver(Ptr<LA::Matrix<Real>> A, ParameterList& opts) : A_(A), options_(opts)
//  {
//    bool inplace = options_.get("Inplace", false);
//    if(options_.get("Factorization", "LU") == "LU")
//      P_ = dynamicPtrCast<MatrixFactorization>(makePtr<LuFactorization>(A_, inplace));
//    else
//      throw Exception::NotImplemented("Only LU factorization implemented for Eigen backend");
//  }

//  void solve(const Ptr<LA::Vector<Real>>&x, const Ptr<LA::Vector<Real>>&b, std::ostream& outStream = std::cout)
//  {
//      auto res = P_.solve(b);
//      *b = *res; // Does this DTRT?
//  };

//  Ptr<LA::Vector<Real>> solve(const Ptr<LA::Vector<Real>>&b, std::ostream& outStream = std::cout)
//  {
//      return P_.solve(b);
//  };

//};

//template<typename> class EigenvalueSolver;

//template<typename Real>
//class EigenvalueProblem {
//private:
 
//  Ptr<LA::Vector<Real>> d_;         // Vector of eigenvalues
//  Ptr<LA::Matrix<Real>> A_,Vl_,Vr_; // Left and right eigenvectors

//public:
 
//  friend class EigenvalueSolver<Real>;

//  EigenvalueProblem();

//  EigenvalueProblem( const Ptr<LA::Matrix<Real>>& A,
//                     const Ptr<LA::Matrix<Real>>& V,	
//                     const Ptr<LA::Vector<Real>>& d );

//  EigenvalueProblem( const Ptr<LA::Matrix<Real>>& A,
//         const Ptr<LA::Matrix<Real>>& Vl,	
//                     const Ptr<LA::Matrix<Real>>& Vr,	
//                     const Ptr<LA::Vector<Real>>& d );

  

//  static Ptr<EigenvaluerProble> create( const Ptr<LA::Matrix<Real>>& A,
//                                        const Ptr<LA::Matrix<Real>>& V,	
//                                        const Ptr<LA::Vector<Real>>& d );

//  static Ptr<EigenvaluerProble> create( const Ptr<LA::Matrix<Real>>& A,
//                                        const Ptr<LA::Matrix<Real>>& Vl,	
//                                        const Ptr<LA::Matrix<Real>>& Vr,	
//                                        const Ptr<LA::Vector<Real>>& d );
  
//};


//template<typename Real>
//class EigenvalueSolver {
//private:

//  Ptr<LA::LinearProblem<Real>> problem_;
//  ParameterList                options_;

//public:
//  void setOptions( ParameterList& opts );
//  void setProblem( const Ptr<LA::EigenvalueProblem<Real>>& problem );
//  void solve( std::ostream& outStream );
//};


} // namespace LA

} // namespace ROL

#endif // ROL_LINEARAPGEBRA_HPP
