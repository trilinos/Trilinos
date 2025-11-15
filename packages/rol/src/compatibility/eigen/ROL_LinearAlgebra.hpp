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
    EVector<Real> v_;  // Owned data
    Eigen::Map<EVector<Real>> map_;  // View of external data
    mutable bool using_map_;  // Flag to track which member to use
    
  public:
    Vector(int size) : v_(size), map_(v_.data(), size), using_map_(false) {}
    Vector() : v_(), map_(nullptr, 0), using_map_(false) {}
    
    Vector(const Vector& other) : v_(), map_(nullptr, 0), using_map_(false) {
      if (other.using_map_) {
        // Copy data from the map into owned storage
        v_ = other.map_;
        map_ = Eigen::Map<EVector<Real>>(v_.data(), v_.size());
      } else {
        v_ = other.v_;  
        map_ = Eigen::Map<EVector<Real>>(v_.data(), v_.size());
      }
    }
    
    Vector(DataAccess access, Real* data, int size) : v_(), map_(nullptr, 0), using_map_(true) {
      if (access == View) {
        // Create map to external data
        map_ = Eigen::Map<EVector<Real>>(data, size);
      } else {
        // Copy data into owned storage
        v_.resize(size);
        using_map_ = false;
        map_ = Eigen::Map<EVector<Real>>(v_.data(), size);
        std::copy(data, data + size, v_.data());
      }
    }

    Real* values() { 
      return using_map_ ? const_cast<Real*>(map_.data()) : const_cast<Real*>(v_.data()); 
    }
    
    void resize(int n) { 
      if (!using_map_) {
        v_.resize(n); 
        map_ = Eigen::Map<EVector<Real>>(v_.data(), n);
      }
      // Cannot resize a mapped vector
    }
    
    void size(int n) { 
      if (!using_map_) {
        v_.resize(n); 
        v_.setZero();
        map_ = Eigen::Map<EVector<Real>>(v_.data(), n);
      }
    }
    
    Real& operator()(int i) {
      return using_map_ ? map_(i) : v_(i);
    }
    
    Real& operator[](int i) { 
      return using_map_ ? map_[i] : v_[i];
    }
    
    Real dot(const Vector<Real>& x) const {
      auto& this_ref = using_map_ ? map_ : v_;
      auto& x_ref = x.using_map_ ? x.map_ : x.v_;
      return this_ref.dot(x_ref);
    }
    
    void operator -=(const Vector<Real>& x) { 
//      auto& this_ref = using_map_ ? map_ : v_;
//      const auto& x_ref = x.using_map_ ? x.map_ : x.v_;
//      this_ref -= x_ref;
      map_ -= x.map_;
    }
    
    void operator +=(const Vector<Real>& x) { 
//      auto& this_ref = using_map_ ? map_ : v_;
//      const auto& x_ref = x.using_map_ ? x.map_ : x.v_;
//      this_ref += x_ref;
      map_ += x.map_;
    }
    
    void scale(Real alpha) { 
      if (using_map_) {
        map_ *= alpha;
      } else {
        v_ *= alpha;
      }
    }
    
    int numRows() const { 
      return using_map_ ? map_.size() : v_.size(); 
    }
    
    int stride() const { 
      return numRows(); // For LAPACK compatibility - return size for vectors
    }
};

template<typename Real>
class Matrix{
  private:
    EMatrix<Real> M_;  // Owned data
    Eigen::Map<EMatrix<Real>> map_;  // View of external data
    bool using_map_;  // Flag to track which member to use
    
  public:
    Matrix() : M_(), map_(nullptr, 0, 0), using_map_(false) {}
    
    Matrix(const Matrix& other) : M_(), map_(nullptr, 0, 0), using_map_(false) {
      if (other.using_map_) {
        // Copy data from the map into owned storage
        M_ = other.map_;
        map_ = Eigen::Map<EMatrix<Real>>(M_.data(), M_.rows(), M_.cols());
      } else {
        M_ = other.M_;  
        map_ = Eigen::Map<EMatrix<Real>>(M_.data(), M_.rows(), M_.cols());
      }
    }
    
    Matrix(int rows, int columns) : M_(rows, columns), map_(M_.data(), rows, columns), using_map_(false) {
      // Create owned matrix and map pointing to it
    }
    
    // New constructor to support DataAccess with external data pointer
    Matrix(DataAccess access, Real* data, int rows, int cols, int stride = 0) 
      : M_(), map_(nullptr, 0, 0), using_map_(true) {
      if (stride == 0) stride = rows;  // Default to column-major
      
      if (access == View) {
        // Create map to external data
        map_ = Eigen::Map<EMatrix<Real>>(data, rows, cols);
      } else { // Copy
        // Create owned storage and copy data
        M_.resize(rows, cols);
        using_map_ = false;
        map_ = Eigen::Map<EMatrix<Real>>(M_.data(), rows, cols);
        
        // Copy data respecting stride (assume column-major external data)
        for (int j = 0; j < cols; j++) {
          for (int i = 0; i < rows; i++) {
            M_(i, j) = data[i + j * stride];
          }
        }
      }
    }
    
    // Existing block constructor - modified to work with dual member approach
    Matrix(DataAccess access, Matrix<Real> A, int rows, int cols, int rowstart = 0, int colstart=0)
      : M_(), map_(nullptr, 0, 0), using_map_(false) {
      if(access == Copy) {
        if (A.using_map_) {
          M_ = A.map_.block(rowstart, colstart, rows, cols).eval();
        } else {
          M_ = A.M_.block(rowstart, colstart, rows, cols).eval();
        }
        map_ = Eigen::Map<EMatrix<Real>>(M_.data(), M_.rows(), M_.cols());
      } else {
        // View case - this is tricky with Map, we'll create a copy for now
        // A proper implementation would need nested Map support
        if (A.using_map_) {
          M_ = A.map_.block(rowstart, colstart, rows, cols).eval();
        } else {
          M_ = A.M_.block(rowstart, colstart, rows, cols).eval();
        }
        map_ = Eigen::Map<EMatrix<Real>>(M_.data(), M_.rows(), M_.cols());
      }
    }
    
    Real* values() { 
      return using_map_ ? const_cast<Real*>(map_.data()) : const_cast<Real*>(M_.data()); 
    }
    
    int stride() { 
      // Return leading dimension for LAPACK compatibility
      return using_map_ ? map_.rows() : M_.rows();
    }
    void reshape(int m, int n) { 
      if (!using_map_) {
        M_.resize(m, n); 
        map_ = Eigen::Map<EMatrix<Real>>(M_.data(), m, n);
      }
      // Cannot reshape a mapped matrix
    }
    
    Real normOne() { 
      return using_map_ ? map_.template lpNorm<1>() : M_.template lpNorm<1>(); 
    }
    
    Eigen::PartialPivLU<EMatrix<Real>> partialPivLu(bool inplace) {
      if (using_map_) {
        if(inplace)
          return Eigen::PartialPivLU<Eigen::Ref<Eigen::Map<EMatrix<Real>>>>(map_);
        else
          return map_.partialPivLu();
      } else {
        if(inplace)
          return Eigen::PartialPivLU<Eigen::Ref<EMatrix<Real>>>(M_);
        else
          return M_.partialPivLu();
      }
    }
    
    Real& operator()(int i, int j) {
      return using_map_ ? map_(i, j) : M_(i, j);
    }

    void multiply(ETransp transa, ETransp transb, Real alpha, const Matrix& A, const Matrix& B, Real beta) {
      // Get the appropriate matrix references
      auto& A_ref = A.using_map_ ? A.map_ : A.M_;
      auto& B_ref = B.using_map_ ? B.map_ : B.M_;
//      auto& this_ref = using_map_ ? map_ : M_;

      EMatrix<Real> AA;
      if(transa == NO_TRANS)
        AA = A_ref;
      else if(transa == TRANS)
        AA = A_ref.transpose();
      else
        AA = A_ref.conjugate();
        
      EMatrix<Real> BB;
      if(transb == NO_TRANS)
        BB = B_ref;
      else if(transb == TRANS)
        BB = B_ref.transpose();
      else
        BB = B_ref.conjugate();
        
      if(beta != Real(1)) map_ *= beta;
      map_.noalias() += alpha * AA * BB;
    }

    int numRows() const { 
      return using_map_ ? map_.rows() : M_.rows(); 
    }
    
    int numCols() const { 
      return using_map_ ? map_.cols() : M_.cols(); 
    }
    
    // Set all elements to scalar value
    void putScalar(Real value) {
      if (using_map_) {
        map_.setConstant(value);
      } else {
        M_.setConstant(value);
      }
    }

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
