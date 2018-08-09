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

#ifndef ROL_SKETCH_H
#define ROL_SKETCH_H

#include "ROL_Vector.hpp"
#include "ROL_LinearAlgebra.hpp"
#include "ROL_LAPACK.hpp"
#include "ROL_Types.hpp"

/** @ingroup func_group
    \class ROL::Sketch
    \brief Provides an interface for randomized sketching.

    ---
*/

namespace ROL {

template <class Real>
class Sketch {
private:
  std::vector<Ptr<Vector<Real>>> Psi_, Y_;

  LA::Matrix<Real> Omega_, W_;

  int ncol_, rank_, k_, l_;

  ROL::LAPACK<int,Real> lapack_;
  int solverType_;

  bool flagQ_, flagX_;

  void mgs(std::vector<Ptr<Vector<Real>>> &Y) const {
    const Real one(1);
    Real rii(0), rij(0);
    for (int i = 0; i < k_; ++i) {
      rii = Y[i]->norm();
      if (rii > ROL_EPSILON<Real>()) { // Ignore update if Y[i] is zero.
        Y[i]->scale(one/rii);
        for (int j = i+1; j < k_; ++j) {
          rij = Y[i]->dot(*Y[j]);
          Y[j]->axpy(-rij,*Y[i]);
        }
      }
    }
  }

  void computeQ(void) {
    if (!flagQ_) {
      mgs(Y_);
      mgs(Y_);
      flagQ_ = true;
    }
  }

  void computeX(void) {
    computeQ();
    if (!flagX_) {
      LA::Matrix<Real> Z(l_,k_);
      for (int i = 0; i < l_; ++i) {
        for (int j = 0; j < k_; ++j) {
          Z(i,j) = Psi_[i]->dot(*Y_[j]);
        }
      }
      // Solve least squares problem using LAPACK
      int M      = l_;
      int N      = k_;
      int NRHS   = ncol_;
      int LDA    = M;
      int LDB    = M;
      std::vector<Real> WORK(1);
      int LWORK  = -1;
      int INFO;
      if (solverType_ == 0 ) { // QR
        char TRANS = 'N';
        lapack_.GELS(TRANS,M,N,NRHS,Z.values(),LDA,W_.values(),LDB,&WORK[0],LWORK,&INFO);
        LWORK = static_cast<int>(WORK[0]);
        WORK.resize(LWORK);
        lapack_.GELS(TRANS,M,N,NRHS,Z.values(),LDA,W_.values(),LDB,&WORK[0],LWORK,&INFO);
      }
      else {                   // SVD
        std::vector<Real> S(N);
        Real RCOND = -1.0;
        int RANK;
        lapack_.GELSS(M,N,NRHS,Z.values(),LDA,W_.values(),LDB,&S[0],RCOND,&RANK,&WORK[0],LWORK,&INFO);
        LWORK = static_cast<int>(WORK[0]);
        WORK.resize(LWORK);
        lapack_.GELSS(M,N,NRHS,Z.values(),LDA,W_.values(),LDB,&S[0],RCOND,&RANK,&WORK[0],LWORK,&INFO);
      }
      // Set flag
      flagX_ = true;
    }
  }

  bool testQ(std::ostream &outStream = std::cout, const int verbosity = 0) {
    const Real one(1), tol(std::sqrt(ROL_EPSILON<Real>()));
    computeQ();
    Real qij(0), err(0), maxerr(0);
    std::ios_base::fmtflags oflags(outStream.flags());
    if (verbosity > 0) {
      outStream << std::scientific << std::setprecision(3);
    }
    if (verbosity > 1) {
      outStream << std::endl
                << " Printing Q'Q...This should be approximately equal to I"
                << std::endl << std::endl;
    }
    for (int i = 0; i < k_; ++i) {
      for (int j = 0; j < k_; ++j) {
        qij    = Y_[i]->dot(*Y_[j]);
        err    = (i==j ? std::abs(qij-one) : std::abs(qij));
        maxerr = (err > maxerr ? err : maxerr);
        if (verbosity > 1) {
          outStream << std::setw(12) << std::right << qij;
        }
      }
      if (verbosity > 1) {
        outStream << std::endl;
      }
    }
    if (verbosity > 0) {
      outStream << std::endl << " TEST ORTHOGONALIZATION: Max Error = "
                << std::setw(12) << std::right << maxerr
                << std::endl;
      outStream.flags(oflags);
    }
    return (maxerr < tol ? true : false);
  }

  void reset(void) {
    const Real one(1);
    // Randomize Psi and Omega and zero W and Y
    for (int i = 0; i < l_; ++i) {
      Psi_[i]->randomize(-one,one);
      for (int j = 0; j < ncol_; ++j) {
        W_(i,j) = static_cast<Real>(0);
      }
    }
    Real a(2), b(-1), x(0);
    for (int i = 0; i < k_; ++i) {
      Y_[i]->zero();
      for (int j = 0; j < ncol_; ++j) {
        x = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
        Omega_(j,i) = a*x + b;
      }
    }
  }

public:
  virtual ~Sketch(void) {}

  Sketch(const Vector<Real> &x, const int ncol, const int rank, const int solverType=0)
    : ncol_(ncol), rank_(rank), solverType_(solverType), flagQ_(false), flagX_(false) {
    // Compute reduced dimensions
    l_ = std::min(4*rank_+3,ncol_);
    k_ = std::min(2*rank_+1,l_);
    // Initialize matrix storage
    Psi_.clear(); Y_.clear();
    Omega_.reshape(ncol_,k_); W_.reshape(l_,ncol_);
    for (int i = 0; i < l_; ++i) {
      Psi_.push_back(x.clone());
    }
    for (int i = 0; i < k_; ++i) {
      Y_.push_back(x.clone());
    }
    // Randomize Psi and Omega and zero W and Y
    reset();
  }

  void setRank(const int rank) {
    rank_ = rank;
    // Compute reduced dimensions
    l_ = std::min(4*rank_+3,ncol_);
    k_ = std::min(2*rank_+1,l_);
    update();
  }

  void update(void) {
    flagQ_ = false;
    flagX_ = false;
    // Randomize Psi and Omega and zero W and Y
    reset();
  }

  void advance(const Real alpha, Vector<Real> &h, const int col, const Real beta = 1.0) {
    // Check to see if col is less than ncol_
    if ( col >= ncol_ ) {
      throw Exception::NotImplemented(">>> ROL::Sketch::advance: Input column index exceeds total number of columns!");
    }
    if (!flagQ_ && !flagX_) {
      // Update Y
      for (int i = 0; i < k_; ++i) {
        Y_[i]->scale(beta);
        Y_[i]->axpy(alpha*Omega_(col,i),h);
      }
      // Update W
      for (int i = 0; i < l_; ++i) {
        for (int j = 0; j < ncol_; ++j) {
          W_(i,j) *= beta;
        }
        W_(i,col) += alpha*h.dot(*Psi_[i]);
      }
    }
    else {
      throw Exception::NotImplemented(">>> ROL::Sketch::advance: Reconstruct has already been called!");
    }
  }

  void reconstruct(Vector<Real> &a, const int col) {
    // Check to see if col is less than ncol_
    if ( col >= ncol_ ) {
      throw Exception::NotImplemented(">>> ROL::Sketch::reconstruct: Input column index exceeds total number of columns!");
    }
    // Compute QR factorization of Y store in Y
    computeQ();
    // Solve (Psi Q)X = W (over determined least squares) store in W
    computeX();
    // Apply Q to col column of X
    a.zero();
    for (int i = 0; i < k_; ++i) {
      a.axpy(W_(i,col),*Y_[i]);
    }
  }

  bool test(const int rank, std::ostream &outStream = std::cout, const int verbosity = 0) {
    const Real one(1), tol(std::sqrt(ROL_EPSILON<Real>()));
    // Initialize low rank factors
    std::vector<Ptr<Vector<Real>>> U(rank);
    LA::Matrix<Real> V(ncol_,rank);
    for (int i = 0; i < rank; ++i) {
      U[i] = Y_[0]->clone();
      U[i]->randomize(-one,one);
      for (int j = 0; j < ncol_; ++j) {
        V(j,i) = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
      }
    }
    // Initialize A and build sketch
    update();
    std::vector<Ptr<Vector<Real>>> A(ncol_);
    for (int i = 0; i < ncol_; ++i) {
      A[i] = Y_[0]->clone(); A[i]->zero();
      for (int j = 0; j < rank; ++j) {
        A[i]->axpy(V(i,j),*U[j]);
      }
      advance(one,*A[i],i,one);
    }
    // Test QR decomposition of Y
    bool flagQ = testQ(outStream, verbosity);
    // Test reconstruction of A
    Real nerr(0), maxerr(0);
    Ptr<Vector<Real>> err = Y_[0]->clone();
    for (int i = 0; i < ncol_; ++i) {
      reconstruct(*err,i);
      err->axpy(-one,*A[i]);
      nerr = err->norm();
      maxerr = (nerr > maxerr ? nerr : maxerr);
    }
    if (verbosity > 0) {
      std::ios_base::fmtflags oflags(outStream.flags());
      outStream << std::scientific << std::setprecision(3);
      outStream << " TEST RECONSTRUCTION:    Max Error = "
                << std::setw(12) << std::right << maxerr
                << std::endl << std::endl;
      outStream.flags(oflags);
    }
    return flagQ & (maxerr < tol ? true : false);
  }

}; // class Sketch

} // namespace ROL

#endif
