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
#include <random>
#include <chrono>

/** @ingroup func_group
    \class ROL::Sketch
    \brief Provides an interface for randomized sketching.

    ---
*/

namespace ROL {

template <class Real>
class Sketch {
private:
  std::vector<Ptr<Vector<Real>>> Upsilon_, Phi_, Y_;
  LA::Matrix<Real> Omega_, Psi_, X_, Z_, C_;

  int maxRank_, ncol_, rank_, k_, s_;

  const Real orthTol_;
  const int  orthIt_;

  const bool truncate_;

  LAPACK<int,Real> lapack_;

  bool flagP_, flagQ_, flagC_;

  Ptr<std::ostream> out_;

  Ptr<Elementwise::NormalRandom<Real>> nrand_;
  Ptr<std::mt19937_64> gen_;
  Ptr<std::normal_distribution<Real>> dist_;

  int computeP(void) {
    int INFO(0);
    if (!flagP_) {
      // Solve least squares problem using LAPACK
      int M      = ncol_;
      int N      = k_;
      int K      = std::min(M,N);
      int LDA    = M;
      std::vector<Real> TAU(K);
      std::vector<Real> WORK(1);
      int LWORK  = -1;
      // Compute QR factorization of X
      lapack_.GEQRF(M,N,X_.values(),LDA,&TAU[0],&WORK[0],LWORK,&INFO);
      LWORK = static_cast<int>(WORK[0]);
      WORK.resize(LWORK);
      lapack_.GEQRF(M,N,X_.values(),LDA,&TAU[0],&WORK[0],LWORK,&INFO);
      // Generate Q
      LWORK = -1;
      lapack_.ORGQR(M,N,K,X_.values(),LDA,&TAU[0],&WORK[0],LWORK,&INFO);
      LWORK = static_cast<int>(WORK[0]);
      WORK.resize(LWORK);
      lapack_.ORGQR(M,N,K,X_.values(),LDA,&TAU[0],&WORK[0],LWORK,&INFO);
      flagP_ = true;
    }
    return INFO;
  }

  void mgs2(std::vector<Ptr<Vector<Real>>> &Y) const {
    const Real one(1);
    Real rjj(0), rij(0);
    std::vector<Real> normQ(k_,0);
    bool flag(true);
    for (int j = 0; j < k_; ++j) {
      rjj = Y[j]->norm();
      if (rjj > ROL_EPSILON<Real>()) { // Ignore update if Y[i] is zero.
        for (int k = 0; k < orthIt_; ++k) {
          for (int i = 0; i < j; ++i) {
            rij = Y[i]->dot(*Y[j]);
            Y[j]->axpy(-rij,*Y[i]);
          }
          normQ[j] = Y[j]->norm();
          flag = true;
          for (int i = 0; i < j; ++i) {
            rij = std::abs(Y[i]->dot(*Y[j]));
            if (rij > orthTol_*normQ[j]*normQ[i]) {
              flag = false;
              break;
            }
          }
          if (flag) {
            break;
          }
        }
      }
      rjj = normQ[j];
      Y[j]->scale(one/rjj);
    }
  }

  int computeQ(void) {
    if (!flagQ_) {
      mgs2(Y_);
      flagQ_ = true;
    }
    return 0;
  }

  int LSsolver(LA::Matrix<Real> &A, LA::Matrix<Real> &B, const bool trans = false) const {
    int flag(0);
    char TRANS = (trans ? 'T' : 'N');
    int M      = A.numRows();
    int N      = A.numCols();
    int NRHS   = B.numCols();
    int LDA    = M;
    int LDB    = std::max(M,N);
    std::vector<Real> WORK(1);
    int LWORK  = -1;
    int INFO;
    lapack_.GELS(TRANS,M,N,NRHS,A.values(),LDA,B.values(),LDB,&WORK[0],LWORK,&INFO);
    flag += INFO;
    LWORK = static_cast<int>(WORK[0]);
    WORK.resize(LWORK);
    lapack_.GELS(TRANS,M,N,NRHS,A.values(),LDA,B.values(),LDB,&WORK[0],LWORK,&INFO);
    flag += INFO;
    return flag;
  }

  int lowRankApprox(LA::Matrix<Real> &A, const int r) const {
    const Real zero(0);
    char JOBU  = 'S';
    char JOBVT = 'S';
    int  M     = A.numRows();
    int  N     = A.numCols();
    int  K     = std::min(M,N);
    int  LDA   = M;
    std::vector<Real> S(K);
    LA::Matrix<Real> U(M,K);
    int  LDU   = M;
    LA::Matrix<Real> VT(K,N);
    int  LDVT  = K;
    std::vector<Real> WORK(1), WORK0(1);
    int  LWORK = -1;
    int  INFO;
    lapack_.GESVD(JOBU,JOBVT,M,N,A.values(),LDA,&S[0],U.values(),LDU,VT.values(),LDVT,&WORK[0],LWORK,&WORK0[0],&INFO);
    LWORK = static_cast<int>(WORK[0]);
    WORK.resize(LWORK);
    lapack_.GESVD(JOBU,JOBVT,M,N,A.values(),LDA,&S[0],U.values(),LDU,VT.values(),LDVT,&WORK[0],LWORK,&WORK0[0],&INFO);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        A(i,j) = zero;
        for (int k = 0; k < r; ++k) {
          A(i,j) += S[k] * U(i,k) * VT(k,j);
        }
      }
    }
    return INFO;
  }

  int computeC(void) {
    int infoP(0), infoQ(0), infoLS1(0), infoLS2(0), infoLRA(0);
    infoP = computeP();
    infoQ = computeQ();
    if (!flagC_) {
      const Real zero(0);
      LA::Matrix<Real> L(s_,k_), R(s_,k_);
      for (int i = 0; i < s_; ++i) {
        for (int j = 0; j < k_; ++j) {
          L(i,j)  = Phi_[i]->dot(*Y_[j]);
          R(i,j)  = zero;
          for (int k = 0; k < ncol_; ++k) {
            R(i,j) += Psi_(k,i) * X_(k,j);
          }
        }
      }
      // Solve least squares problems using LAPACK
      infoLS1 = LSsolver(L,Z_,false);
      LA::Matrix<Real> Z(s_,k_);
      for (int i = 0; i < k_; ++i) {
        for (int j = 0; j < s_; ++j) {
          Z(j,i) = Z_(i,j);
        }
      }
      infoLS2 = LSsolver(R,Z,false);
      for (int i = 0; i < k_; ++i) {
        for (int j = 0; j < k_; ++j) {
          C_(j,i) = Z(i,j);
        }
      }
      // Compute best rank r approximation
      if (truncate_) {
        infoLRA = lowRankApprox(C_,rank_);
      }
      // Set flag
      flagC_ = true;
    }
    return std::abs(infoP)+std::abs(infoQ)+std::abs(infoLS1)
                          +std::abs(infoLS2)+std::abs(infoLRA);
  }

  void reset(void) {
    const Real zero(0);
    // Randomize Upsilon, Omega, Psi and Phi, and zero X, Y and Z
    X_.scale(zero); Z_.scale(zero); C_.scale(zero);
    for (int i = 0; i < s_; ++i) {
      Phi_[i]->applyUnary(*nrand_);
      for (int j = 0; j < ncol_; ++j) {
        Psi_(j,i) = (*dist_)(*gen_);
      }
    }
    for (int i = 0; i < k_; ++i) {
      Y_[i]->zero();
      Upsilon_[i]->applyUnary(*nrand_);
      for (int j = 0; j < ncol_; ++j) {
        Omega_(j,i) = (*dist_)(*gen_);
      }
    }
  }

//  void reset(void) {
//    const Real zero(0), one(1), a(2), b(-1);
//    Real x(0);
//    // Randomize Upsilon, Omega, Psi and Phi, and zero X, Y and Z
//    X_.scale(zero); Z_.scale(zero); C_.scale(zero);
//    for (int i = 0; i < s_; ++i) {
//      Phi_[i]->randomize(-one,one);
//      for (int j = 0; j < ncol_; ++j) {
//        x = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
//        Psi_(j,i) = a*x + b; 
//      }
//    }
//    for (int i = 0; i < k_; ++i) {
//      Y_[i]->zero();
//      Upsilon_[i]->randomize(-one,one);
//      for (int j = 0; j < ncol_; ++j) {
//        x = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
//        Omega_(j,i) = a*x + b;
//      }
//    }
//  }

public:
  virtual ~Sketch(void) {}

  Sketch(const Vector<Real> &x, const int ncol, const int rank,
         const Real orthTol = 1e-8, const int orthIt = 2,
         const bool truncate = false,
         const unsigned dom_seed = 0, const unsigned rng_seed = 0)
    : ncol_(ncol), orthTol_(orthTol), orthIt_(orthIt),
      truncate_(truncate), flagP_(false), flagQ_(false), flagC_(false),
      out_(nullPtr) {
    Real mu(0), sig(1);
    nrand_ = makePtr<Elementwise::NormalRandom<Real>>(mu,sig,dom_seed);
    unsigned seed = rng_seed;
    if (seed == 0) {
      seed = std::chrono::system_clock::now().time_since_epoch().count();
    }
    gen_  = makePtr<std::mt19937_64>(seed);
    dist_ = makePtr<std::normal_distribution<Real>>(mu,sig);
    // Compute reduced dimensions
    maxRank_ = std::min(ncol_, x.dimension());
    rank_ = std::min(rank, maxRank_);
    k_    = std::min(2*rank_+1, maxRank_);
    s_    = std::min(2*k_   +1, maxRank_);
    // Initialize matrix storage
    Upsilon_.clear(); Phi_.clear(); Y_.clear();
    Omega_.reshape(ncol_,k_); Psi_.reshape(ncol_,s_);
    X_.reshape(ncol_,k_); Z_.reshape(s_,s_); C_.reshape(k_,k_);
    for (int i = 0; i < k_; ++i) {
      Y_.push_back(x.clone());
      Upsilon_.push_back(x.clone());
    }
    for (int i = 0; i < s_; ++i) {
      Phi_.push_back(x.clone());
    }
    // Randomize Psi and Omega and zero W and Y
    reset();
  }

  void setStream(Ptr<std::ostream> &out) {
    out_ = out;
  }

  void setRank(const int rank) {
    rank_ = std::min(rank, maxRank_);
    // Compute reduced dimensions
    int sold = s_, kold = k_;
    k_ = std::min(2*rank_+1, maxRank_);
    s_ = std::min(2*k_   +1, maxRank_);
    Omega_.reshape(ncol_,k_); Psi_.reshape(ncol_,s_);
    X_.reshape(ncol_,k_); Z_.reshape(s_,s_); C_.reshape(k_,k_);
    if (s_ > sold) {
      for (int i = sold; i < s_; ++i) {
        Phi_.push_back(Phi_[0]->clone());
      }
    }
    if (k_ > kold) {
      for (int i = kold; i < k_; ++i) {
        Y_.push_back(Y_[0]->clone());
        Upsilon_.push_back(Upsilon_[0]->clone());
      }
    }
    update();
    if ( out_ != nullPtr ) {
      *out_ << std::string(80,'=')            << std::endl;
      *out_ << "  ROL::Sketch::setRank"       << std::endl;
      *out_ << "    **** Rank    = " << rank_ << std::endl;
      *out_ << "    **** k       = " << k_    << std::endl;
      *out_ << "    **** s       = " << s_    << std::endl;
      *out_ << std::string(80,'=')            << std::endl;
    }
  }

  void update(void) {
    flagP_ = false;
    flagQ_ = false;
    flagC_ = false;
    // Randomize Psi and Omega and zero W and Y
    reset();
  }

  int advance(const Real nu, Vector<Real> &h, const int col, const Real eta = 1.0) {
    // Check to see if col is less than ncol_
    if ( col >= ncol_ || col < 0 ) {
      // Input column index  out of range!
      return 1;
    }
    if (!flagP_ && !flagQ_ && !flagC_) {
      for (int i = 0; i < k_; ++i) {
        // Update X
        for (int j = 0; j < ncol_; ++j) {
          X_(j,i) *= eta;
        }
        X_(col,i) += nu*h.dot(*Upsilon_[i]);
        // Update Y
        Y_[i]->scale(eta);
        Y_[i]->axpy(nu*Omega_(col,i),h);
      }
      // Update Z
      Real hphi(0);
      for (int i = 0; i < s_; ++i) {
        hphi = h.dot(*Phi_[i]);
        for (int j = 0; j < s_; ++j) {
          Z_(i,j) *= eta;
          Z_(i,j) += nu*Psi_(col,j)*hphi;
        }
      }
      if ( out_ != nullPtr ) {
        *out_ << std::string(80,'=')               << std::endl;
        *out_ << "  ROL::Sketch::advance"          << std::endl;
        *out_ << "    **** col     = " << col      << std::endl;
        *out_ << "    **** norm(h) = " << h.norm() << std::endl;
        *out_ << std::string(80,'=')               << std::endl;
      }
    }
    else {
      // Reconstruct has already been called!
      return 1;
    }
    return 0;
  }

  int reconstruct(Vector<Real> &a, const int col) {
    // Check to see if col is less than ncol_
    if ( col >= ncol_ || col < 0 ) {
      // Input column index out of range!
      return 2;
    }
    const Real zero(0);
    int flag(0);
    // Compute QR factorization of X store in X
    flag = computeP();
    if (flag > 0 ) {
      return 3;
    }
    // Compute QR factorization of Y store in Y
    flag = computeQ();
    if (flag > 0 ) {
      return 4;
    }
    // Compute (Phi Q)\Z/(Psi P)* store in C
    flag = computeC();
    if (flag > 0 ) {
      return 5;
    }
    // Recover sketch
    a.zero();
    Real coeff(0);
    for (int i = 0; i < k_; ++i) {
      coeff = zero;
      for (int j = 0; j < k_; ++j) {
        coeff += C_(i,j) * X_(col,j);
      }
      a.axpy(coeff,*Y_[i]);
    }
    if ( out_ != nullPtr ) {
      *out_ << std::string(80,'=')               << std::endl;
      *out_ << "  ROL::Sketch::reconstruct"      << std::endl;
      *out_ << "    **** col     = " << col      << std::endl;
      *out_ << "    **** norm(a) = " << a.norm() << std::endl;
      *out_ << std::string(80,'=')               << std::endl;
    }
    return 0;
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
    // Test QR decomposition of X
    bool flagP = testP(outStream, verbosity);
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
      outStream << std::scientific << std::setprecision(3) << std::endl;
      outStream << " TEST RECONSTRUCTION:    Max Error = "
                << std::setw(12) << std::right << maxerr
                << std::endl << std::endl;
      outStream.flags(oflags);
    }
    return flagP & flagQ & (maxerr < tol ? true : false);
  }

private:

  // Test functions
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
      if (maxerr > tol) {
        break;
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

  bool testP(std::ostream &outStream = std::cout, const int verbosity = 0) {
    const Real zero(0), one(1), tol(std::sqrt(ROL_EPSILON<Real>()));
    computeP();
    Real qij(0), err(0), maxerr(0);
    std::ios_base::fmtflags oflags(outStream.flags());
    if (verbosity > 0) {
      outStream << std::scientific << std::setprecision(3);
    }
    if (verbosity > 1) {
      outStream << std::endl
                << " Printing P'P...This should be approximately equal to I"
                << std::endl << std::endl;
    }
    for (int i = 0; i < k_; ++i) {
      for (int j = 0; j < k_; ++j) {
        qij    = zero;
        for (int k = 0; k < ncol_; ++k) {
          qij += X_(k,i) * X_(k,j);
        }
        err    = (i==j ? std::abs(qij-one) : std::abs(qij));
        maxerr = (err > maxerr ? err : maxerr);
        if (verbosity > 1) {
          outStream << std::setw(12) << std::right << qij;
        }
      }
      if (verbosity > 1) {
        outStream << std::endl;
      }
      if (maxerr > tol) {
        break;
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

}; // class Sketch

} // namespace ROL

#endif
