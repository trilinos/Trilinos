// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_GMRES_H
#define ROL_GMRES_H

/** \class ROL::GMRES
    \brief Preconditioned GMRES solver.
*/

#include "ROL_Krylov.hpp"
#include "ROL_Types.hpp"
#include "ROL_LAPACK.hpp"
#include "ROL_LinearAlgebra.hpp"


namespace ROL {

template<class Real>
class VectorArray {
private:
  std::vector<Ptr<Vector<Real>>> data_;
  unsigned size_;

public:
  VectorArray() : size_(0u) {}

  void allocate(const Vector<Real> &x, unsigned numVec) {
    data_.clear();
    for (unsigned i = 0u; i < numVec; ++i)
      data_.push_back(x.clone());
    size_ = numVec;
  }

  const Ptr<Vector<Real>> get(const Vector<Real> &x, unsigned ind) {
    if (ind < size_) {
      return data_[ind];
    }
    else if (ind == size_) {
      data_.push_back(x.clone());
      size_++;
      return data_[ind];
    }
    else {
      throw Exception::NotImplemented(">>> GMRES : Index out of range!");
      return nullPtr;
    }
  }

  const Ptr<Vector<Real>> get(unsigned ind) {
    if (size_ > 0u) {
      return get(*data_[0],ind);
    }
    else {
      throw Exception::NotImplemented(">>> GMRES : No vectors allocated!");
      return nullPtr;
    }
  }
};

template<class Real>
class GMRES : public Krylov<Real> {

  typedef LA::Matrix<Real> SDMatrix;
  typedef LA::Vector<Real> SDVector;

private:

  Ptr<Vector<Real>> r_;
  Ptr<Vector<Real>> z_;
  Ptr<Vector<Real>> w_;

  Ptr<SDMatrix> H_;      // quasi-Hessenberg matrix
  Ptr<SDVector> cs_;     // Givens Rotations cosine components
  Ptr<SDVector> sn_;     // Givens Rotations sine components
  Ptr<SDVector> s_;
  Ptr<SDVector> y_;
  Ptr<SDVector> cnorm_;

  Ptr<std::vector<Real>> res_;
  Ptr<VectorArray<Real>> V_, Z_;

  bool isInitialized_;
  bool useInexact_;
  bool useInitialGuess_;    // If false, inital x will be ignored and zero vec used
  bool printIters_;
  Ptr<std::ostream> outStream_;

  LAPACK<int,Real> lapack_;

public:

  GMRES( ParameterList &parlist ) : Krylov<Real>(parlist), isInitialized_(false), printIters_(false) {

    using std::vector;

    Real zero(0);

    ROL::ParameterList &gList = parlist.sublist("General");
    ROL::ParameterList &kList = gList.sublist("Krylov");

    useInexact_      = gList.get("Inexact Hessian-Times-A-Vector",false);
    useInitialGuess_ = kList.get("Use Initial Guess",false);
    int maxit = Krylov<Real>::getMaximumIteration();

    H_     = makePtr<SDMatrix>( maxit+1, maxit );
    cs_    = makePtr<SDVector>( maxit );
    sn_    = makePtr<SDVector>( maxit );
    s_     = makePtr<SDVector>( maxit+1 );
    y_     = makePtr<SDVector>( maxit+1 );
    cnorm_ = makePtr<SDVector>( maxit );
    res_   = makePtr<std::vector<Real>>(maxit+1,zero);

    V_ = makePtr<VectorArray<Real>>();
    Z_ = makePtr<VectorArray<Real>>();
  }

  Real run( Vector<Real> &x, LinearOperator<Real> &A, const Vector<Real> &b,
            LinearOperator<Real> &M, int &iter, int &flag ) {

    Real absTol = Krylov<Real>::getAbsoluteTolerance();
    Real relTol = Krylov<Real>::getRelativeTolerance();
    int maxit = Krylov<Real>::getMaximumIteration();

    flag = 0;

    Real zero(0), one(1);

    if ( !isInitialized_ ) {
      r_  = b.clone();
      w_  = b.clone();
      z_  = x.clone();

      V_->allocate(b,std::min(20u,static_cast<unsigned>(maxit))+1u);
      Z_->allocate(x,std::min(20u,static_cast<unsigned>(maxit)));

      isInitialized_ = true;
    }

    Real itol  = std::sqrt(ROL_EPSILON<Real>());

    // Compute initial residual
    if(useInitialGuess_) {

      A.apply(*r_,x,itol);
      r_->scale(-one);
      r_->plus(b);       // r = b-Ax

    }
    else {
      x.zero();
      r_->set(b);
    }

    Real temp  = 0;

    (*res_)[0] = r_->norm();

    if (printIters_)
      *outStream_ << "GMRES Iteration " << 0 << ", Residual = " << (*res_)[0] << "\n";

    // This should be a tolerance check
    Real rtol = std::min(absTol,relTol*(*res_)[0]);
    if ((*res_)[0] <= rtol) {
      iter = 0;
      flag = 0;
      return (*res_)[0];
    }

    V_->get(0)->set(*r_);
    V_->get(0)->scale(one/(*res_)[0]);

    (*s_)(0) = (*res_)[0];

    for( iter=0; iter<maxit; ++iter ) {

      if( useInexact_ ) itol = rtol/(maxit*(*res_)[iter]);

      // Apply right preconditioner
      M.applyInverse(*(Z_->get(iter)),*(V_->get(iter)),itol);

      // Apply operator
      A.apply(*w_,*(Z_->get(iter)),itol);

      // Evaluate coefficients and orthogonalize using Gram-Schmidt
      for( int k=0; k<=iter; ++k ) {
        (*H_)(k,iter) = w_->dot(*(V_->get(k)));
        w_->axpy( -(*H_)(k,iter), *(V_->get(k)) );
      }

      (*H_)(iter+1,iter) = w_->norm();

      (V_->get(iter+1))->set(*w_);
      (V_->get(iter+1))->scale(one/((*H_)(iter+1,iter)));

      // Apply Givens rotations
      for( int k=0; k<=iter-1; ++k ) {
        temp            =  (*cs_)(k)*(*H_)(k,iter) + (*sn_)(k)*(*H_)(k+1,iter);
        (*H_)(k+1,iter) = -(*sn_)(k)*(*H_)(k,iter) + (*cs_)(k)*(*H_)(k+1,iter);
        (*H_)(k,iter)   = temp;
      }

      // Form i-th rotation matrix
      if( (*H_)(iter+1,iter) == zero ) {
        (*cs_)(iter) = one;
        (*sn_)(iter) = zero;
      }
      else if ( std::abs((*H_)(iter+1,iter)) > std::abs((*H_)(iter,iter)) ) {
        temp = (*H_)(iter,iter) / (*H_)(iter+1,iter);
        (*sn_)(iter) = one / std::sqrt( one + temp*temp );
        (*cs_)(iter) = temp*(*sn_)(iter);
      }
      else {
        temp = (*H_)(iter+1,iter) / (*H_)(iter,iter);
        (*cs_)(iter) = one / std::sqrt( one + temp*temp );
        (*sn_)(iter) = temp*(*cs_)(iter);
      }

      // Approximate residual norm
      temp               = (*cs_)(iter)*(*s_)(iter);
      (*s_)(iter+1)      = -(*sn_)(iter)*(*s_)(iter);
      (*s_)(iter)        = temp;
      (*H_)(iter,iter)   = (*cs_)(iter)*(*H_)(iter,iter) + (*sn_)(iter)*(*H_)(iter+1,iter);
      (*H_)(iter+1,iter) = zero;
      (*res_)[iter+1]    = std::abs((*s_)(iter+1));

      if (printIters_) {
        *outStream_ << "GMRES Iteration " << iter+1 << ", Residual = " << (*res_)[iter+1] << "\n";
      }

      // Update solution approximation.
      const char uplo = 'U';
      const char trans = 'N';
      const char diag = 'N';
      const char normin = 'N';
      Real scaling = zero;
      int info = 0;
      *y_ = *s_;
      lapack_.LATRS(uplo, trans, diag, normin, iter+1, H_->values(), maxit+1, y_->values(), &scaling, cnorm_->values(), &info);

      z_->zero();

      for( int k=0; k<=iter;++k ) {
        z_->axpy((*y_)(k),*(Z_->get(k)));
      }

      if( (*res_)[iter+1] <= rtol ) {
        // Update solution vector
        x.plus(*z_);
        break;
      }

    } // loop over iter

    if(iter == maxit) {
      flag = 1;
      x.plus(*z_);
      return (*res_)[iter];
    }

    return (*res_)[iter+1];
  }

  void enableOutput(std::ostream & outStream)  {
    printIters_ = true;
    outStream_ = ROL::makePtrFromRef(outStream);;
  }

  void disableOutput() {printIters_ = false;}

}; // class GMRES

} // namespace ROL

#endif // ROL_GMRES_H

