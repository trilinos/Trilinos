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

#ifndef ROL_TEUCHOS_LINEAROPERATOR_H
#define ROL_TEUCHOS_LINEAROPERATOR_H

#include "ROL_Ptr.hpp"
#include "ROL_TeuchosVector.hpp"
#include "ROL_LinearOperator.hpp"
#include "ROL_LAPACK.hpp"

#include "Teuchos_SerialDenseVector.hpp"

/** @ingroup func_group
    \class ROL::TeuchosLinearOperator
    \brief Provides the interface to apply a linear operator defined by a
           Teuchos::SerialDenseMatrix.

    ---
*/


namespace ROL {

template <class Ordinal, class Real>
class TeuchosLinearOperator : public LinearOperator<Real> {
private:
  const Ptr<const Teuchos::SerialDenseMatrix<Ordinal,Real>> A_;
  LAPACK<Ordinal,Real> lapack_;
  Ptr<Teuchos::SerialDenseMatrix<Ordinal,Real>> Awork_;
  std::vector<Ordinal> ipiv_;

  Ptr<const Teuchos::SerialDenseVector<Ordinal,Real>> getVector( const Vector<Real>& x ) const {
    return dynamic_cast<const TeuchosVector<Ordinal,Real>&>(x).getVector();
  }
    
  Ptr<Teuchos::SerialDenseVector<Ordinal,Real>> getVector( Vector<Real>& x ) const {
    return dynamic_cast<TeuchosVector<Ordinal,Real>&>(x).getVector();
  }

  void inverse(Teuchos::SerialDenseVector<Ordinal,Real> &x, const char TRANS) const {
    Ordinal  N     = A_->numRows();
    Ordinal  NRHS  = 1;
    Ordinal  LDA   = N;
    Ordinal  LDB   = N;
    Ordinal  INFO  = 0;
    lapack_.GETRS(TRANS,N,NRHS,Awork_->values(),LDA,&ipiv_[0],x.values(),LDB,&INFO);
  }

public:

  TeuchosLinearOperator(const Ptr<const Teuchos::SerialDenseMatrix<Ordinal,Real>> &A)
    : A_(A) {
    Ordinal nrows = A_->numRows();
    Ordinal ncols = A_->numCols();
    Awork_ = makePtr<Teuchos::SerialDenseMatrix<Ordinal,Real>>(nrows,ncols);
    ipiv_.resize(std::min(nrows,ncols));
    for (Ordinal i = 0; i < nrows; ++i) {
      for (Ordinal j = 0; j < ncols; ++j) {
        (*Awork_)(i,j)  = (*A_)(i,j);
      }
    }
    Ordinal info(0);
    lapack_.GETRF(nrows,ncols,Awork_->values(),nrows,&ipiv_[0],&info);
  }

  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Ptr<Teuchos::SerialDenseVector<Ordinal,Real>>       Hv_data = getVector(Hv);
    Ptr<const Teuchos::SerialDenseVector<Ordinal,Real>>  v_data = getVector(v);
    Ordinal nrows = A_->numRows();
    Ordinal ncols = A_->numCols();
    for (Ordinal i = 0; i < nrows; ++i) {
      (*Hv_data)(i) = static_cast<Real>(0);
      for (Ordinal j = 0; j < ncols; ++j) {
        (*Hv_data)(i) += (*A_)(i,j)*(*v_data)(j);
      }
    }
  }

  virtual void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Ptr<Teuchos::SerialDenseVector<Ordinal,Real>>       Hv_data = getVector(Hv);
    Ptr<const Teuchos::SerialDenseVector<Ordinal,Real>>  v_data = getVector(v);
    Ordinal nrows = A_->numRows();
    for (Ordinal i = 0; i < nrows; ++i) {
      (*Hv_data)(i) = (*v_data)(i);
    }
    inverse(*Hv_data,'N');
  }

  void applyAdjoint( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Ptr<Teuchos::SerialDenseVector<Ordinal,Real>>       Hv_data = getVector(Hv);
    Ptr<const Teuchos::SerialDenseVector<Ordinal,Real>>  v_data = getVector(v);
    Ordinal nrows = A_->numRows();
    Ordinal ncols = A_->numCols();
    for (Ordinal j = 0; j < ncols; ++j) {
      (*Hv_data)(j) = static_cast<Real>(0);
      for (Ordinal i = 0; i < nrows; ++i) {
        (*Hv_data)(j) += (*A_)(i,j)*(*v_data)(i);
      }
    }
  }

  virtual void applyAdjointInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Ptr<Teuchos::SerialDenseVector<Ordinal,Real>>       Hv_data = getVector(Hv);
    Ptr<const Teuchos::SerialDenseVector<Ordinal,Real>>  v_data = getVector(v);
    Ordinal nrows = A_->numRows();
    for (Ordinal i = 0; i < nrows; ++i) {
      (*Hv_data)(i) = (*v_data)(i);
    }
    inverse(*Hv_data,'T');
  }

}; // class TeuchosLinearOperator

} // namespace ROL

#endif
