// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
