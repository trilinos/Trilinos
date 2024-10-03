// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_STD_COVARIANCE_OPERATOR_DEF_HPP
#define ROL_OED_STD_COVARIANCE_OPERATOR_DEF_HPP

namespace ROL {
namespace OED {

template<typename Real>
void StdMomentOperator<Real>::initialize(int nfactors) {
  startTimer("initialize");
  M_.shape(nfactors,nfactors);
  Minv_.shape(nfactors,nfactors);
  U_.shape(nfactors,nfactors);
  V_.shape(nfactors,nfactors);
  sval_.resize(nfactors);
  sval_.assign(nfactors,static_cast<Real>(0));
  stopTimer("initialize");
}

template<typename Real>
void StdMomentOperator<Real>::build(const Vector<Real> &p) {
  startTimer("build");
  if (!isBuilt_) {
    const Real zero(0), one(1);
    const int nsamples = Xdata_[0].numRows();
    const int nfactors = Xdata_[0].numCols();
    M_.putScalar(zero);
    LA::Matrix<Real> X(nsamples,nfactors);
    LA::Matrix<Real> M0(nfactors,nfactors);
    LA::Matrix<Real> Mtmp(nfactors,nfactors);
    Real wt(0);
    for (int k = 0; k < nobs_; ++k) {
      X = Xdata_[k];
      for (int i = 0; i < nsamples; ++i) {
        wt = std::sqrt(get(p,i)*getNoise(i));
        for (int j = 0; j < nfactors; ++j) X(i,j) *= wt;
      }
      M0.putScalar(zero);
      blas_->SYRK(Teuchos::UPPER_TRI,Teuchos::TRANS,nfactors,nsamples,
                  one,X.values(),nsamples,zero,M0.values(),nfactors);
      Mtmp += M0;
    }
    sumAll(p,Mtmp.values(),M_.values(),nfactors*nfactors);
    if (isPset_) blas_->AXPY(nfactors*nfactors,one,P_.values(),1,M_.values(),1);
    isBuilt_      = true;
    isFactorized_ = false;
    useSVD_       = false;
  }
  stopTimer("build");
}

template<typename Real>
void StdMomentOperator<Real>::factorize(const Vector<Real> &p) {
  startTimer("factorize");
  if (isBuilt_) {
    if (!isFactorized_) {
      const int nfactors = M_.numRows();
      Minv_.assign(M_);
      int info;
      // Factorize matrix
      lapack_->POTRF('U',nfactors,Minv_.values(),nfactors,&info);
      if (info != 0) {
        useSVD_ = true;
        //std::stringstream ss;
        //ss << ">>> OED::StdMomentOperator::factorize : "
        //   << "POTRF exited with INFO equal to " << info << "!" << std::endl
        //   << "  < 0: if INFO = -i, the i-th argument had an illegal value" << std::endl
        //   << "  > 0: if INFO = i, the leading minor of order i is not positive definite";
        //throw Exception::NotImplemented(ss.str());
      }
      if (useSVD_) {
        Minv_.assign(M_);
        int lwork = -1;
        std::vector<Real> work(1), rwork(5*nfactors);
        lapack_->GESVD('A','A',nfactors,nfactors,Minv_.values(),nfactors,
                       &sval_[0],U_.values(),nfactors,V_.values(),nfactors,
                       &work[0],lwork,&rwork[0],&info);
        if (info != 0) {
          std::stringstream ss;
          ss << ">>> OED::StdMomentOperator::factorize : "
             << "GESVD exited with INFO equal to " << info << "!" << std::endl
             << "  < 0: if INFO = -i, the i-th argument had an illegal value" << std::endl
             << "  > 0: if DBDSQR did not converge, INFO specifies how many" << std::endl
             << "       superdiagonals of an intermediate bidiagonal form B" << std::endl
             << "       did not converge to zero";
          throw Exception::NotImplemented(ss.str());
        }
        lwork = work[0];
        work.resize(lwork); work.assign(lwork,static_cast<Real>(0));
        lapack_->GESVD('A','A',nfactors,nfactors,Minv_.values(),nfactors,
                       &sval_[0],U_.values(),nfactors,V_.values(),nfactors,
                       &work[0],lwork,&rwork[0],&info);
        if (info != 0) {
          std::stringstream ss;
          ss << ">>> OED::StdMomentOperator::factorize : "
             << "GESVD exited with INFO equal to " << info << "!" << std::endl
             << "  < 0: if INFO = -i, the i-th argument had an illegal value" << std::endl
             << "  > 0: if DBDSQR did not converge, INFO specifies how many" << std::endl
             << "       superdiagonals of an intermediate bidiagonal form B" << std::endl
             << "       did not converge to zero";
          throw Exception::NotImplemented(ss.str());
        }
      }
      isFactorized_ = true;
    }
  }
  else {
    std::stringstream ss;
    ss << ">>> OED::StdMomentOperator::factorize : "
       << "Matrix has not been built!";
    throw Exception::NotImplemented(ss.str());
  }
  stopTimer("factorize");
}

/***************************************************************************/
/* Begin Accessor Functions                                                */
/***************************************************************************/
template<typename Real>
std::vector<Real>& StdMomentOperator<Real>::getData(Vector<Real> &x) const {
  return *dynamic_cast<StdVector<Real>&>(x).getVector();
}

template<typename Real>
const std::vector<Real>& StdMomentOperator<Real>::getConstData(const Vector<Real> &x) const {
  return *dynamic_cast<const StdVector<Real>&>(x).getVector();
}
/***************************************************************************/
/* End Accessor Functions                                                  */
/***************************************************************************/

template<typename Real>
StdMomentOperator<Real>::StdMomentOperator(RegressionType regType,
                  bool homNoise,
                  const Ptr<Noise<Real>> &noise)
  : MomentOperator<Real>(regType,homNoise,noise),
    lapack_(makePtr<LAPACK<int,Real>>()),
    blas_(makePtr<Teuchos::BLAS<int,Real>>()),
    isBuilt_(false), isFactorized_(false), isSet_(false),
    isFullSet_(false), useSVD_(false), isPset_(false) {
  ProfiledClass<Real,std::string>::rename("OED::StdMomentOperator");
}

template<typename Real>
Ptr<MomentOperator<Real>> StdMomentOperator<Real>::clone() const {
  RegressionType type;
  bool hom;
  Ptr<Noise<Real>> noise;
  MomentOperator<Real>::getRegressionInfo(type,hom,noise);
  return makePtr<StdMomentOperator<Real>>(type,hom,noise);
}

template<typename Real>
void StdMomentOperator<Real>::update(const Vector<Real> &p, UpdateType type, int iter) {
  if ( type == UpdateType::Initial
    || type == UpdateType::Trial
    || type == UpdateType::Temp ) {
    isBuilt_ = false;
  }
}

template<typename Real>
void StdMomentOperator<Real>::applyInverse(Vector<Real> &Mx,
            const Vector<Real> &x,
            const Vector<Real> &p) {
  startTimer("applyInverse");
  if (isSet_) {
    build(p);
    factorize(p);

    Mx.set(x);
    std::vector<Real>       &Mxdata = getData(Mx);
    const std::vector<Real> & xdata = getConstData(x);

    const int nfactors = M_.numRows();
    int info;

    if (!useSVD_) {
      // Solve triangular systems
      lapack_->POTRS('U',nfactors,1,Minv_.values(),nfactors,
                     &Mxdata[0],nfactors,&info);
      if (info != 0) {
        std::stringstream ss;
        ss << ">>> OED::StdMomentOperator::applyInverse : "
           << "POTRS exited with INFO equal to " << info << "!" << std::endl
           << "  < 0: if INFO = -i, the i-th argument had an illegal value";
        throw Exception::NotImplemented(ss.str());
      }
      // Refine solution
      std::vector<Real> ferr(1), berr(1), work(3*nfactors);
      std::vector<int> lwork(nfactors);
      lapack_->PORFS('U',nfactors,1,M_.values(),nfactors,Minv_.values(),
                     nfactors,&xdata[0],nfactors,&Mxdata[0],nfactors,
                     &ferr[0],&berr[0],&work[0],&lwork[0],&info);
      if (info != 0) {
        std::stringstream ss;
        ss << ">>> OED::StdMomentOperator::applyInverse : "
           << "PORFS exited with INFO equal to " << info << "!" << std::endl
           << "  < 0: if INFO = -i, the i-th argument had an illegal value";
        throw Exception::NotImplemented(ss.str());
      }
    }
    else {
      const Real zero(0), one(1);
      std::vector<Real> Ux(nfactors);
      blas_->GEMV(Teuchos::TRANS,nfactors,nfactors,one,U_.values(),nfactors,
                  &xdata[0],1,zero,&Ux[0],1);
      Real maxSV(-1);
      for (int i = 0; i < nfactors; ++i) {
        maxSV = (maxSV < sval_[i] ? sval_[i] : maxSV);
      }
      const Real tol = ROL_EPSILON<Real>()*static_cast<Real>(nfactors) * maxSV;
      for (int i = 0; i < nfactors; ++i) {
        if (sval_[i] > tol) {
          Ux[i] /= sval_[i];
        }
        else {
          Ux[i] = zero;
        }
      }
      blas_->GEMV(Teuchos::TRANS,nfactors,nfactors,one,V_.values(),nfactors,
                  &Ux[0],1,zero,&Mxdata[0],1);
    }
  }
  else {
    std::stringstream ss;
    ss << ">>> OED::StdMomentOperator::applyInverse : "
       << "Factors have not been set!";
    throw Exception::NotImplemented(ss.str());
  }
  stopTimer("applyInverse");
}

template<typename Real>
void StdMomentOperator<Real>::apply(Vector<Real> &Mc,
           const Vector<Real> &c,
           const Vector<Real> &p) {
  startTimer("apply");
  if (isSet_) {
    build(p);
    const Real zero(0), one(1);
    const int M = M_.numRows();
    std::vector<Real>       &Mcdata = getData(Mc);
    const std::vector<Real> & cdata = getConstData(c);
    blas_->SYMM(Teuchos::LEFT_SIDE,Teuchos::UPPER_TRI,M,1,
                one,M_.values(),M,&cdata[0],M,zero,
                &Mcdata[0],M);
  }
  else {
    std::stringstream ss;
    ss << ">>> OED::StdMomentOperator::apply : "
       << "Factors have not been set!";
    throw Exception::NotImplemented(ss.str());
  }
  stopTimer("apply");
}

template<typename Real>
void StdMomentOperator<Real>::applyDeriv(Vector<Real> &Mc,
                const Vector<Real> &c,
                const Vector<Real> &p) {
  startTimer("applyDeriv");
  if (isSet_) {
    Mc.zero();
    const Real zero(0), one(1);
    const int M = Xdata_[0].numRows();
    const int N = Xdata_[0].numCols();
    std::vector<Real> Xctmp(M), mMctmp(N), mMcdata(N);
    std::vector<Real>       &Mcdata = getData(Mc);
    const std::vector<Real> & cdata = getConstData(c);
    for (int k = 0; k < nobs_; ++k) {
      blas_->GEMV(Teuchos::NO_TRANS,M,N,one,Xdata_[k].values(),M,
                  &cdata[0],1,zero,&Xctmp[0],1);
      for (int i = 0; i < M; ++i) Xctmp[i] *= get(p,i)*getNoise(i);
      blas_->GEMV(Teuchos::TRANS,M,N,one,Xdata_[k].values(),M,
                  &Xctmp[0],1,zero,&mMctmp[0],1);
      for (int i = 0; i < N; ++i) mMcdata[i] += mMctmp[i];
    }
    sumAll(p,&mMcdata[0],&Mcdata[0],N);
  }
  else {
    std::stringstream ss;
    ss << ">>> OED::StdMomentOperator::applyDeriv : "
       << "Factors have not been set!";
    throw Exception::NotImplemented(ss.str());
  }
  stopTimer("applyDeriv");
}

template<typename Real>
void StdMomentOperator<Real>::applySampleMatrices(Vector<Real> &uXv, const Vector<Real> &u, const Vector<Real> &v) {
  startTimer("applySampleMatrices");
  if (isSet_) {
    const Real zero(0), one(1);
    const int M = Xdata_[0].numRows();
    const int N = Xdata_[0].numCols();
    uXv.zero();
    std::vector<Real>        Xutmp(M,0), Xvtmp(M,0);
    std::vector<Real>       &Xudata = getLocalDesign(uXv);
    const std::vector<Real> &udata  = getConstData(u);
    const std::vector<Real> &vdata  = getConstData(v);
    for (int k = 0; k < nobs_; ++k) {
      blas_->GEMV(Teuchos::NO_TRANS,M,N,one,Xdata_[k].values(),M,
                  &udata[0],1,zero,&Xutmp[0],1);
      blas_->GEMV(Teuchos::NO_TRANS,M,N,one,Xdata_[k].values(),M,
                  &vdata[0],1,zero,&Xvtmp[0],1);
      for (int i = 0; i < M; ++i) Xudata[i] += Xutmp[i]*Xvtmp[i];
    }
    for (int i = 0; i < M; ++i) Xudata[i] *= getNoise(i);
  }
  else {
    std::stringstream ss;
    ss << ">>> OED::StdMomentOperator::applySampleMatrices : "
       << "Factors have not been set!";
    throw Exception::NotImplemented(ss.str());
  }
  stopTimer("applySampleMatrices");
}

template<typename Real>
void StdMomentOperator<Real>::setFactors(const Ptr<Factors<Real>> &factors) {
  startTimer("setFactors");
  MomentOperator<Real>::setFactors(factors);
  int nrows = factors->numMySamples(), ncols = factors->numFactors();
  initialize(ncols);
  Ptr<Vector<Real>> obs = factors->createObservationVector(true);
  nobs_ = obs->dimension();
  Xdata_.resize(nobs_);
  Xpred_.shape(nrows,ncols);
  if (nobs_ == 1) {
    Xdata_[0].shape(nrows,ncols);
    for (int i = 0; i < nrows; ++i) {
      const std::vector<Real> &xdata = getConstData(*factors->get(i));
      for (int j = 0; j < ncols; ++j) Xdata_[0](i,j) = xdata[j];
    }
    Xpred_ = Xdata_[0];
  }
  else {
    Ptr<Vector<Real>> c = obs->clone();
    factors->getPredictionVector(*c);
    const std::vector<Real> &cdata = getConstData(*c);
    for (int k = 0; k < nobs_; ++k) {
      Xdata_[k].shape(nrows,ncols);
      factors->setPredictionVector(*obs->basis(k));
      for (int i = 0; i < nrows; ++i) {
        const std::vector<Real> &xdata = getConstData(*factors->get(i));
        for (int j = 0; j < ncols; ++j) Xdata_[k](i,j) = xdata[j];
      }
      blas_->AXPY(ncols*nrows,cdata[k],Xdata_[k].values(),1,Xpred_.values(),1);
      //Xpred_ += cdata[k]*Xdata_[k];
    }
    factors->setPredictionVector(*c);
  }
  if (isPset_) {
    P_.shape(ncols,ncols);
    Ptr<Vector<Real>> ei  = factors->createParameterVector();
    Ptr<Vector<Real>> Pei = factors->createParameterVector();
    for (int i = 0; i < ncols; ++i) {
      Pei->zero();
      ei->set(*Pei->basis(i));
      applyPerturbation(*Pei,*ei);
      const std::vector<Real> &Pei_data = getConstData(*Pei);
      // Build upper triangle of perturbation
      for (int j = 0; j <= i; ++j) P_(j,i) = Pei_data[j];
    }
  }
  isSet_ = true;
  stopTimer("setFactors");
}

template<typename Real>
void StdMomentOperator<Real>::setPerturbation(const Ptr<LinearOperator<Real>> &pOp) {
  MomentOperator<Real>::setPerturbation(pOp);
  isPset_ = true;
}

template<typename Real>
Real StdMomentOperator<Real>::logDeterminant(const Vector<Real> &z) {
  startTimer("logDeterminant");
  Real val(0);
  if (isSet_) {
    build(z);
    factorize(z);
    const int nfactors = M_.numRows();
    for (int j = 0; j < nfactors; ++j) {
      val += std::log(!useSVD_ ? Minv_(j,j) : sval_[j]);
    }
  }
  else {
    std::stringstream ss;
    ss << ">>> OED::StdMomentOperator::logDeterminant : "
       << "Factors have not been set!";
    throw Exception::NotImplemented(ss.str());
  }
  if (!useSVD_) {
    val *= static_cast<Real>(2);
  }
  stopTimer("logDeterminant");
  return val;
}

} // End OED Namespace
} // End ROL Namespace

#endif
