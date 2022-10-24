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

#ifndef ROL_OED_COVARIANCE_OPERATOR_DEF_HPP
#define ROL_OED_COVARIANCE_OPERATOR_DEF_HPP

namespace ROL {
namespace OED {

/***************************************************************************/
/* Begin Accessor Functions                                                */
/***************************************************************************/
template<typename Real>
Real MomentOperator<Real>::get(const Vector<Real> &x, int i) const {
  return dynamic_cast<const ProbabilityVector<Real>&>(x).getProbability(i);
}

template<typename Real>
void MomentOperator<Real>::set(Vector<Real> &x, int i, Real val) const {
  dynamic_cast<ProbabilityVector<Real>&>(x).setProbability(i,val);
}

template<typename Real>
std::vector<Real>& MomentOperator<Real>::getLocalDesign(Vector<Real> &p) const {
  return *dynamic_cast<StdVector<Real>&>(p).getVector();
}

template<typename Real>
const std::vector<Real>& MomentOperator<Real>::getConstLocalDesign(const Vector<Real> &p) const {
  return *dynamic_cast<const StdVector<Real>&>(p).getVector();
}

template<typename Real>
void MomentOperator<Real>::sumAll(const Vector<Real> &p, Real* input, Real* output, int size) const {
  dynamic_cast<const ProbabilityVector<Real>&>(p).getBatchManager()->sumAll(input,output,size);
}

template<typename Real>
void MomentOperator<Real>::applyPerturbation(Vector<Real> &Px, const Vector<Real> &x) const {
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  pOp_->apply(Px,x,tol);
}
/***************************************************************************/
/* End Accessor Functions                                                  */
/***************************************************************************/

template<typename Real>
MomentOperator<Real>::MomentOperator(RegressionType regType,
               bool homNoise,
               const Ptr<Noise<Real>> &noise)
  : ProfiledClass<Real,std::string>("OED::MomentOperator"),
    regType_(regType), homNoise_(homNoise), noise_(noise), matNum_(0),
    sum_(nullPtr), setSolver_(true) {}

template<typename Real>
Ptr<MomentOperator<Real>> MomentOperator<Real>::clone() const {
  RegressionType type;
  bool hom;
  Ptr<Noise<Real>> noise;
  MomentOperator<Real>::getRegressionInfo(type,hom,noise);
  return makePtr<MomentOperator<Real>>(type,hom,noise);
}

template<typename Real>
void MomentOperator<Real>::setMatrixNumber(int matNum) {
  matNum_ = matNum;
}

template<typename Real>
void MomentOperator<Real>::setFactors(const Ptr<Factors<Real>> &factors) {
  factors_ = factors;
}

template<typename Real>
void MomentOperator<Real>::setPerturbation(const Ptr<LinearOperator<Real>> &pOp) {
  pOp_ = pOp;
}

// Compute M(p)x where M(p) = p_1 X_1 S_1 X_1' + ... + p_N X_N S_N X_N'
template<typename Real>
void MomentOperator<Real>::apply(Vector<Real> &Mx,
                   const Vector<Real> &x,
                   const Vector<Real> &p) {
  startTimer("apply");
  if (sum_ == nullPtr) sum_ = Mx.clone();
  int nsamples = factors_->numMySamples();
  if (nsamples > 0) {
    sum_->zero();
    for (int i = 0; i < nsamples; ++i) {
      factors_->applyProduct(Mx,x,i);
      sum_->axpy(get(p,i)*getNoise(i),Mx);
    }
    Mx.zero();
    dynamic_cast<const ProbabilityVector<Real>&>(p).getBatchManager()->sumAll(*sum_,Mx);
    if (pOp_ != nullPtr) {
      sum_->zero();
      applyPerturbation(*sum_,x);
      Mx.plus(*sum_);
    }
  }
  else {
    throw Exception::NotImplemented(">>> OED::MomentOperator::apply : Factors have not been set!");
  }
  stopTimer("apply");
}

// Compute M(q)x where M(q) = q_1 X_1 S_1 X_1' + ... + q_N X_N S_N X_N'
// This function is distinct from apply to allow the user to store M(p)
template<typename Real>
void MomentOperator<Real>::applyDeriv(Vector<Real> &Mx,
                        const Vector<Real> &x,
                        const Vector<Real> &q) {
  startTimer("applyDeriv");
  if (sum_ == nullPtr) sum_ = Mx.clone();
  int nsamples = factors_->numMySamples();
  if (nsamples > 0) {
    sum_->zero();
    for (int i = 0; i < nsamples; ++i) {
      factors_->applyProduct(Mx,x,i);
      sum_->axpy(get(q,i)*getNoise(i),Mx);
    }
    Mx.zero(); 
    dynamic_cast<const ProbabilityVector<Real>&>(q).getBatchManager()->sumAll(*sum_,Mx);
  }
  else {
    throw Exception::NotImplemented(">>> OED::MomentOperator::applyDeriv : Factors have not been set!");
  }
  stopTimer("applyDeriv");
}

// Compute inv(M(p))x where M(p) = p_1 X_1 S_1 X_1' + ... + p_N X_N S_N X_N'
template<typename Real>
void MomentOperator<Real>::applyInverse(Vector<Real> &Mx,
                          const Vector<Real> &x,
                          const Vector<Real> &p) {
  startTimer("applyInverse");
  // Set up Krylov solver
  if (setSolver_) {
    krylov_ = makePtr<ConjugateResiduals<Real>>(1e-10,1e0,1000);
    M_ = makePtr<moment>(factors_,pOp_);
    P_ = makePtr<precond>();
    setSolver_ = false;
  }
  M_->setProbabilityVector(static_cast<const ProbabilityVector<Real>&>(p));
  int iter(0), flag(0);
  krylov_->run(Mx,*M_,x,*P_,iter,flag);
  stopTimer("applyInverse");
}

template<typename Real>
void MomentOperator<Real>::applySampleMatrices(Vector<Real> &uXv, const Vector<Real> &u, const Vector<Real> &v) {
  startTimer("applySampleMatrices");
  int nsamples = factors_->numMySamples();
  if (nsamples > 0) {
    Real val(0);
    for (int i = 0; i < nsamples; ++i) {
      val = getNoise(i)*factors_->applyProduct2(u,v,i);
      set(uXv,i,val);
    }
  }
  else {
    throw Exception::NotImplemented(">>> OED::MomentOperator::applySampleMatrices : Factors have not been set!");
  }
  stopTimer("applySampleMatrices");
}

template<typename Real>
void MomentOperator<Real>::setNoise(const Ptr<Noise<Real>> &noise, bool isHom) {
  noise_ = noise;
  homNoise_ = isHom;
}

template<typename Real>
Real MomentOperator<Real>::getNoise(int k) const {
  Real val(1);
  if (noise_ != nullPtr && !homNoise_) {
    const Real one(1);
    Real noise = noise_->evaluate(factors_->getSample(k));
    switch(regType_) {
      case(LEASTSQUARES)  : val = (matNum_==1 ? one       : std::pow(noise,2));                        break;
      case(QUANTILE)      : val = (matNum_==1 ? one/noise : one);                                      break;
      case(SUPERQUANTILE) : val = (matNum_==1 ? one/noise : one);                                      break;
      case(EXPONENTIAL)   : val = (matNum_==1 ? one       : noise_->mgf2(factors_->getSample(k))-one); break;
      case(REGTYPE_LAST)  : val = one;
    }
  }
  return val;
}

template<typename Real>
void MomentOperator<Real>::getRegressionInfo(RegressionType &regType, bool &homNoise,
                       Ptr<Noise<Real>> &noise) const {
  regType  = regType_;
  homNoise = homNoise_;
  noise    = noise_;
}

template<typename Real>
Real MomentOperator<Real>::logDeterminant(const Vector<Real> &z) {
  throw Exception::NotImplemented(">>> OED::MomentOperator::logDeterminant is not implemented!");
}

} // End OED Namespace
} // End ROL Namespace

#endif
