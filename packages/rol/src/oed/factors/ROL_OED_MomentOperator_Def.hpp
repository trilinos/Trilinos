// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
                                     const Ptr<Noise<Real>>& noise)
  : ProfiledClass<Real,std::string>("OED::MomentOperator"),
    regType_(regType),
    homNoise_(noise!=nullPtr ? noise->isHomoscedastic()&&!noise->isCorrelated() : true),
    noise_(noise), matNum_(0), sum_(nullPtr), setSolver_(true) {
  if (!isValidRegressionType(regType)) {
    std::stringstream ss;
    ss << ">>> OED::MomentOperator::MomentOperator : "
       << regType << " is not a valid regression type!" << std::endl;
    throw Exception::NotImplemented(ss.str());
  }
}
    
template<typename Real>
MomentOperator<Real>::MomentOperator(std::string regType,
                                     const Ptr<Noise<Real>>& noise)
  : MomentOperator<Real>(StringToRegressionType(regType),noise) {}

template<typename Real>
Ptr<MomentOperator<Real>> MomentOperator<Real>::clone() const {
  RegressionType type;
  bool hom;
  Ptr<Noise<Real>> noise;
  MomentOperator<Real>::getRegressionInfo(type,hom,noise);
  return makePtr<MomentOperator<Real>>(type,noise);
}

template<typename Real>
void MomentOperator<Real>::setMatrixNumber(int matNum) {
  matNum_ = matNum;
}

template<typename Real>
void MomentOperator<Real>::setFactors(const Ptr<Factors<Real>>& factors) {
  factors_ = factors;
  Jx_  = factors_->createObservationVector(false);
  SJx_ = factors_->createObservationVector(true);
  sum_ = factors_->createParameterVector(true);
  const int odim = Jx_->dimension();
  if (odim>1 && regType_!=LEASTSQUARES && regType_!=WEIGHTEDLEASTSQUARES) {
    std::stringstream ss;
    ss << ">>> OED::MomentOperator::setFactors : "
       << RegressionTypeToString(regType_)
       << "  moment operator is not defined for factor range "
       << "dimension greater than one!" << std::endl
       << "  Input factor range dimension is " << odim << "!" << std::endl;
    throw Exception::NotImplemented(ss.str());
  }
}

template<typename Real>
void MomentOperator<Real>::generateFactors(const Ptr<Constraint<Real>>      &model,
                                           const Ptr<Vector<Real>>          &theta,
                                           const Ptr<Vector<Real>>          &obs,
                                           const Ptr<SampleGenerator<Real>> &sampler) {
  auto factors = makePtr<Factors<Real>>(model,theta,obs,sampler);
  setFactors(factors);
}

template<typename Real>
void MomentOperator<Real>::generateFactors(const Ptr<Objective<Real>>       &model,
                                           const Ptr<Vector<Real>>          &theta,
                                           const Ptr<SampleGenerator<Real>> &sampler) {
  auto factors = makePtr<Factors<Real>>(model,theta,sampler);
  setFactors(factors);
}

template<typename Real>
void MomentOperator<Real>::setPerturbation(const Ptr<LinearOperator<Real>>& pOp) {
  pOp_ = pOp;
}

template<typename Real>
const Ptr<Factors<Real>> MomentOperator<Real>::getFactors() const {
  return factors_;
}

template<typename Real>
const Ptr<LinearOperator<Real>> MomentOperator<Real>::getPerturbation() const {
  return pOp_;
}

// Compute M(p)x where M(p) = p_1 X_1 S_1 X_1' + ... + p_N X_N S_N X_N'
template<typename Real>
void MomentOperator<Real>::apply(Vector<Real>& Mx,
                           const Vector<Real>& x,
                           const Vector<Real>& p) {
  startTimer("apply");
  const int nsamples = (factors_==nullPtr ? 0 : factors_->numMySamples());
  if (nsamples > 0) {
    sum_->zero();
    for (int i = 0; i < nsamples; ++i) {
      factors_->apply(*Jx_, x, i);          // Compute J_i*x
      applyNoise(*SJx_, *Jx_, i);           // Compute S_i*J_i*x
      factors_->applyAdjoint(Mx, *SJx_, i); // Compute J_i'*S_i*J_i*x
      sum_->axpy(get(p,i), Mx);
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
void MomentOperator<Real>::applyDeriv(Vector<Real>& Mx,
                                const Vector<Real>& x,
                                const Vector<Real>& q) {
  startTimer("applyDeriv");
  int nsamples = factors_->numMySamples();
  if (nsamples > 0) {
    sum_->zero();
    for (int i = 0; i < nsamples; ++i) {
      factors_->apply(*Jx_, x, i);          // Compute J_i*x
      applyNoise(*SJx_, *Jx_, i);           // Compute S_i*J_i*x
      factors_->applyAdjoint(Mx, *SJx_, i); // Compute J_i'*S_i*J_i*x
      sum_->axpy(get(q,i), Mx);
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
void MomentOperator<Real>::applyInverse(Vector<Real>& Mx,
                                  const Vector<Real>& x,
                                  const Vector<Real>& p) {
  startTimer("applyInverse");
  // Set up Krylov solver
  if (setSolver_) {
    krylov_ = makePtr<ConjugateResiduals<Real>>(1e-10,1e0,1000);
    M_ = makePtr<moment>(factors_,regType_,matNum_,homNoise_,noise_,pOp_);
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
      factors_->apply(*Jx_, v, i); // Compute J_i*v
      applyNoise(*SJx_, *Jx_, i);  // Compute S_i*J_i*v
      factors_->apply(*Jx_, u, i); // Compute J_i*u
      val = SJx_->apply(*Jx_);
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
void MomentOperator<Real>::applyNoise(Vector<Real>& Nx, const Vector<Real>& x, int i) const {
  if (noise_ != nullPtr && !homNoise_) {
    auto pt = factors_->getSample(i);
    const Real one(1);
    if (regType_==LEASTSQUARES) {
      if (matNum_==1) Nx.set(x.dual());
      else noise_->apply(Nx,x,pt);
    }
    else if (regType_==WEIGHTEDLEASTSQUARES) {
      noise_->applyInverse(Nx,x,pt);
    }
    else {
      Real val(1);
      if (regType_==QUANTILE || regType_==SUPERQUANTILE) {
        Real noise = noise_->evaluate(pt);
        val = (matNum_==1 ? one/noise : one);
      }
      else if(regType_==EXPONENTIAL) {
        val = (matNum_==1 ? one : noise_->mgf2(pt)-one);
      }
      Nx.set(x.dual());
      Nx.scale(val);
    }
  }
  else {
    Nx.set(x);
  }
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
