// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_COVARIANCE_OPERATOR_HPP
#define ROL_OED_COVARIANCE_OPERATOR_HPP

#include "ROL_Ptr.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_ProbabilityVector.hpp"
#include "ROL_ConjugateResiduals.hpp"
#include "ROL_UpdateType.hpp"
#include "ROL_Types.hpp"
#include "ROL_OED_Factors.hpp"
#include "ROL_OED_Noise.hpp"
#include "ROL_OED_ProfiledClass.hpp"
#include <vector>

namespace ROL::OED {

enum RegressionType {
  LEASTSQUARES = 0,
  WEIGHTEDLEASTSQUARES,
  QUANTILE,
  SUPERQUANTILE,
  EXPONENTIAL,
  REGTYPE_LAST
};

inline std::string RegressionTypeToString(RegressionType tr) {
  std::string retString;
  switch(tr) {
    case LEASTSQUARES:         retString = "Least Squares";          break;
    case WEIGHTEDLEASTSQUARES: retString = "Weighted Least Squares"; break;
    case QUANTILE:             retString = "Quantile";               break;
    case SUPERQUANTILE:        retString = "Superquantile";          break;
    case EXPONENTIAL:          retString = "Exponential";            break;
    case REGTYPE_LAST:         retString = "Last Type (Dummy)";      break;
    default:                   retString = "INVALID RegressionType";
  }
  return retString;
}

inline int isValidRegressionType(RegressionType ls) {
  return( (ls == LEASTSQUARES)         ||
          (ls == WEIGHTEDLEASTSQUARES) ||
          (ls == QUANTILE)             ||
          (ls == SUPERQUANTILE)        ||
          (ls == EXPONENTIAL) );
}

inline RegressionType & operator++(RegressionType &type) {
  return type = static_cast<RegressionType>(type+1);
}

inline RegressionType operator++(RegressionType &type, int) {
  RegressionType oldval = type;
  ++type;
  return oldval;
}

inline RegressionType & operator--(RegressionType &type) {
  return type = static_cast<RegressionType>(type-1);
}

inline RegressionType operator--(RegressionType &type, int) {
  RegressionType oldval = type;
  --type;
  return oldval;
}

inline RegressionType StringToRegressionType(std::string s) {
  s = removeStringFormat(s);
  for ( RegressionType st = LEASTSQUARES; st < REGTYPE_LAST; ++st ) {
    if ( !s.compare(removeStringFormat(RegressionTypeToString(st))) ) {
      return st;
    }
  }
  return REGTYPE_LAST;
}

template<typename Real>
class MomentOperator : public ProfiledClass<Real,std::string> {
private:
  RegressionType regType_;
  bool homNoise_;
  Ptr<Noise<Real>> noise_;
  int matNum_;
  Ptr<Vector<Real>> sum_, Jx_, SJx_;
  Ptr<Krylov<Real>> krylov_;
  bool setSolver_;
  Ptr<LinearOperator<Real>> pOp_;

  class moment : public LinearOperator<Real> {
  private:
    const Ptr<Factors<Real>> factors_;
    const RegressionType regType_;
    const bool homNoise_;
    const Ptr<Noise<Real>> noise_;
    const int matNum_;
    const Ptr<LinearOperator<Real>> pOp_;

    Ptr<Vector<Real>> sum_, Jx_, SJx_;
    Ptr<const ProbabilityVector<Real>> p_;

    void applyNoise(Vector<Real>& Nx, const Vector<Real>& x, int i) const {
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

  public:
    moment(const Ptr<Factors<Real>>& factors,
           RegressionType regType, int matNum,
           bool homNoise,
           const Ptr<Noise<Real>>& noise,
           const Ptr<LinearOperator<Real>>& pOp)
      : factors_(factors), regType_(regType), homNoise_(homNoise),
        noise_(noise), matNum_(matNum), pOp_(pOp) {
      Jx_ = factors_->createObservationVector(false);
      SJx_ = factors_->createObservationVector(true);
      sum_ = factors_->createParameterVector(true);
    }
    void setProbabilityVector(const ProbabilityVector<Real> &p) {
      p_ = makePtrFromRef(p);
    }
    void apply(Vector<Real> &Mx, const Vector<Real> &x, Real &tol) const {
      int nsamples = factors_->numMySamples();
      sum_->zero();
      for (int i = 0; i < nsamples; ++i) {
        factors_->apply(*Jx_, x, i);          // Compute J_i*x
        applyNoise(*SJx_, *Jx_, i);           // Compute S_i*J_i*x
        factors_->applyAdjoint(Mx, *SJx_, i); // Compute J_i'*S_i*J_i*x
        sum_->axpy(p_->getProbability(i),Mx);
      }
      Mx.zero();
      p_->getBatchManager()->sumAll(*sum_,Mx);
      if (pOp_ != nullPtr) {
        sum_->zero();
        pOp_->apply(*sum_,x,tol);
        Mx.plus(*sum_);
      }
    }
  };

  class precond : public LinearOperator<Real> {
  public:
    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      Hv.set(v.dual());
    }
    void applyInverse(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      Hv.set(v.dual());
    }
  };

  Ptr<moment>  M_;
  Ptr<precond> P_;

protected:
  Ptr<Factors<Real>> factors_;

  /***************************************************************************/
  /* Begin Accessor Functions                                                */
  /***************************************************************************/
  Real get(const Vector<Real> &x, int i) const;
  void set(Vector<Real> &x, int i, Real val) const;
  std::vector<Real>& getLocalDesign(Vector<Real> &p) const;
  const std::vector<Real>& getConstLocalDesign(const Vector<Real> &p) const;
  void sumAll(const Vector<Real> &p, Real* input, Real* output, int size) const;
  void applyPerturbation(Vector<Real> &Px, const Vector<Real> &x) const;
  /***************************************************************************/
  /* End Accessor Functions                                                  */
  /***************************************************************************/

  using ProfiledClass<Real,std::string>::startTimer;
  using ProfiledClass<Real,std::string>::stopTimer;

public:
  virtual ~MomentOperator() {}
  MomentOperator(RegressionType regType=LEASTSQUARES,
                 const Ptr<Noise<Real>>& noise=nullPtr);
  MomentOperator(std::string regType, const Ptr<Noise<Real>>& noise = nullPtr);

  virtual Ptr<MomentOperator<Real>> clone() const;
  void setMatrixNumber(int matNum);
  virtual void update(const Vector<Real> &p, UpdateType type, int iter = -1) {}

  virtual void generateFactors(const Ptr<Constraint<Real>>      &model,
                               const Ptr<Vector<Real>>          &theta,
                               const Ptr<Vector<Real>>          &obs,
                               const Ptr<SampleGenerator<Real>> &sampler);
  virtual void generateFactors(const Ptr<Objective<Real>>       &model,
                               const Ptr<Vector<Real>>          &theta,
                               const Ptr<SampleGenerator<Real>> &sampler);

  virtual void setFactors(const Ptr<Factors<Real>> &factors);
  virtual void setPerturbation(const Ptr<LinearOperator<Real>> &pOp);
  void setNoise(const Ptr<Noise<Real>> &noise, bool isHom = false);

  void getRegressionInfo(RegressionType &regType, bool &homNoise,
                         Ptr<Noise<Real>> &noise) const;
  virtual const Ptr<Factors<Real>> getFactors() const;
  virtual const Ptr<LinearOperator<Real>> getPerturbation() const;

  // Compute M(p)x where M(p) = p_1 X_1 S_1 X_1' + ... + p_N X_N S_N X_N'
  virtual void apply(Vector<Real> &Mx,
                     const Vector<Real> &x,
                     const Vector<Real> &p);

  // Compute M(q)x where M(q) = q_1 X_1 S_1 X_1' + ... + q_N X_N S_N X_N'
  // This function is distinct from apply to allow the user to store M(p)
  virtual void applyDeriv(Vector<Real> &Mx,
                          const Vector<Real> &x,
                          const Vector<Real> &q);

  // Compute inv(M(p))x where M(p) = p_1 X_1 S_1 X_1' + ... + p_N X_N S_N X_N'
  virtual void applyInverse(Vector<Real> &Mx,
                            const Vector<Real> &x,
                            const Vector<Real> &p);

  virtual void applySampleMatrices(Vector<Real> &uXv, const Vector<Real> &u, const Vector<Real> &v);

  void applyNoise(Vector<Real>& Nx, const Vector<Real>& x, int i) const;

  virtual Real logDeterminant(const Vector<Real> &z);



}; // class MomentOperator

} // End ROL::OED Namespace

#include "ROL_OED_MomentOperator_Def.hpp"

#endif
