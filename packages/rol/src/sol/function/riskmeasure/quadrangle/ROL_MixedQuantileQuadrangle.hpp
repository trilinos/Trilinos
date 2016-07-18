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

#ifndef ROL_MIXEDQUANTILEQUADRANGLE_HPP
#define ROL_MIXEDQUANTILEQUADRANGLE_HPP

#include "ROL_RiskMeasure.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_RiskVector.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"

/** @ingroup risk_group
    \class ROL::MixedQuantileQuadrangle
    \brief Provides an interface for a convex combination of
           conditional value-at-risks.

    The risk measure associated with the mixed-quantile quadrangle is defined
    as
    \f[
       \mathcal{R}(X) = \lambda_1 \mathrm{CVaR}_{\beta_1}(X)
         + \ldots + \lambda_n \mathrm{CVaR}_{\beta_n}(X)
    \f]
    where \f$0 \le \beta_1 \le \cdots \le \beta_n < 1\f$ and
    \f$0 \le \lambda_i\f$, \f$i=1,\ldots,n\f$, satisfies
    \f[
       \lambda_1 + \ldots + \lambda_n = 1.
    \f]
    Here, the conditional value-at-risk (CVaR) with confidence level
    \f$0\le \beta < 1\f$ is
    \f[
       \mathrm{CVaR}_\beta(X) = \inf_{t\in\mathbb{R}} \left\{
         t + \frac{1}{1-\beta} \mathbb{E}\left[(X-t)_+\right]
         \right\}
    \f]
    where \f$(x)_+ = \max\{0,x\}\f$.  If the distribution of \f$X\f$ is
    continuous, then \f$\mathrm{CVaR}_{\beta}(X)\f$ is the conditional
    expectation of \f$X\f$ exceeding the \f$\beta\f$-quantile of \f$X\f$ and
    the optimal \f$t\f$ is the \f$\beta\f$-quantile.
    Additionally, \f$\mathcal{R}\f$ is a law-invariant coherent risk measure.

    When using derivative-based optimization, the user can provide a smooth
    approximation of \f$(\cdot)_+\f$ using the ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class MixedQuantileQuadrangle : public RiskMeasure<Real> {
private:
  Teuchos::RCP<PlusFunction<Real> > plusFunction_;

  Teuchos::Array<Real> prob_;
  Teuchos::Array<Real> coeff_;

  Teuchos::RCP<Vector<Real> > dualVector_;
  std::vector<Real> xvar_;
  std::vector<Real> vvar_;

  std::vector<Real> vec_;

  int size_;

  bool firstReset_;

  void checkInputs(void) const {
    int pSize = prob_.size(), cSize = coeff_.size();
    TEUCHOS_TEST_FOR_EXCEPTION((pSize!=cSize),std::invalid_argument,
      ">>> ERROR (ROL::MixedQuantileQuadrangle): Probability and coefficient arrays have different sizes!");
    Real sum(0), zero(0), one(1);
    for (int i = 0; i < pSize; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION((prob_[i]>one || prob_[i]<zero), std::invalid_argument,
        ">>> ERROR (ROL::MixedQuantileQuadrangle): Element of probability array out of range!");
      TEUCHOS_TEST_FOR_EXCEPTION((coeff_[i]>one || coeff_[i]<zero), std::invalid_argument,
        ">>> ERROR (ROL::MixedQuantileQuadrangle): Element of coefficient array out of range!");
      sum += coeff_[i];
    }
    TEUCHOS_TEST_FOR_EXCEPTION((std::abs(sum-one) > std::sqrt(ROL_EPSILON<Real>())),std::invalid_argument,
      ">>> ERROR (ROL::MixedQuantileQuadrangle): Coefficients do not sum to one!");
    TEUCHOS_TEST_FOR_EXCEPTION(plusFunction_ == Teuchos::null, std::invalid_argument,
      ">>> ERROR (ROL::MixedQuantileQuadrangle): PlusFunction pointer is null!");
  }

  void initialize(void) {
    size_ = prob_.size();
    // Initialize temporary storage
    Real zero(0);
    xvar_.clear(); xvar_.resize(size_,zero);
    vvar_.clear(); vvar_.resize(size_,zero);
    vec_.clear();  vec_.resize(size_,zero);
  }

public:

  MixedQuantileQuadrangle( Teuchos::ParameterList &parlist )
    : RiskMeasure<Real>(), firstReset_(true) {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mixed-Quantile Quadrangle");
    // Grab probability and coefficient arrays
    prob_  = Teuchos::getArrayFromStringParameter<Real>(list,"Probability Array");
    coeff_ = Teuchos::getArrayFromStringParameter<Real>(list,"Coefficient Array");
    plusFunction_ = Teuchos::rcp(new PlusFunction<Real>(list));
    // Check inputs
    checkInputs();
    initialize();
  }

  MixedQuantileQuadrangle(const std::vector<Real> &prob,
                          const std::vector<Real> &coeff,
                          const Teuchos::RCP<PlusFunction<Real> > &pf )
    : RiskMeasure<Real>(), plusFunction_(pf), prob_(prob), coeff_(coeff), firstReset_(true) {
    checkInputs();
    initialize();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    Real zero(0);
    RiskMeasure<Real>::reset(x0,x);
    for (int i = 0; i < size_; i++) {
      xvar_[i] = Teuchos::dyn_cast<const RiskVector<Real> >(
                 Teuchos::dyn_cast<const Vector<Real> >(x)).getStatistic(i);
      vec_[i]  = zero;
    }
    if ( firstReset_ ) {
      dualVector_ = (x0->dual()).clone();
      firstReset_ = false;
    }
    dualVector_->zero();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const RiskVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(v)).getVector());
    for (int i = 0; i < size_; i++) {
      vvar_[i] = Teuchos::dyn_cast<const RiskVector<Real> >(
                 Teuchos::dyn_cast<const Vector<Real> >(v)).getStatistic(i);
    }
  }

  void update(const Real val, const Real weight) {
    Real pf(0), one(1);
    for (int i = 0; i < size_; i++) {
      pf = plusFunction_->evaluate(val-xvar_[i],0);
      RiskMeasure<Real>::val_ += weight*coeff_[i]/(one-prob_[i])*pf;
    }
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real pf(0), c(0), one(1);
    for (int i = 0; i < size_; i++) {
      pf = plusFunction_->evaluate(val-xvar_[i],1);
      c  = weight*coeff_[i]/(one-prob_[i])*pf;
      vec_[i] -= c;
      RiskMeasure<Real>::g_->axpy(c,g);
    }
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    Real pf1(0), pf2(0), c(0), one(1);
    for (int i = 0; i < size_; i++) {
      pf1 = plusFunction_->evaluate(val-xvar_[i],1);
      pf2 = plusFunction_->evaluate(val-xvar_[i],2);
      c   = weight*coeff_[i]/(one-prob_[i])*pf2*(gv-vvar_[i]);
      vec_[i] -= c;
      //c  *= (gv-vvar_[i]);
      RiskMeasure<Real>::hv_->axpy(c,g);
      c = weight*coeff_[i]/(one-prob_[i])*pf1;
      RiskMeasure<Real>::hv_->axpy(c,hv);
    }
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val  = RiskMeasure<Real>::val_, cvar(0);
    sampler.sumAll(&val,&cvar,1);
    for (int i = 0; i < size_; i++) {
      cvar += coeff_[i]*xvar_[i];
    }
    return cvar;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    RiskVector<Real> &gs = Teuchos::dyn_cast<RiskVector<Real> >(Teuchos::dyn_cast<Vector<Real> >(g));
    std::vector<Real> var(size_);
    sampler.sumAll(&vec_[0],&var[0],size_);
    
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector_);
    for (int i = 0; i < size_; i++) {
      var[i] += coeff_[i];
    }
    gs.setStatistic(var);
    gs.setVector(*(Teuchos::rcp_dynamic_cast<Vector<Real> >(dualVector_))); 
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    RiskVector<Real> &hs = Teuchos::dyn_cast<RiskVector<Real> >(Teuchos::dyn_cast<Vector<Real> >(hv));
    std::vector<Real> var(size_);
    sampler.sumAll(&vec_[0],&var[0],size_);

    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector_);
//    for (int i = 0; i < size_; i++) {
//      var[i] *= coeff_[i]/(1.0-prob_[i]);
//    }
    hs.setStatistic(var);
    hs.setVector(*(Teuchos::rcp_dynamic_cast<Vector<Real> >(dualVector_)));
  }
};

}

#endif
