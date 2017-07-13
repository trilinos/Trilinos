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

#ifndef ROL_MEANDEVIATIONFROMTARGET_HPP
#define ROL_MEANDEVIATIONFROMTARGET_HPP

#include "ROL_RiskMeasure.hpp"
#include "ROL_PositiveFunction.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_AbsoluteValue.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

/** @ingroup risk_group
    \class ROL::MeanDeviationFromTarget
    \brief Provides an interface for the mean plus a sum of arbitrary order
    deviations from targets.

    The mean plus deviations from targets risk measure is
    \f[
       \mathcal{R}(X) = \mathbb{E}[X]
        + \sum_{k=1}^n c_k \mathbb{E}[\wp(X-t_k)^{p_k}]^{1/p_k}
    \f]
    where \f$\wp:\mathbb{R}\to[0,\infty)\f$ is either the absolute value
    or \f$(x)_+ = \max\{0,x\}\f$, \f$c_k > 0\f$ and \f$p_k\in\mathbb{N}\f$.
    In general, \f$\mathcal{R}\f$ is law-invariant, but not coherent since
    it violates translation equivariance.

    When using derivative-based optimization, the user can
    provide a smooth approximation of \f$(\cdot)_+\f$ using the
    ROL::PositiveFunction class.
*/

namespace ROL {

template<class Real>
class MeanDeviationFromTarget : public RiskMeasure<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  Teuchos::RCP<PositiveFunction<Real> > positiveFunction_;

  Teuchos::RCP<Vector<Real> > dualVector1_;
  Teuchos::RCP<Vector<Real> > dualVector2_;
  Teuchos::RCP<Vector<Real> > dualVector3_;
  Teuchos::RCP<Vector<Real> > dualVector4_;

  std::vector<Real> target_;
  std::vector<Real> order_;
  std::vector<Real> coeff_;
  uint NumMoments_;

  std::vector<Real> pval_; 
  std::vector<Real> pgv_; 

  std::vector<Teuchos::RCP<Vector<Real> > > pg0_;
  std::vector<Teuchos::RCP<Vector<Real> > > pg_;
  std::vector<Teuchos::RCP<Vector<Real> > > phv_;

  bool firstReset_;

  void initialize(void) {
    // Initialize additional storage
    pg_.clear(); pg0_.clear(); phv_.clear(); pval_.clear(); pgv_.clear();
    pg_.resize(NumMoments_);
    pg0_.resize(NumMoments_);
    phv_.resize(NumMoments_);
    pval_.resize(NumMoments_);
    pgv_.resize(NumMoments_);
  }

  void checkInputs(void) const {
    int oSize = order_.size(), cSize = coeff_.size(), tSize = target_.size();
    TEUCHOS_TEST_FOR_EXCEPTION((oSize!=cSize),std::invalid_argument,
      ">>> ERROR (ROL::MeanDeviationFromTarget): Order and coefficient arrays have different sizes!");
    TEUCHOS_TEST_FOR_EXCEPTION((oSize!=tSize),std::invalid_argument,
      ">>> ERROR (ROL::MeanDeviationFromTarget): Order and target arrays have different sizes!");
    Real zero(0), two(2);
    for (int i = 0; i < oSize; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION((order_[i] < two), std::invalid_argument,
        ">>> ERROR (ROL::MeanDeviationFromTarget): Element of order array out of range!");
      TEUCHOS_TEST_FOR_EXCEPTION((coeff_[i] < zero), std::invalid_argument,
        ">>> ERROR (ROL::MeanDeviationFromTarget): Element of coefficient array out of range!");
    }
    TEUCHOS_TEST_FOR_EXCEPTION(positiveFunction_ == Teuchos::null, std::invalid_argument,
      ">>> ERROR (ROL::MeanDeviationFromTarget): PositiveFunction pointer is null!");
  }

public:
  /** \brief Constructor.

      @param[in]     target  is the target scalar value
      @param[in]     order   is the deviation order
      @param[in]     coeff   is the weight for deviation term
      @param[in]     pf      is the plus function or an approximation

      This constructor produces a mean plus deviation from target risk measure
      with a single deviation.
  */
  MeanDeviationFromTarget( const Real target, const Real order, const Real coeff,
                           const Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf), firstReset_(true) {
    order_.clear();  order_.push_back(order);
    coeff_.clear();  coeff_.push_back(coeff);
    target_.clear(); target_.push_back(target);
    checkInputs();
    NumMoments_ = order_.size();
    initialize();
  }

  /** \brief Constructor.

      @param[in]     target  is a vector of targets
      @param[in]     order   is a vector of deviation orders
      @param[in]     coeff   is a vector of weights for the deviation terms
      @param[in]     pf      is the plus function or an approximation

      This constructor produces a mean plus deviation from target risk measure
      with an arbitrary number of deviations.
  */
  MeanDeviationFromTarget( const std::vector<Real> &target,
                           const std::vector<Real> &order,
                           const std::vector<Real> &coeff, 
                           const Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf), firstReset_(true) {
    target_.clear(); order_.clear(); coeff_.clear();
    for ( uint i = 0; i < target.size(); i++ ) {
      target_.push_back(target[i]);
    }
    for ( uint i = 0; i < order.size(); i++ ) {
      order_.push_back(order[i]);
    }
    for ( uint i = 0; i < coeff.size(); i++ ) {
      coeff_.push_back(coeff[i]);
    }
    checkInputs();
    NumMoments_ = order_.size();
    initialize();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"Mean Plus Deviation From Target" and
      within the "Mean Plus Deviation From Target" sublist should have the following parameters
      \li "Targets" (array of scalars)
      \li "Orders" (array of unsigned integers)
      \li "Coefficients" (array of positive scalars)
      \li "Deviation Type" (eighter "Upper" or "Absolute")
      \li A sublist for positive function information.
  */
  MeanDeviationFromTarget( Teuchos::ParameterList &parlist )
    : RiskMeasure<Real>(), firstReset_(true) {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Deviation From Target");
    // Get data from parameter list
    Teuchos::Array<Real> target
      = Teuchos::getArrayFromStringParameter<double>(list,"Targets");
    target_ = target.toVector();
    Teuchos::Array<Real> order
      = Teuchos::getArrayFromStringParameter<double>(list,"Orders");
    order_ = order.toVector();
    Teuchos::Array<Real> coeff
      = Teuchos::getArrayFromStringParameter<double>(list,"Coefficients");
    coeff_ = coeff.toVector();
    // Build (approximate) positive function
    std::string type = list.get<std::string>("Deviation Type");
    if ( type == "Upper" ) {
      positiveFunction_ = Teuchos::rcp(new PlusFunction<Real>(list));
    }
    else if ( type == "Absolute" ) {
      positiveFunction_ = Teuchos::rcp(new AbsoluteValue<Real>(list));
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> (ROL::MeanDeviation): Deviation type is not recoginized!");
    }
    // Check inputs
    checkInputs();
    NumMoments_ = order.size();
    initialize();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    Real zero(0);
    RiskMeasure<Real>::reset(x0,x);
    if (firstReset_) {
      for ( uint p = 0; p < NumMoments_; p++ ) {
        pg0_[p] = (x0->dual()).clone();
        pg_[p]  = (x0->dual()).clone();
        phv_[p] = (x0->dual()).clone();
      }
      dualVector1_ = (x0->dual()).clone();
      dualVector2_ = (x0->dual()).clone();
      dualVector3_ = (x0->dual()).clone();
      dualVector4_ = (x0->dual()).clone();
      firstReset_  = false;
    }
    for ( uint p = 0; p < NumMoments_; p++ ) {
      pg0_[p]->zero(); pg_[p]->zero(); phv_[p]->zero();
      pval_[p] = zero; pgv_[p] = zero;
    }
    dualVector1_->zero(); dualVector2_->zero();
    dualVector3_->zero(); dualVector4_->zero();
  }
    
  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const RiskVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(v)).getVector());
  }
  
  void update(const Real val, const Real weight) {
    Real diff(0), pf0(0);
    RiskMeasure<Real>::val_ += weight * val;
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val-target_[p];
      pf0  = positiveFunction_->evaluate(diff,0);
      pval_[p] += weight * std::pow(pf0,order_[p]);
    }
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real diff(0), pf0(0), pf1(0), c(0), one(1);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val-target_[p];
      pf0 = positiveFunction_->evaluate(diff,0);
      pf1 = positiveFunction_->evaluate(diff,1);
      c    = std::pow(pf0,order_[p]-one) * pf1;
      (pg_[p])->axpy(weight * c,g);
      pval_[p] += weight * std::pow(pf0,order_[p]);
    }
    RiskMeasure<Real>::g_->axpy(weight,g);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    Real diff(0), pf0(0), pf1(0), pf2(0), p0(0), p1(0), p2(0), c(0), one(1), two(2);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val - target_[p];
      pf0 = positiveFunction_->evaluate(diff,0);
      pf1 = positiveFunction_->evaluate(diff,1);
      pf2 = positiveFunction_->evaluate(diff,2);
      p0   = std::pow(pf0,order_[p]);
      p1   = std::pow(pf0,order_[p]-one);
      p2   = std::pow(pf0,order_[p]-two);
      c    = -(order_[p]-one)*p1*pf1;
      pg0_[p]->axpy(weight*c,g);
      c    = gv*((order_[p]-one)*p2*pf1*pf1 + p1*pf2);
      pg_[p]->axpy(weight*c,g);
      c    = p1*pf1;
      phv_[p]->axpy(weight*c,hv);
      pval_[p] += weight*p0;
      pgv_[p]  += weight*p1*pf1*gv;
    }
    RiskMeasure<Real>::hv_->axpy(weight,hv);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val = RiskMeasure<Real>::val_, dev(0), one(1);
    sampler.sumAll(&val,&dev,1);
    std::vector<Real> pval_sum(NumMoments_);
    sampler.sumAll(&(pval_)[0],&pval_sum[0],NumMoments_);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      dev += coeff_[p] * std::pow(pval_sum[p],one/order_[p]);
    }
    return dev;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    Real zero(0), one(1);
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector1_);
    std::vector<Real> pval_sum(NumMoments_);
    sampler.sumAll(&(pval_)[0],&pval_sum[0],NumMoments_);
    Teuchos::RCP<Vector<Real> > pg;
    for ( uint p = 0; p < NumMoments_; p++ ) {
      if ( pval_sum[p] > zero ) {
        pg = (pg_[p])->clone();
        sampler.sumAll(*(pg_[p]),*pg);
        dualVector1_->axpy(coeff_[p]/std::pow(pval_sum[p],one-one/order_[p]),*pg);
      }
    }
    // Set RiskVector
    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setVector(*dualVector1_);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    Real zero(0), one(1), two(2);
    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector1_);
    std::vector<Real> pval_sum(NumMoments_);
    sampler.sumAll(&(pval_)[0],&pval_sum[0],NumMoments_);
    std::vector<Real> pgv_sum(NumMoments_);
    sampler.sumAll(&(pgv_)[0],&pgv_sum[0],NumMoments_);
    Real c(0);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      if ( pval_sum[p] > zero ) {
        sampler.sumAll(*(pg_[p]),*dualVector2_);
        sampler.sumAll(*(pg0_[p]),*dualVector3_);
        sampler.sumAll(*(phv_[p]),*dualVector4_);
        c = coeff_[p]*(pgv_sum[p]/std::pow(pval_sum[p],two-one/order_[p]));
        dualVector1_->axpy(c,*dualVector3_);
        c = coeff_[p]/std::pow(pval_sum[p],one-one/order_[p]);
        dualVector1_->axpy(c,*dualVector2_);
        dualVector1_->axpy(c,*dualVector4_);
      }
    }
    // Set RiskVector
    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setVector(*dualVector1_);
  }
};

}

#endif
