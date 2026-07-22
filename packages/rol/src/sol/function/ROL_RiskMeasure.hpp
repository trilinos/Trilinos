// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RISKMEASURE_HPP
#define ROL_RISKMEASURE_HPP

#include "ROL_RiskVector.hpp"

/** @ingroup stochastic_group 
    \class ROL::RiskMeasure
    \brief Provides the interface to implement risk measures.

    Let \f$(\Omega,\mathcal{F},\mathbb{P})\f$ be a complete space.
    Here, \f$\Omega\f$ is the set of outcomes,
    \f$\mathcal{F}\subseteq 2^\Omega\f$ is a \f$\sigma\f$-algebra of events and
    \f$\mathbb{P}:\mathcal{F}\to[0,1]\f$ is a probability measure.  Moreover,
    let \f$\mathcal{X}\f$ be a class of random variables.
    A risk measure is an extended real-valued functional that associates
    numerical values to random variables, i.e.,
    \f$\mathcal{R}:\mathcal{X}\to\mathbb{R}\cup\{+\infty\}\f$.  In most cases,
    \f$\mathcal{X} = L^p(\Omega,\mathcal{F},\mathbb{P})\f$ for some
    \f$1\le p\le \infty\f$.

    There are a number of classifications for risk measures.  One important
    class are the coherent risk measures.  \f$\mathcal{R}\f$ is coherent
    if it satisfies the following four axioms:
    \li Convexity: \f$\mathcal{R}(\lambda X + (1-\lambda)X')
                            \le \lambda \mathcal{R}(X)
                               + (1-\lambda)\mathcal{R}(X')\f$
        for all \f$0 \le \lambda \le 1\f$;
    \li Translation Equivariance: \f$\mathcal{R}(X+C)
                                  = \mathcal{R}(X) + C\f$
        for all constants \f$C\f$;
    \li Monotonicity: If \f$X \le X'\f$ \f$\mathbb{P}\f$ almost everywhere
        then \f$\mathcal{R}(X) \le \mathcal{R}(X')\f$;
    \li Positive Homogeneity: \f$\mathcal{R}(tX) = t \mathcal{R}(X)\f$
        for all \f$t \ge 0\f$.

    Another useful characterization is law invariance.  \f$\mathcal{R}\f$ is
    law invariant if \f$\mathcal{R}(X) = \mathcal{R}(X')\f$ whenever
    \f$\mathbb{P}(X\le t) = \mathbb{P}(X'\le t)\f$ for all \f$t\in\mathbb{R}\f$.
    Law invariant risk measures are only functions of the distribution of the
    input random variable.

    ROL's risk measure base class is written in a way to exploit parallel
    sampling.  General risk measures may depend on global information such as
    the expected value of a random variable, \f$\mathbb{E}[X]\f$.  Thus,
    ROL::RiskMeasure contains functions to update intermediate information
    and to compute desired quantities such as risk values, gradients and
    Hessians applied to vectors.
*/

namespace ROL {

template<class Real>
class RiskMeasure {
protected:
  Real val_;
  Real gv_;
  ROL::Ptr<Vector<Real> > g_;
  ROL::Ptr<Vector<Real> > hv_;
  ROL::Ptr<Vector<Real> > dualVector_;
  bool firstReset_;

  int comp_;
  int index_;

public:
  virtual ~RiskMeasure() {}

  RiskMeasure(void) : val_(0), gv_(0), firstReset_(true),
                      comp_(0), index_(0) {}

  void setRiskVectorInfo(const int comp, const int index) {
    comp_ = comp;
    index_ = index;
  }

  int getComponent(void) const {
    return comp_;
  }

  int getIndex(void) const {
    return index_;
  }

  /** \brief Reset internal risk measure storage.
             Called for value and gradient computation.

             @param[out]  x0 is a user-provided optimization vector
             @param[in]   x  is a (potentially) augmented risk vector

             On input, \f$x\f$ carries \f$x_0\f$ and any statistics (scalars)
             associated with the risk measure. 
  */
  virtual void reset(ROL::Ptr<Vector<Real> > &x0, const Vector<Real> &x) {
    x0 = ROL::constPtrCast<Vector<Real> >(
         dynamic_cast<const RiskVector<Real>&>(x).getVector());
    // Create memory for class members
    if ( firstReset_ ) {
      g_  = (x0->dual()).clone();
      hv_ = (x0->dual()).clone();
      dualVector_ = (x0->dual()).clone();
      firstReset_ = false;
    }
    // Zero member variables
    const Real zero(0);
    val_ = zero; gv_ = zero;
    g_->zero(); hv_->zero(); dualVector_->zero();
  }

  /** \brief Reset internal risk measure storage.
             Called for Hessian-times-a-vector computation.

             @param[out]  x0 is a user-provided optimization vector
             @param[in]   x  is a (potentially) augmented risk vector
             @param[out]  v0 is a user-provided direction vector
             @param[in]   v  is a (potentially) augmented risk vector

             On input, \f$x\f$ carries \f$x_0\f$ and any statistics (scalars)
             associated with the risk measure.  Similarly, \f$v\f$ carries
             \f$v_0\f$ and any statistics (scalars) associated with the risk
             measure.
  */
  virtual void reset(ROL::Ptr<Vector<Real> > &x0, const Vector<Real> &x,
                     ROL::Ptr<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    // Get vector component of v.  This is important for CVaR.
    v0 = ROL::constPtrCast<Vector<Real> >(
         dynamic_cast<const RiskVector<Real>&>(v).getVector());
  }

  /** \brief Update internal risk measure storage for value computation.

      @param[in]    val      is the value of the random variable objective
                             function at the current sample point
      @param[in]    weight   is the weight associated with the current
                             sample point
  */
  virtual void update(const Real val, const Real weight) {
    val_ += weight * val;
  }

  /** \brief Update internal risk measure storage for gradient computation.

      @param[in]    val      is the value of the random variable objective
                             function at the current sample point
      @param[in]    g        is the gradient of the random variable objective
                             function at the current sample point
      @param[in]    weight   is the weight associated with the current
                             sample point
  */
  virtual void update(const Real val, const Vector<Real> &g, const Real weight) {
    g_->axpy(weight,g);
  }

  /** \brief Update internal risk measure storage for Hessian-time-a-vector computation.

      @param[in]    val      is the value of the random variable objective
                             function at the current sample point
      @param[in]    g        is the gradient of the random variable objective
                             function at the current sample point
      @param[in]    gv       is the gradient of the random variable objective
                             function at the current sample point applied to
                             the vector v0
      @param[in]    hv       is the Hessian of the random variable objective
                             function at the current sample point applied to
                             the vector v0
      @param[in]    weight   is the weight associated with the current
                             sample point
  */
  virtual void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
                      const Real weight) {
    hv_->axpy(weight,hv);
  }

  /** \brief Return risk measure value.

      @param[in]    sampler is the ROL::SampleGenerator used to sample
                    the objective function

      Upon return, getValue returns \f$\mathcal{R}(f(x_0))\f$ where \f$f(x_0)\f$
      denotes the random variable objective function evaluated at \f$x_0\f$.
  */
  virtual Real getValue(SampleGenerator<Real> &sampler) {
    Real val(0);
    sampler.sumAll(&val_,&val,1);
    return val;
  }

  /** \brief Return risk measure (sub)gradient.

      @param[out]   g       is the (sub)gradient of the risk measure
      @param[in]    sampler is the ROL::SampleGenerator used to sample
                    the objective function

      Upon return, getGradient returns \f$\theta\in\partial\mathcal{R}(f(x_0))\f$
      where \f$f(x_0)\f$ denotes the random variable objective function evaluated
      at \f$x_0\f$ and \f$\partial\mathcal{R}(X)\f$ denotes the subdifferential
      of \f$\mathcal{R}\f$ at \f$X\f$.
  */
  virtual void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    sampler.sumAll(*g_,*dualVector_);
    (dynamic_cast<RiskVector<Real>&>(g)).setVector(*dualVector_);
  }

  /** \brief Return risk measure Hessian-times-a-vector.

      @param[out]   hv      is the Hessian-times-a-vector of the risk measure
      @param[in]    sampler is the ROL::SampleGenerator used to sample
                    the objective function

      Upon return, getHessVec returns \f$\nabla^2 \mathcal{R}(f(x_0))v_0\f$
      (if available)
      where \f$f(x_0)\f$ denotes the random variable objective function evaluated
      at \f$x_0\f$.
  */
  virtual void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    sampler.sumAll(*hv_,*dualVector_);
    (dynamic_cast<RiskVector<Real>&>(hv)).setVector(*dualVector_);
  }
};

}

#endif
