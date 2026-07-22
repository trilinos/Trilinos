// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_EXPECTATIONQUADERROR_HPP
#define ROL_EXPECTATIONQUADERROR_HPP

#include "ROL_ExpectationQuad.hpp"
#include "ROL_RandVarFunctional.hpp"
#include "ROL_Types.hpp"

/** @ingroup risk_group
    \class ROL::ExpectationQuadError
    \brief Provides a general interface for error measures generated through
           the expectation risk quadrangle.

    The expectation risk quadrangle is a specialization of the general
    risk quadrangle that provides a rigorous connection between risk-averse
    optimization and statistical estimation.  The risk quadrangle provides
    fundamental relationships between measures of risk, regret, error and
    deviation.  An expectation risk quadrangle is defined through scalar
    regret and error functions.  The scalar regret function,
    \f$v:\mathbb{R}\to(-\infty,\infty]\f$, must be proper, closed, convex
    and satisfy \f$v(0)=0\f$ and \f$v(x) > x\f$ for all \f$x\neq 0\f$.
    Similarly, the scalar error function,
    \f$e:\mathbb{R}\to[0,\infty]\f$, must be proper, closed, convex
    and satisfy \f$e(0)=0\f$ and \f$e(x) > 0\f$ for all \f$x\neq 0\f$.
    \f$v\f$ and \f$e\f$ are obtained from one another through the relations
    \f[
       v(x) = e(x) + x \quad\text{and}\quad e(x) = v(x) - x.
    \f]
    Given \f$v\f$ (or equivalently \f$e\f$), the associated risk measure
    is
    \f[
       \mathcal{R}(X) = \inf_{t\in\mathbb{R}} \left\{
         t + \mathbb{E}\left[v(X-t)\right]
         \right\}.
    \f]
    In general, \f$\mathcal{R}\f$ is convex and translation equivariant.
    Moreover, \f$\mathcal{R}\f$ is monotonic if \f$v\f$ is increasing
    and \f$\mathcal{R}\f$ is positive homogeneous if \f$v\f$ is.
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$t\f$, then minimizes jointly for \f$(x_0,t)\f$.
*/


namespace ROL {

template<class Real>
class ExpectationQuadError : public RandVarFunctional<Real> {
private:
  Ptr<ExpectationQuad<Real>> eq_;

  using RandVarFunctional<Real>::val_;
  using RandVarFunctional<Real>::gv_;
  using RandVarFunctional<Real>::g_;
  using RandVarFunctional<Real>::hv_;
  using RandVarFunctional<Real>::dualVector_;

  using RandVarFunctional<Real>::point_;
  using RandVarFunctional<Real>::weight_;

  using RandVarFunctional<Real>::computeValue;
  using RandVarFunctional<Real>::computeGradient;
  using RandVarFunctional<Real>::computeGradVec;
  using RandVarFunctional<Real>::computeHessVec;

public:
  ExpectationQuadError(const Ptr<ExpectationQuad<Real>> & eq)
    : RandVarFunctional<Real>(), eq_(eq) {}

  /** \brief Run derivative tests for the scalar error function.
  */
  void checkRegret(void) {
    eq_->check();
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    val_ += weight_ * eq_->error(val,0);
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real r   = eq_->error(val,1);
    if (std::abs(r) >= ROL_EPSILON<Real>()) {
      computeGradient(*dualVector_,obj,x,tol);
      g_->axpy(weight_*r,*dualVector_);
    }
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real r1  = eq_->error(val,1);
    Real r2  = eq_->error(val,2);
    if (std::abs(r2) >= ROL_EPSILON<Real>()) {
      Real gv = computeGradVec(*dualVector_,obj,v,x,tol);
      hv_->axpy(weight_*r2*gv,*dualVector_);
    }
    if (std::abs(r1) >= ROL_EPSILON<Real>()) {
      computeHessVec(*dualVector_,obj,v,x,tol);
      hv_->axpy(weight_*r1,*dualVector_);
    }
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    Real val(0);
    sampler.sumAll(&val_,&val,1);
    return val;
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    sampler.sumAll(*g_,g);
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    sampler.sumAll(*hv_,hv);
  }
};

}

#endif
