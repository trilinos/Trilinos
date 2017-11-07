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

#ifndef ROL_EXPECTATIONQUAD_HPP
#define ROL_EXPECTATIONQUAD_HPP

#include "ROL_RiskVector.hpp"
#include "ROL_RiskMeasure.hpp"
#include "ROL_Types.hpp"

/** @ingroup risk_group
    \class ROL::ExpectationQuad
    \brief Provides a general interface for risk measures generated through
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
class ExpectationQuad : public RiskMeasure<Real> {
private:

  Teuchos::RCP<Vector<Real> > dualVector_;

  Real xstat_;
  Real vstat_;

  bool firstReset_;

public:
  ExpectationQuad(void) : RiskMeasure<Real>(), xstat_(0), vstat_(0), firstReset_(true) {}

  /** \brief Evaluate the scalar regret function at x.

      @param[in]   x      is the scalar input
      @param[in]   deriv  is the derivative order

      This function returns \f$v(x)\f$ or a derivative of \f$v(x)\f$.
  */
  virtual Real regret(Real x, int deriv = 0) = 0;

  /** \brief Run default derivative tests for the scalar regret function.
  */
  virtual void checkRegret(void) {
    Real zero(0), half(0.5), two(2), one(1), oem3(1.e-3), fem4(5.e-4), p1(0.1);
    // Check v(0) = 0
    Real x = zero;
    Real vx = regret(x,0);
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v(0) = 0? \n";
    std::cout << std::right << std::setw(20) << "v(0)" << "\n";
    std::cout << std::scientific << std::setprecision(11) << std::right 
              << std::setw(20) << std::abs(vx) 
              << "\n";
    std::cout << "\n";
    // Check v(x) > x
    Real scale = two;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: x < v(x) for |x| > 0? \n";
    std::cout << std::right << std::setw(20) << "x"
              << std::right << std::setw(20) << "v(x)"
              << "\n";
    for (int i = 0; i < 10; i++) {
      x = scale*(Real)rand()/(Real)RAND_MAX - scale*half;
      vx = regret(x,0);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << x 
                << std::setw(20) << vx 
                << "\n";
      scale *= two;
    }
    std::cout << "\n";
    // Check v(x) is convex
    Real y = zero;
    Real vy = zero;
    Real z = zero;
    Real vz = zero;
    Real l = zero; 
    scale = two;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v(x) is convex? \n";
    std::cout << std::right << std::setw(20) << "v(l*x+(1-l)*y)" 
                            << std::setw(20) << "l*v(x)+(1-l)*v(y)" 
                            << "\n";
    for (int i = 0; i < 10; i++) {
      x = scale*(Real)rand()/(Real)RAND_MAX - scale*half;
      vx = regret(x,0);
      y = scale*(Real)rand()/(Real)RAND_MAX - scale*half;
      vy = regret(y,0);
      l = (Real)rand()/(Real)RAND_MAX;
      z = l*x + (one-l)*y;
      vz = regret(z,0);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << vz 
                << std::setw(20) << l*vx + (one-l)*vy 
                << "\n";
      scale *= two;
    }
    std::cout << "\n";
    // Check v'(x)
    x = oem3*(Real)rand()/(Real)RAND_MAX - fem4;
    vx = regret(x,0);
    Real dv = regret(x,1);
    Real t = one;
    Real diff = zero;
    Real err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v'(x) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v'(x)"
                            << std::setw(20) << "(v(x+t)-v(x))/t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      y = x + t;
      vy = regret(y,0);
      diff = (vy-vx)/t;
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right 
                << std::setw(20) << t
                << std::setw(20) << dv 
                << std::setw(20) << diff 
                << std::setw(20) << err 
                << "\n";
      t *= p1;
    }
    std::cout << "\n";
    // Check v''(x)
    x = oem3*(Real)rand()/(Real)RAND_MAX - fem4;
    vx = regret(x,1);
    dv = regret(x,2);
    t = one;
    diff = zero;
    err = zero;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v''(x) is correct? \n";
    std::cout << std::right << std::setw(20) << "t"
                            << std::setw(20) << "v''(x)"
                            << std::setw(20) << "(v'(x+t)-v'(x))/t"
                            << std::setw(20) << "Error"
                            << "\n";
    for (int i = 0; i < 13; i++) {
      y = x + t;
      vy = regret(y,1);
      diff = (vy-vx)/t;
      err = std::abs(diff-dv);
      std::cout << std::scientific << std::setprecision(11) << std::right 
                << std::setw(20) << t
                << std::setw(20) << dv 
                << std::setw(20) << diff 
                << std::setw(20) << err 
                << "\n";
      t *= p1;
    }
    std::cout << "\n";
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::reset(x0,x);
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    xstat_ = (*Teuchos::dyn_cast<const RiskVector<Real> >(
               Teuchos::dyn_cast<const Vector<Real> >(x)).getStatistic(comp,index))[0];
    if (firstReset_) {
      dualVector_            = (x0->dual()).clone();
      firstReset_ = false;
    }
    dualVector_->zero();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x, 
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const RiskVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(v)).getVector());
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    vstat_ = (*Teuchos::dyn_cast<const RiskVector<Real> >(
               Teuchos::dyn_cast<const Vector<Real> >(v)).getStatistic(comp,index))[0];
  }

  void update(const Real val, const Real weight) {
    Real r = regret(val-xstat_,0);
    RiskMeasure<Real>::val_ += weight * r;
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real r = regret(val-xstat_,1);
    RiskMeasure<Real>::val_ -= weight * r;
    RiskMeasure<Real>::g_->axpy(weight*r,g);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv, 
                      const Real weight) {
    Real r1 = regret(val-xstat_,1);
    Real r2 = regret(val-xstat_,2);
    RiskMeasure<Real>::val_ += weight * r2 * (vstat_ - gv);
    RiskMeasure<Real>::hv_->axpy(weight*r2*(gv-vstat_),g);
    RiskMeasure<Real>::hv_->axpy(weight*r1,hv);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val  = RiskMeasure<Real>::val_, gval(0);
    sampler.sumAll(&val,&gval,1);
    gval += xstat_;
    return gval;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    RiskVector<Real> &gs = Teuchos::dyn_cast<RiskVector<Real> >(Teuchos::dyn_cast<Vector<Real> >(g));
    Real stat  = RiskMeasure<Real>::val_, gstat(0), one(1);
    sampler.sumAll(&stat,&gstat,1);
    gstat += one;
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    gs.setStatistic(gstat,comp,index);

    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector_);
    gs.setVector(*dualVector_);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    RiskVector<Real> &hs = Teuchos::dyn_cast<RiskVector<Real> >(Teuchos::dyn_cast<Vector<Real> >(hv));
    Real stat  = RiskMeasure<Real>::val_, gstat(0);
    sampler.sumAll(&stat,&gstat,1);
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    hs.setStatistic(gstat,comp,index);

    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector_);
    hs.setVector(*dualVector_);
  }
};

}

#endif
