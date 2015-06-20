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

#include "ROL_CVaRVector.hpp"
#include "ROL_RiskMeasure.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class ExpectationQuad : public RiskMeasure<Real> {
private:
  Real xstat_;
  Real vstat_;

public:
  ExpectationQuad(void) : xstat_(0.0), vstat_(0.0) {}

  virtual Real regret(Real x, int deriv = 0) = 0;

  virtual void checkRegret(void) {
    // Check v(0) = 0
    Real x = 0.0;
    Real vx = regret(x,0);
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v(0) = 0? \n";
    std::cout << std::right << std::setw(20) << "v(0)" << "\n";
    std::cout << std::scientific << std::setprecision(11) << std::right 
              << std::setw(20) << std::abs(vx) 
              << "\n";
    std::cout << "\n";
    // Check v(x) > x
    Real scale = 2.0;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: x < v(x) for |x| > 0? \n";
    std::cout << std::right << std::setw(20) << "x"
              << std::right << std::setw(20) << "v(x)"
              << "\n";
    for (int i = 0; i < 10; i++) {
      x = scale*(Real)rand()/(Real)RAND_MAX - scale*0.5;
      vx = regret(x,0);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << x 
                << std::setw(20) << vx 
                << "\n";
      scale *= 2.0;
    }
    std::cout << "\n";
    // Check v(x) is convex
    Real y = 0.0;
    Real vy = 0.0;
    Real z = 0.0;
    Real vz = 0.0;
    Real l = 0.0; 
    scale = 2.0;
    std::cout << std::right << std::setw(20) << "CHECK REGRET: v(x) is convex? \n";
    std::cout << std::right << std::setw(20) << "v(l*x+(1-l)*y)" 
                            << std::setw(20) << "l*v(x)+(1-l)*v(y)" 
                            << "\n";
    for (int i = 0; i < 10; i++) {
      x = scale*(Real)rand()/(Real)RAND_MAX - scale*0.5;
      vx = regret(x,0);
      y = scale*(Real)rand()/(Real)RAND_MAX - scale*0.5;
      vy = regret(y,0);
      l = (Real)rand()/(Real)RAND_MAX;
      z = l*x + (1.0-l)*y;
      vz = regret(z,0);
      std::cout << std::scientific << std::setprecision(11) << std::right
                << std::setw(20) << vz 
                << std::setw(20) << l*vx + (1.0-l)*vy 
                << "\n";
      scale *= 2.0;
    }
    std::cout << "\n";
    // Check v'(x)
    x = 0.001*(Real)rand()/(Real)RAND_MAX - 0.0005;
    vx = regret(x,0);
    Real dv = regret(x,1);
    Real t = 1.0;
    Real diff = 0.0;
    Real err = 0.0;
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
      t *= 0.1;
    }
    std::cout << "\n";
    // Check v''(x)
    x = 0.001*(Real)rand()/(Real)RAND_MAX - 0.0005;
    vx = regret(x,1);
    dv = regret(x,2);
    t = 1.0;
    diff = 0.0;
    err = 0.0;
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
      t *= 0.1;
    }
    std::cout << "\n";
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    x0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const CVaRVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(x)).getVector());
    xstat_ = Teuchos::dyn_cast<const CVaRVector<Real> >(
               Teuchos::dyn_cast<const Vector<Real> >(x)).getVaR();
    RiskMeasure<Real>::val_ = 0.0;
    RiskMeasure<Real>::g_  = x0->clone(); RiskMeasure<Real>::g_->zero();
    RiskMeasure<Real>::hv_ = x0->clone(); RiskMeasure<Real>::hv_->zero();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x, 
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    x0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const CVaRVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(x)).getVector());
    xstat_ = Teuchos::dyn_cast<const CVaRVector<Real> >(
               Teuchos::dyn_cast<const Vector<Real> >(x)).getVaR();
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const CVaRVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(v)).getVector());
    vstat_ = Teuchos::dyn_cast<const CVaRVector<Real> >(
               Teuchos::dyn_cast<const Vector<Real> >(v)).getVaR();
    RiskMeasure<Real>::val_ = 0.0;
    RiskMeasure<Real>::g_  = x0->clone(); RiskMeasure<Real>::g_->zero();
    RiskMeasure<Real>::hv_ = x0->clone(); RiskMeasure<Real>::hv_->zero();
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
    Real val  = RiskMeasure<Real>::val_;
    Real gval = 0.0;
    sampler.sumAll(&val,&gval,1);
    gval += xstat_;
    return gval;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    CVaRVector<Real> &gs = Teuchos::dyn_cast<CVaRVector<Real> >(Teuchos::dyn_cast<Vector<Real> >(g));
    Real stat  = RiskMeasure<Real>::val_;
    Real gstat = 0.0;
    sampler.sumAll(&stat,&gstat,1);
    gstat += 1.0;
    gs.setVaR(gstat);

    Teuchos::RCP<Vector<Real> > gz = RiskMeasure<Real>::g_->clone();
    sampler.sumAll(*(RiskMeasure<Real>::g_),*gz);
    gs.setVector(*gz);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    CVaRVector<Real> &hs = Teuchos::dyn_cast<CVaRVector<Real> >(Teuchos::dyn_cast<Vector<Real> >(hv));
    Real stat  = RiskMeasure<Real>::val_;
    Real gstat = 0.0;
    sampler.sumAll(&stat,&gstat,1);
    hs.setVaR(gstat);

    Teuchos::RCP<Vector<Real> > hz = RiskMeasure<Real>::hv_->clone();
    sampler.sumAll(*(RiskMeasure<Real>::hv_),*hz);
    hs.setVector(*hz);
  }
};

}

#endif
