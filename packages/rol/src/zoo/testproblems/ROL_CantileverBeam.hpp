// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
 *  \brief Contains definitions for Cantilevered Beam example in
 *         G. N. Vanderplaats, Numerical Optimization Techniques for Engineering
 *         Design: With Applications (1984).  This problem contains bound and
 *         inequality constraints.
 */

#ifndef ROL_CANTILEVERBEAM_HPP
#define ROL_CANTILEVERBEAM_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"

namespace ROL {
namespace ZOO {

template<class Real>
class Objective_CantileverBeam : public StdObjective<Real> {
private:
  int nseg_;
  Real len_;
  std::vector<Real> l_;

public:
  Objective_CantileverBeam(const int nseg = 5)
    : nseg_(nseg), len_(500) {
    l_.clear(); l_.resize(nseg, len_/static_cast<Real>(nseg_));
  }

  Real value( const std::vector<Real> &x, Real &tol ) {
    Real val(0);
    for (int i = 0; i < nseg_; ++i) {
      val += l_[i]*x[i]*x[i+nseg_];
    }
    return val;
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    for (int i = 0; i < nseg_; ++i) {
      g[i]       = l_[i]*x[i+nseg_];
      g[i+nseg_] = l_[i]*x[i];
    }
  }  

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    for (int i = 0; i < nseg_; ++i) {
      hv[i]       = l_[i]*v[i+nseg_];
      hv[i+nseg_] = l_[i]*v[i];
    }
  }
}; // class Objective_CantileverBeam


template<class Real> 
class Constraint_CantileverBeam : public StdConstraint<Real>  {
private:
  int nseg_;
  Real len_, P_, E_, Sigma_max_, tip_max_;
  std::vector<Real> l_, suml_, M_, dyp_;

public:
  Constraint_CantileverBeam(const int nseg = 5)
    : nseg_(nseg), len_(500), P_(50e3), E_(200e5), Sigma_max_(14e3), tip_max_(2.5), suml_(0) {
    const Real half(0.5);
    l_.clear();
    l_.resize(nseg_, len_/static_cast<Real>(nseg_));
    suml_.resize(nseg_);
    M_.resize(nseg_);
    dyp_.resize(nseg_);
    for (int i = 0; i < nseg_; ++i) {
      suml_[i] = static_cast<Real>(i+1)*l_[i];
      M_[i]    = P_*(len_ - suml_[i] + l_[i]);
      dyp_[i]  = (len_ - suml_[i] + half*l_[i])*static_cast<Real>(nseg_-i-1);
    }
  }

  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real one(1), two(2), three(3), twelve(12), twenty(20);
    std::vector<Real> y1(nseg_), yp(nseg_);
    Real Inertia(0), sigma(0), sumy1(0), sumypl(0);
    for (int i = 0; i < nseg_; ++i) {
      // stress constraint
      Inertia    = x[i]*std::pow(x[i+nseg_],3)/twelve;
      sigma      = M_[i]*(x[i+nseg_]/two)/Inertia;
      c[i]       = sigma / Sigma_max_ - one;
      // lateral slope constraint
      c[i+nseg_] = x[i+nseg_] - twenty*x[i];
      // tip deflection constraint
      y1[i]      = P_*std::pow(l_[i],2)/(two*E_*Inertia) * (len_ - suml_[i] + two/three*l_[i]);
      yp[i]      = P_*l_[i]/(E_*Inertia) * (len_ - suml_[i] + l_[i]/two);
      sumy1     += y1[i];
    }
    // tip deflection constraint
    for (int i = 1; i < nseg_; ++i) {
      yp[i]  += yp[i-1];
      sumypl += yp[i-1]*l_[i];
    }
    Real yN = sumy1 + sumypl;
    c[2*nseg_] = yN / tip_max_ - one;
  }

  void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, 
                      const std::vector<Real> &x, Real &tol ) {
    const Real two(2), three(3), six(6), twelve(12), twenty(20);
    // const Real one(1); // Unused
    Real Inertia(0), dyN(0), sumW(0), sumH(0);
    for (int i = 0; i < nseg_; ++i) {
      // stress constraint
      jv[i]       = -six/Sigma_max_*M_[i]*v[i]/std::pow(x[i]*x[i+nseg_],2)
                    -twelve/Sigma_max_*M_[i]*v[i+nseg_]/(x[i]*std::pow(x[i+nseg_],3));
      // lateral slope constraint
      jv[i+nseg_] = v[i+nseg_] - twenty*v[i];
      // tip deflection constraint
      Inertia     = x[i]*std::pow(x[i+nseg_],3)/twelve;
      dyN         = -P_/E_ * std::pow(l_[i],2)/Inertia*((len_ - suml_[i] + two/three*l_[i])/two + dyp_[i]);
      sumW       += v[i]*(dyN/x[i])/tip_max_;
      sumH       += three*v[i+nseg_]*(dyN/x[i+nseg_])/tip_max_;
    }
    // tip deflection constraint
    jv[2*nseg_] = sumW+sumH;
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, 
                             const std::vector<Real> &x, Real &tol ) {
    const Real two(2), three(3), six(6), twelve(12), twenty(20);
//    const Real one(1);  // Unused
    Real Inertia(0), dyN(0);
    for (int i = 0; i < nseg_; ++i) {
      // stress constraint
      ajv[i]        = -six/Sigma_max_*M_[i]*v[i]/std::pow(x[i]*x[i+nseg_],2);
      ajv[i+nseg_]  = -twelve/Sigma_max_*M_[i]*v[i]/(x[i]*std::pow(x[i+nseg_],3));
      // lateral slope constraint
      ajv[i]       += -twenty*v[i+nseg_];
      ajv[i+nseg_] += v[i+nseg_];
      // tip deflection constraint
      Inertia       = x[i]*std::pow(x[i+nseg_],3)/twelve;
      dyN           = -P_/E_ * std::pow(l_[i],2)/Inertia*((len_ - suml_[i] + two/three*l_[i])/two + dyp_[i]);
      ajv[i]       += v[2*nseg_]*(dyN/x[i])/tip_max_;
      ajv[i+nseg_] += three*v[2*nseg_]*(dyN/x[i+nseg_])/tip_max_;
    }
  }

//  void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u,
//                            const std::vector<Real> &v, const std::vector<Real> &x,
//                            Real &tol ) {
//  }

}; // class Constraint_CantileverBeam



template<class Real>
class getCantileverBeam : public TestProblem<Real> {
private:
  int nseg_;

public:
  getCantileverBeam(int nseg = 5) : nseg_(nseg) {}

  Ptr<Objective<Real> > getObjective( void ) const {
    return makePtr<Objective_CantileverBeam<Real>>(nseg_);
  }

  Ptr<Constraint<Real> > getInequalityConstraint( void ) const {
    return makePtr<Constraint_CantileverBeam<Real>>(nseg_);
  }

  Ptr<BoundConstraint<Real>> getBoundConstraint( void ) const {
    // Lower bound is zero  
    Ptr<std::vector<Real>> lp = makePtr<std::vector<Real>>(2*nseg_);
    for (int i = 0; i < nseg_; ++i) {
      (*lp)[i]       = static_cast<Real>(1);
      (*lp)[i+nseg_] = static_cast<Real>(5);
    }
    
    // No upper bound
    Ptr<std::vector<Real>> up = makePtr<std::vector<Real>>(2*nseg_,ROL_INF<Real>());
   
    Ptr<Vector<Real>> l = makePtr<StdVector<Real>>(lp);
    Ptr<Vector<Real>> u = makePtr<StdVector<Real>>(up);
  
    return makePtr<Bounds<Real>>(l,u);
  }

  Ptr<Vector<Real>> getInitialGuess( void ) const {
    Ptr<std::vector<Real> > x0p = makePtr<std::vector<Real>>(2*nseg_,0);
    for (int i = 0; i < nseg_; ++i) {
      (*x0p)[i]       = static_cast<Real>(5);
      (*x0p)[i+nseg_] = static_cast<Real>(40);
    }
  
    return makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real>> getSolution( const int i = 0 ) const {
    //Ptr<std::vector<Real> > xp = makePtr<std::vector<Real>>(2);
    //(*xp)[0] = 3.0;
    //(*xp)[1] = std::sqrt(3.0);

    Ptr<std::vector<Real>> xp = makePtr<std::vector<Real>>(2*nseg_);
    for (int i = 0; i < nseg_; ++i) {
      (*xp)[i]       = 3.0;
      (*xp)[i+nseg_] = std::sqrt(3.0);
    }

    return makePtr<StdVector<Real>>(xp);
  }

  Ptr<Vector<Real>> getInequalityMultiplier( void ) const {
    Ptr<std::vector<Real> > lp = makePtr<std::vector<Real>>(2*nseg_+1,0.0);
    return makePtr<StdVector<Real>>(lp);
  }

  Ptr<BoundConstraint<Real>> getSlackBoundConstraint(void) const {
    // No lower bound
    Ptr<std::vector<Real> > lp = makePtr<std::vector<Real>>(2*nseg_+1,ROL_NINF<Real>());
    
    // Upper bound is zero
    Ptr<std::vector<Real> > up = makePtr<std::vector<Real>>(2*nseg_+1,0);
   
    Ptr<Vector<Real> > l = makePtr<StdVector<Real>>(lp);
    Ptr<Vector<Real> > u = makePtr<StdVector<Real>>(up);
  
    return makePtr<Bounds<Real>>(l,u);
  }
};


} // namespace ZOO
} // namespace ROL

#endif // ROL_CANTILEVERBEAM_HPP

