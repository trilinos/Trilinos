// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for W. Hock and K. Schittkowski 25th test function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_HS25_HPP
#define ROL_HS25_HPP

#include "ROL_ScaledStdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace ZOO {

/** \brief W. Hock and K. Schittkowski 25th test function.
 */
template<class Real>
class Objective_HS25 : public Objective<Real> {

typedef typename std::vector<Real>::size_type uint;

private:
  std::vector<Real> u_vec_;
  uint u_size_;

public:
  Objective_HS25() {
    u_size_ = 99;
    for ( uint i = 0; i < u_size_; i++ ) {
      u_vec_.push_back(static_cast<Real>(25)
        + std::pow((static_cast<Real>(-50)
         *std::log(static_cast<Real>(0.01)*static_cast<Real>(i+1))),
          static_cast<Real>(2)/static_cast<Real>(3)));
    }
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    Ptr<const std::vector<Real> > ex
      = dynamic_cast<const PrimalScaledStdVector<Real>&>(x).getVector();

    Real val(0), f(0), u(0);
    Real x1 = (*ex)[0], x2 = (*ex)[1], x3 = (*ex)[2];
    for ( uint i = 0; i < u_size_; i++ ) {
      u = u_vec_[i];
      f = -static_cast<Real>(0.01)*static_cast<Real>(i+1)
        + std::exp(-std::pow(u-x2,x3)/x1);
      val += f*f;
    }
    return val;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    Ptr<std::vector<Real> > eg
      = dynamic_cast<DualScaledStdVector<Real>&>(g).getVector();
    Ptr<const std::vector<Real> > ex
      = dynamic_cast<const PrimalScaledStdVector<Real>&>(x).getVector();
    g.zero();

    Real f(0), df1(0), df2(0), df3(0);
    Real u(0), tmp(0), tmp0(0), tmp1(0);
    Real x1 = (*ex)[0], x2 = (*ex)[1], x3 = (*ex)[2];
    Real x1sqr = x1*x1;
    for ( uint i = 0; i < u_size_; i++ ) {
      u    = u_vec_[i];
      tmp0 = std::pow(u-x2,x3);
      tmp1 = std::pow(u-x2,x3-static_cast<Real>(1));
      tmp  = std::exp(-tmp0/x1);

      f    = -static_cast<Real>(0.01)*static_cast<Real>(i+1) + tmp;

      df1  = tmp*tmp0/x1sqr;
      df2  = tmp*x3*tmp1/x1;
      df3  = tmp*tmp0*std::log(u-x2)/x1;

      (*eg)[0] += static_cast<Real>(2)*f*df1;
      (*eg)[1] += static_cast<Real>(2)*f*df2;
      (*eg)[2] += static_cast<Real>(2)*f*df3;
    }
  }
#if USE_HESSVEC
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Ptr<std::vector<Real> > ehv
      = dynamic_cast<DualScaledStdVector<Real>&>(hv).getVector();
    Ptr<const std::vector<Real> > ev
      = dynamic_cast<const PrimalScaledStdVector<Real>&>(v).getVector();
    Ptr<const std::vector<Real> > ex
      = dynamic_cast<const PrimalScaledStdVector<Real>&>(x).getVector();
    hv.zero();

    Real f(0);
    Real df1(0),  df2(0),  df3(0);
    Real df11(0), df12(0), df13(0);
    Real df21(0), df22(0), df23(0);
    Real df31(0), df32(0), df33(0);
    Real u(0), tmp(0), tmp0(0), tmp1(0), tmp2(0), tmp3(0), tmp4(0);
    Real x1 = (*ex)[0], x2 = (*ex)[1], x3 = (*ex)[2];
    Real v1 = (*ev)[0], v2 = (*ev)[1], v3 = (*ev)[2];
    Real x1sqr = x1*x1, x1cub = x1sqr*x1, x1quar = x1cub*x1;
    for ( uint i = 0; i < u_size_; i++ ) {
      u = u_vec_[i];
      tmp0 = std::pow(u-x2,x3);
      tmp1 = std::pow(u-x2,x3-static_cast<Real>(1));
      tmp2 = std::pow(u-x2,static_cast<Real>(2)*(x3-static_cast<Real>(1)));
      tmp3 = std::pow(u-x2,x3-static_cast<Real>(2));
      tmp4 = std::pow(u-x2,static_cast<Real>(2)*x3-static_cast<Real>(1));
      tmp  = std::exp(-tmp0/x1);

      f = -static_cast<Real>(0.01)*static_cast<Real>(i+1) + tmp;

      df1 = tmp*tmp0/x1sqr;
      df2 = tmp*x3*tmp1/x1;
      df3 = tmp*tmp0*std::log(u-x2)/x1;

      df11 = tmp0*tmp*(tmp0-static_cast<Real>(2)*x1)/x1quar;
      df12 = x3*tmp1*tmp*(tmp0-x1)/x1cub;
      df13 = tmp0*std::log(u-x2)*tmp*(x1-tmp0)/x1cub;

      df21 = df12;
      df22 = x3*x3*tmp2*tmp/(x1*x1)
             -(x3-static_cast<Real>(1))*x3*tmp3*tmp/x1;
      df23 = -x3*tmp4*std::log(u-x2)*tmp/x1sqr
             +tmp1*tmp/x1 + x3*tmp1*std::log(u-x2)*tmp/x1;

      df31 = df13;
      df32 = df23;
      df33 = tmp0*std::pow(std::log(u-x2),2)*tmp*(tmp0-x1)/x1sqr;

      (*ehv)[0] += static_cast<Real>(2)*(f*(df11*v1 + df12*v2 + df13*v3)
                   + df1*(df1*v1 + df2*v2 + df3*v3));
      (*ehv)[1] += static_cast<Real>(2)*(f*(df21*v1 + df22*v2 + df23*v3)
                   + df2*(df1*v1 + df2*v2 + df3*v3));
      (*ehv)[2] += static_cast<Real>(2)*(f*(df31*v1 + df32*v2 + df33*v3)
                   + df3*(df1*v1 + df2*v2 + df3*v3));
    }
  }
#endif
};

template<class Real>
class getHS25 : public TestProblem<Real> {
private:
  int n_;
  Ptr<std::vector<Real> > scale_;

public:
  getHS25(void) {
    // Problem dimension 
    n_ = 3;
    // Set up vector scaling
    scale_ = makePtr<std::vector<Real>>(n_,0);
    (*scale_)[0] = static_cast<Real>(1.e-4);
    (*scale_)[1] = static_cast<Real>(1.e-3);
    (*scale_)[2] = static_cast<Real>(0.5);
  }

  Ptr<Objective<Real>> getObjective(void) const {
    // Instantiate Objective Function
    return makePtr<Objective_HS25<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Get Initial Guess
    Ptr<std::vector<Real> > x0p = makePtr<std::vector<Real>>(n_,0);
    (*x0p)[0] = static_cast<Real>(100);
    (*x0p)[1] = static_cast<Real>(12.5);
    (*x0p)[2] = static_cast<Real>(3);
    return makePtr<PrimalScaledStdVector<Real>>(x0p,scale_);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Get Solution
    Ptr<std::vector<Real> > xp = makePtr<std::vector<Real>>(n_,0);
    (*xp)[0] = static_cast<Real>(50);
    (*xp)[1] = static_cast<Real>(25);
    (*xp)[2] = static_cast<Real>(1.5);
    return makePtr<PrimalScaledStdVector<Real>>(xp,scale_);
  }

  Ptr<BoundConstraint<Real>> getBoundConstraint(void) const {
    // Instantiate BoundConstraint
    Ptr<std::vector<Real> > lp = makePtr<std::vector<Real>>(n_,0);
    (*lp)[0] = static_cast<Real>(0.1);
    (*lp)[1] = static_cast<Real>(0);
    (*lp)[2] = static_cast<Real>(0);
    Ptr<Vector<Real> > l = makePtr<StdVector<Real>>(lp);
    Ptr<std::vector<Real> > up = makePtr<std::vector<Real>>(n_,0);
    (*up)[0] = static_cast<Real>(100);
    (*up)[1] = static_cast<Real>(25.6);
    (*up)[2] = static_cast<Real>(5);
    Ptr<Vector<Real> > u = makePtr<StdVector<Real>>(up);
    return makePtr<Bounds<Real>>(l,u);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
