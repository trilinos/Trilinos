// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_EXPECTATIONQUAD_HPP
#define ROL_EXPECTATIONQUAD_HPP

#include "ROL_Types.hpp"

/** @ingroup risk_group
    \class ROL::ExpectationQuad
    \brief Provides a general interface for risk and error measures generated
           through the expectation risk quadrangle.

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
class ExpectationQuad {
public:
  virtual ~ExpectationQuad(void) {}
  ExpectationQuad(void) {}

  /** \brief Evaluate the scalar regret function at x.

      @param[in]   x      is the scalar input
      @param[in]   deriv  is the derivative order

      This function returns \f$v(x)\f$ or a derivative of \f$v(x)\f$.
  */
  virtual Real regret(Real x, int deriv = 0) = 0;

  /** \brief Evaluate the scalar error function at x.
      
      @param[in]   x      is the scalar input
      @param[in]   deriv  is the derivative order

      This function returns \f$e(x)\f$ or a derivative of \f$e(x)\f$.
  */
  virtual Real error(Real x, int deriv = 0) {
    const Real one(1), zero(0);
    Real X = (deriv==0 ? x : (deriv==1 ? one : zero));
    return regret(x,deriv) - X;
  }

  /** \brief Run default derivative tests for the scalar regret function.
  */
  virtual void check(void) {
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
};

}

#endif
