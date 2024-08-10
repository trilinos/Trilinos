// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BRENTS_H
#define ROL_BRENTS_H

/** \class ROL::Brents
    \brief Implements a Brent's method line search.
*/

#include "ROL_LineSearch.hpp"
#include "ROL_BackTracking.hpp"

namespace ROL { 

template<class Real>
class Brents : public LineSearch<Real> {
private:
  Real tol_;
  int niter_;
  bool test_;

  ROL::Ptr<Vector<Real> > xnew_; 
//  ROL::Ptr<LineSearch<Real> > btls_;

public:

  virtual ~Brents() {}

  // Constructor
  Brents( ROL::ParameterList &parlist ) : LineSearch<Real>(parlist) {
    Real oem10(1.e-10);
    ROL::ParameterList &list
      = parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").sublist("Brent's");
    tol_ = list.get("Tolerance",oem10);
    niter_ = list.get("Iteration Limit",1000);
    test_ = list.get("Run Test Upon Initialization",true);
//    tol_ = parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").get("Bracketing Tolerance",1.e-8);
//    btls_ = ROL::makePtr<BackTracking<Real>>(parlist);
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g,
                   Objective<Real> &obj, BoundConstraint<Real> &con ) {
    LineSearch<Real>::initialize(x,s,g,obj,con);
    xnew_ = x.clone();
//    btls_->initialize(x,s,g,obj,con);

    if ( test_ ) {
      if ( test_brents() ) {
        std::cout << "Brent's Test Passed!\n";
      }
      else {
        std::cout << "Brent's Test Failed!\n";
      }
    }
  }

  // Find the minimum of phi(alpha) = f(x + alpha*s) using Brent's method
  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj, BoundConstraint<Real> &con ) {
    ls_neval = 0; ls_ngrad = 0;

    // Get initial line search parameter
    alpha = LineSearch<Real>::getInitialAlpha(ls_neval,ls_ngrad,fval,gs,x,s,obj,con);

    // TODO: Bracketing

    // Run Brents
    ROL::Ptr<typename LineSearch<Real>::ScalarFunction> phi
      = ROL::makePtr<typename LineSearch<Real>::Phi>(*xnew_,x,s,obj,con);
    int neval = 0;
    Real A(0), B = alpha;
    run_brents(neval, fval, alpha, *phi, A, B);
    ls_neval += neval;
  }

private:
  void run_brents(int &neval, Real &fval, Real &alpha,
                  typename LineSearch<Real>::ScalarFunction &phi,
                  const Real A, const Real B) const {
    neval = 0;
    // ---> Set algorithmic constants
    const Real zero(0), half(0.5), one(1), two(2), three(3), five(5);
    const Real c   = half*(three - std::sqrt(five));
    const Real eps = std::sqrt(ROL_EPSILON<Real>());
    // ---> Set end points and initial guess
    Real a = A, b = B;
    alpha = a + c*(b-a);
    // ---> Evaluate function
    Real fx = phi.value(alpha);
    neval++;
    // ---> Initialize algorithm storage
    Real v = alpha, w = v, u(0), fu(0);
    Real p(0), q(0), r(0), d(0), e(0);
    Real fv = fx, fw = fx, tol(0), t2(0), m(0);
    for (int i = 0; i < niter_; i++) {
      m = half*(a+b);
      tol = eps*std::abs(alpha) + tol_; t2 = two*tol;
      // Check stopping criterion
      if (std::abs(alpha-m) <= t2 - half*(b-a)) {
        break;
      }
      p = zero; q = zero; r = zero;
      if ( std::abs(e) > tol ) {
        // Fit parabola
        r = (alpha-w)*(fx-fv);         q = (alpha-v)*(fx-fw);
        p = (alpha-v)*q - (alpha-w)*r; q = two*(q-r);
        if ( q > zero ) {
          p *= -one;
        }
        q = std::abs(q);
        r = e; e = d;
      }
      if ( std::abs(p) < std::abs(half*q*r) && p > q*(a-alpha) && p < q*(b-alpha) ) {
        // A parabolic interpolation step
        d = p/q; u = alpha + d;
        // f must not be evaluated too close to a or b
        if ( (u - a) < t2 || (b - u) < t2 ) {
          d = (alpha < m) ? tol : -tol;
        }
      }
      else {
        // A golden section step
        e = ((alpha < m) ? b : a) - alpha; d = c*e;
      }
      // f must not be evaluated too close to alpha
      u  = alpha + ((std::abs(d) >= tol) ? d : ((d > zero) ? tol : -tol));
      fu = phi.value(u);
      neval++;
      // Update a, b, v, w, and alpha
      if ( fu <= fx ) {
        if ( u < alpha ) {
          b = alpha;
        }
        else {
          a = alpha;
        }
        v = w; fv = fw; w = alpha; fw = fx; alpha = u; fx = fu;
      }
      else {
        if ( u < alpha ) {
          a = u;
        }
        else {
          b = u;
        }
        if ( fu <= fw || w == alpha ) {
          v = w; fv = fw; w = u; fw = fu;
        }
        else if ( fu <= fv || v == alpha || v == w ) {
          v = u; fv = fu;
        }
      }
    }
    fval = fx;
  }

  class testFunction : public LineSearch<Real>::ScalarFunction {
  public:
    Real value(const Real x) {
      Real val(0), I(0), two(2), five(5);
      for (int i = 0; i < 20; i++) {
        I = (Real)(i+1);
        val += std::pow((two*I - five)/(x-(I*I)),two);
      }
      return val;
    }
  };

  bool test_brents(void) const {
    ROL::Ptr<typename LineSearch<Real>::ScalarFunction> phi
       = ROL::makePtr<testFunction>();
    Real A(0), B(0), alpha(0), fval(0);
    Real error(0), error_i(0);
    Real zero(0), two(2), three(3);
    int neval = 0;
    std::vector<Real> fvector(19,zero), avector(19,zero);
    fvector[0]  = 3.6766990169; avector[0]  =   3.0229153;
    fvector[1]  = 1.1118500100; avector[1]  =   6.6837536;
    fvector[2]  = 1.2182217637; avector[2]  =  11.2387017;
    fvector[3]  = 2.1621103109; avector[3]  =  19.6760001;
    fvector[4]  = 3.0322905193; avector[4]  =  29.8282273;
    fvector[5]  = 3.7583856477; avector[5]  =  41.9061162;
    fvector[6]  = 4.3554103836; avector[6]  =  55.9535958;
    fvector[7]  = 4.8482959563; avector[7]  =  71.9856656;
    fvector[8]  = 5.2587585400; avector[8]  =  90.0088685;
    fvector[9]  = 5.6036524295; avector[9]  = 110.0265327;
    fvector[10] = 5.8956037976; avector[10] = 132.0405517;
    fvector[11] = 6.1438861542; avector[11] = 156.0521144;
    fvector[12] = 6.3550764593; avector[12] = 182.0620604;
    fvector[13] = 6.5333662003; avector[13] = 210.0711010;
    fvector[14] = 6.6803639849; avector[14] = 240.0800483;
    fvector[15] = 6.7938538365; avector[15] = 272.0902669;
    fvector[16] = 6.8634981053; avector[16] = 306.1051233;
    fvector[17] = 6.8539024631; avector[17] = 342.1369454;
    fvector[18] = 6.6008470481; avector[18] = 380.2687097;
    for ( int i = 0; i < 19; i++ ) {
      A = std::pow((Real)(i+1),two);
      B = std::pow((Real)(i+2),two);
      run_brents(neval, fval, alpha, *phi, A, B);
      error_i = std::max(std::abs(fvector[i]-fval)/fvector[i],
                         std::abs(avector[i]-alpha)/avector[i]);
      error = std::max(error,error_i);
    }
    return (error < three*(std::sqrt(ROL_EPSILON<Real>())*avector[18]+tol_)) ? true : false;
  }

};

}

#endif
