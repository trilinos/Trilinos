// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for Poisson optimal control.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_POISSONCONTROL_HPP
#define ROL_POISSONCONTROL_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"

namespace ROL {
namespace ZOO {

/** \brief Poisson distributed control.
 */
template<class Real>
class Objective_PoissonControl : public Objective<Real> {

typedef std::vector<Real>  vector;
typedef Vector<Real>       V;
typedef StdVector<Real>    SV;  

typedef typename vector::size_type uint;

private:
  Real alpha_;

  ROL::Ptr<const vector> getVector( const V& x ) {
    
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::Ptr<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

public:

  Objective_PoissonControl(Real alpha = 1.e-4) : alpha_(alpha) {}

  void apply_mass(Vector<Real> &Mz, const Vector<Real> &z ) {

    
    ROL::Ptr<const vector> zp = getVector(z);
    ROL::Ptr<vector> Mzp = getVector(Mz);

    uint  n = zp->size();
    Real h = 1.0/((Real)n+1.0);
    for (uint i=0; i<n; i++) {
      if ( i == 0 ) {
        (*Mzp)[i] = h/6.0*(4.0*(*zp)[i] + (*zp)[i+1]);
      }
      else if ( i == n-1 ) {
        (*Mzp)[i] = h/6.0*((*zp)[i-1] + 4.0*(*zp)[i]);
      }
      else {
        (*Mzp)[i] = h/6.0*((*zp)[i-1] + 4.0*(*zp)[i] + (*zp)[i+1]);
      }
    }
  }

  void solve_poisson(Vector<Real> & u, const Vector<Real> & z) {

    
    

    ROL::Ptr<vector> up = getVector(u);

    uint  n = up->size();
    Real h = 1.0/((Real)n+1.0);
    SV b( ROL::makePtr<vector>(n,0.0) );
    ROL::Ptr<vector> bp = getVector(b);
    apply_mass(b,z);

    Real d   =  2.0/h;
    Real o   = -1.0/h;
    Real m   = 0.0;
    vector c(n,o);
    c[0]     = c[0]/d;
    (*up)[0] = (*bp)[0]/d;
    for ( uint i = 1; i < n; i++ ) {
      m        = 1.0/(d - o*c[i-1]);
      c[i]     = c[i]*m;
      (*up)[i] = ( (*bp)[i] - o*(*up)[i-1] )*m;
    }
    for ( uint i = n-1; i > 0; i-- ) {
      (*up)[i-1] = (*up)[i-1] - c[i-1]*(*up)[i];
    }
  }

  Real evaluate_target(Real x) {
    Real val = 1.0/3.0*std::pow(x,4.0) - 2.0/3.0*std::pow(x,3.0) + 1.0/3.0*x + 8.0*alpha_;
    return val;
  }

  Real value( const Vector<Real> &z, Real &tol ) {

    
    
    ROL::Ptr<const vector> zp = getVector(z);
    uint n    = zp->size();
    Real h    = 1.0/((Real)n+1.0);
    // SOLVE STATE EQUATION
    SV u( ROL::makePtr<vector>(n,0.0) );
    solve_poisson(u,z);
    ROL::Ptr<vector> up = getVector(u);

    Real val  = 0.0;
    Real res  = 0.0;
    Real res1 = 0.0;
    Real res2 = 0.0;
    Real res3 = 0.0;
    for (uint i=0; i<n; i++) {
      res = alpha_*(*zp)[i];
      if ( i == 0 ) {
        res *= h/6.0*(4.0*(*zp)[i] + (*zp)[i+1]);
        res1 = (*up)[i]-evaluate_target((Real)(i+1)*h);
        res2 = (*up)[i+1]-evaluate_target((Real)(i+2)*h);
        res += h/6.0*(4.0*res1 + res2)*res1;
      }
      else if ( i == n-1 ) {
        res *= h/6.0*((*zp)[i-1] + 4.0*(*zp)[i]);
        res1 = (*up)[i-1]-evaluate_target((Real)(i)*h);
        res2 = (*up)[i]-evaluate_target((Real)(i+1)*h);
        res += h/6.0*(res1 + 4.0*res2)*res2;
      }
      else {
        res *= h/6.0*((*zp)[i-1] + 4.0*(*zp)[i] + (*zp)[i+1]);
        res1 = (*up)[i-1]-evaluate_target((Real)(i)*h);
        res2 = (*up)[i]-evaluate_target((Real)(i+1)*h);
        res3 = (*up)[i+1]-evaluate_target((Real)(i+2)*h);
        res += h/6.0*(res1 + 4.0*res2 + res3)*res2;
      }
      val += 0.5*res;
    }
    return val;
 }

  void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {

    
    
    ROL::Ptr<const vector> zp = getVector(z);
    ROL::Ptr<vector> gp = getVector(g);

    uint  n = zp->size();
    Real h = 1.0/((Real)n+1.0);

    // SOLVE STATE EQUATION
    SV u( ROL::makePtr<vector>(n,0.0) );
    solve_poisson(u,z);
    ROL::Ptr<vector> up = getVector(u);

    // SOLVE ADJOINT EQUATION
    StdVector<Real> res( ROL::makePtr<std::vector<Real>>(n,0.0) );
    ROL::Ptr<vector> rp = getVector(res);

    for (uint i=0; i<n; i++) {
      (*rp)[i] = -((*up)[i]-evaluate_target((Real)(i+1)*h));
    }

    SV p( ROL::makePtr<vector>(n,0.0) );
    solve_poisson(p,res);
    ROL::Ptr<vector> pp = getVector(p);

    Real res1 = 0.0;
    Real res2 = 0.0;
    Real res3 = 0.0;
    for (uint i=0; i<n; i++) {
      if ( i == 0 ) {
        res1 = alpha_*(*zp)[i] - (*pp)[i];
        res2 = alpha_*(*zp)[i+1] - (*pp)[i+1];
        (*gp)[i] = h/6.0*(4.0*res1 + res2);
      }
      else if ( i == n-1 ) {
        res1 = alpha_*(*zp)[i-1] - (*pp)[i-1];
        res2 = alpha_*(*zp)[i] - (*pp)[i];
        (*gp)[i] = h/6.0*(res1 + 4.0*res2);
      }
      else {
        res1 = alpha_*(*zp)[i-1] - (*pp)[i-1];
        res2 = alpha_*(*zp)[i] - (*pp)[i];
        res3 = alpha_*(*zp)[i+1] - (*pp)[i+1];
        (*gp)[i] = h/6.0*(res1 + 4.0*res2 + res3);
      }
    }
  }
#if USE_HESSVEC
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {

    
    
    ROL::Ptr<const vector> zp = getVector(z);
    ROL::Ptr<const vector> vp = getVector(v);
    ROL::Ptr<vector> hvp = getVector(hv);

    uint  n = zp->size();
    Real h = 1.0/((Real)n+1.0);

    // SOLVE STATE EQUATION
    SV u( ROL::makePtr<vector>(n,0.0) );
    solve_poisson(u,v);
    ROL::Ptr<vector> up = getVector(u);

    // SOLVE ADJOINT EQUATION
    SV p( ROL::makePtr<vector>(n,0.0) );

    solve_poisson(p,u);
    ROL::Ptr<vector> pp = getVector(p);  

    Real res1 = 0.0;
    Real res2 = 0.0;
    Real res3 = 0.0;
    for (uint i=0; i<n; i++) {
      if ( i == 0 ) {
        res1 = alpha_*(*vp)[i] + (*pp)[i];
        res2 = alpha_*(*vp)[i+1] + (*pp)[i+1];
        (*hvp)[i] = h/6.0*(4.0*res1 + res2);
      }
      else if ( i == n-1 ) {
        res1 = alpha_*(*vp)[i-1] + (*pp)[i-1];
        res2 = alpha_*(*vp)[i] + (*pp)[i];
        (*hvp)[i] = h/6.0*(res1 + 4.0*res2);
      }
      else {
        res1 = alpha_*(*vp)[i-1] + (*pp)[i-1];
        res2 = alpha_*(*vp)[i] + (*pp)[i];
        res3 = alpha_*(*vp)[i+1] + (*pp)[i+1];
        (*hvp)[i] = h/6.0*(res1 + 4.0*res2 + res3);
      }
    }
  }
#endif
};

template<class Real>
class getPoissonControl : public TestProblem<Real> {
public:
  getPoissonControl(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    // Instantiate Objective Function
    return ROL::makePtr<Objective_PoissonControl<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem dimension
    int n = 512;
    // Get Initial Guess
    ROL::Ptr<std::vector<Real> > x0p = ROL::makePtr<std::vector<Real>>(n,0.0);
    for (int i=0; i<n; i++) {
      (*x0p)[i] = 0.0;
    }
    return ROL::makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem dimension
    int n = 512;
    // Get Solution
    ROL::Ptr<std::vector<Real> > xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    Real h = 1.0/((Real)n+1.0), pt = 0.0;
    for( int i = 0; i < n; i++ ) {
      pt = (Real)(i+1)*h;
      (*xp)[i] = 4.0*pt*(1.0-pt);
    }
    return ROL::makePtr<StdVector<Real>>(xp);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
