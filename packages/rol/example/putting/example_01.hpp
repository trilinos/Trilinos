// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.hpp
    \brief Provides definitions of equality constraint and objective for
           example_01.
*/

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_StdVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_KrylovFactory.hpp"

template<class Real>
class PuttingConstraint : public ROL::Constraint_SimOpt<Real> {
/******************************************************************************/
/* DATA LAYOUT FOR SIMULATION VARIABLE:                                       */
/*   0:n            X-position                                                */
/*   n+1:2n+1       Y-position                                                */
/*   2n+2:3n+2      Z-position                                                */
/*   3n+3:4n+3      X-velocity                                                */
/*   4n+4:5n+4      Y-velocity                                                */
/*   5n+5:6n+5      Z-velocity                                                */
/*   6n+6:7n+6      X-acceleration                                            */
/*   7n+7:8n+7      Y-acceleration                                            */
/*   8n+8:9n+8      Z-acceleration                                            */
/*   9n+9           Final time                                                */
/* DATA LAYOUT FOR CONSTRAINT VARIABLE:                                       */
/*   0:n            X-velocity definition                                     */
/*   n+1:2n+1       Y-velocity definition                                     */
/*   2n+2:3n+2      Z-velocity definition                                     */
/*   3n+3:4n+3      X-acceleration definition                                 */
/*   4n+4:5n+4      Y-acceleration definition                                 */
/*   5n+5:6n+5      Z-acceleration definition                                 */
/*   6n+6:7n+6      X-Newton's 2nd Law                                        */
/*   7n+7:8n+7      Y-Newton's 2nd Law                                        */
/*   8n+8:9n+7      Z-position definition                                     */
/*   9n+8           X-final position                                          */
/*   9n+9           Y-final position                                          */
/******************************************************************************/

private:
  Real g_;  // Acceleration due to gravity
  Real m_;  // Mass of a golf ball
  Real x0_; // x-coordinate of starting point
  Real y0_; // y-coordinate of starting point
  Real xn_; // x-coordinate of ending point
  Real yn_; // y-coordinate of ending point
  Real mu_; // Friction coefficient
  int n_;   // Number of time steps

  int getSize(const std::vector<Real> &u) const {
    return static_cast<int>(u.size()-10)/9;
  }

  void parseState(std::vector<Real> &x1, std::vector<Real> &x2, std::vector<Real> &x3,
                  std::vector<Real> &v1, std::vector<Real> &v2, std::vector<Real> &v3,
                  std::vector<Real> &a1, std::vector<Real> &a2, std::vector<Real> &a3,
                  Real &T, const std::vector<Real> &u) const {
    int n = getSize(u);
    T = u[9*n+9];
    for (int i = 0; i < n+1; ++i) {
      x1[i] = u[i];       x2[i] = u[n+1+i];   x3[i] = u[2*n+2+i];
      v1[i] = u[3*n+3+i]; v2[i] = u[4*n+4+i]; v3[i] = u[5*n+5+i];
      a1[i] = u[6*n+6+i]; a2[i] = u[7*n+7+i]; a3[i] = u[8*n+8+i];
    }
  }

  void parseResidual(std::vector<Real> &cv1, std::vector<Real> &cv2, std::vector<Real> &cv3,
                     std::vector<Real> &ca1, std::vector<Real> &ca2, std::vector<Real> &ca3,
                     std::vector<Real> &cn1, std::vector<Real> &cn2, std::vector<Real> &cx3,
                     Real &cx1, Real &cx2, const std::vector<Real> &c) const {
    int n = getSize(c);
    cx1 = c[9*n+8];
    cx2 = c[9*n+9];
    for (int i = 0; i < n+1; ++i) {
      cv1[i] = c[i];       cv2[i] = c[n+1+i];   cv3[i] = c[2*n+2+i];
      ca1[i] = c[3*n+3+i]; ca2[i] = c[4*n+4+i]; ca3[i] = c[5*n+5+i];
      cn1[i] = c[6*n+6+i]; cn2[i] = c[7*n+7+i];
    }
    for (int i = 0; i < n; ++i) {
      cx3[i] = c[8*n+8+i];
    }
  }

  ROL::Ptr<const std::vector<Real> > getConstVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<std::vector<Real> > getVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
  }
  
  Real S(const Real x, const Real y) const {
    return static_cast<Real>(-0.3)*std::atan(y) + static_cast<Real>(0.05)*(x+y);
  }

  Real dSdx(const Real x, const Real y) const {
    return static_cast<Real>(0.05);
  }

  Real dSdy(const Real x, const Real y) const {
    return static_cast<Real>(-0.3)/(static_cast<Real>(1)+y*y) + static_cast<Real>(0.05);
  }

  Real d2Sdx2(const Real x, const Real y) const {
    return static_cast<Real>(0);
  }

  Real d2Sdy2(const Real x, const Real y) const {
    return static_cast<Real>(0.6)*y/std::pow(static_cast<Real>(1)+y*y,2);
  }

  Real d2Sdxdy(const Real x, const Real y) const {
    return static_cast<Real>(0);
  }

  void N(Real &n1, Real &n2, Real &n3,
         const Real x1, const Real x2, const Real x3,
         const Real a1, const Real a2, const Real a3) const {
    Real dZdx = dSdx(x1,x2), dZdy = dSdy(x1,x2), one(1);
    n3 = m_*(g_-a1*dZdx-a2*dZdy+a3)/(dZdx*dZdx+dZdy*dZdy+one);
    n1 = -dZdx*n3;
    n2 = -dZdy*n3;
  }

  void dNdx(Real &dv1, Real &dv2, Real &dv3,
            const Real v1, const Real v2, const Real v3,
            const Real x1, const Real x2, const Real x3,
            const Real a1, const Real a2, const Real a3,
            const bool trans = false) const {
    Real one(1), two(2);
    Real dZdx = dSdx(x1,x2), dZdy = dSdy(x1,x2);
    Real d2Zdx2 = d2Sdx2(x1,x2), d2Zdxdy = d2Sdxdy(x1,x2), d2Zdy2 = d2Sdy2(x1,x2);
    Real n1(0), n2(0), n3(0);
    N(n1,n2,n3,x1,x2,x3,a1,a2,a3);

    Real dn3dx1 = m_*(-a1*d2Zdx2-a2*d2Zdxdy)/(dZdx*dZdx+dZdy*dZdy+one)
                 -m_*(g_-a1*dZdx-a2*dZdy+a3)*two*(dZdx*d2Zdx2+dZdy*d2Zdxdy)/std::pow((dZdx*dZdx+dZdy*dZdy+one),2);
    Real dn3dx2 = m_*(-a1*d2Zdxdy-a2*d2Zdy2)/(dZdx*dZdx+dZdy*dZdy+one)
                 -m_*(g_-a1*dZdx-a2*dZdy+a3)*two*(dZdx*d2Zdxdy+dZdy*d2Zdy2)/std::pow((dZdx*dZdx+dZdy*dZdy+one),2);
    Real dn3dx3 = static_cast<Real>(0);

    Real dn1dx1 = -d2Zdx2*n3 - dZdx*dn3dx1;
    Real dn1dx2 = -d2Zdxdy*n3 - dZdx*dn3dx2;
    Real dn1dx3 = static_cast<Real>(0);

    Real dn2dx1 = -d2Zdxdy*n3 - dZdy*dn3dx1;
    Real dn2dx2 = -d2Zdy2*n3 - dZdy*dn3dx2;
    Real dn2dx3 = static_cast<Real>(0);

    if (!trans) {
      dv1 = dn1dx1*v1 + dn1dx2*v2 + dn1dx3*v3;
      dv2 = dn2dx1*v1 + dn2dx2*v2 + dn2dx3*v3;
      dv3 = dn3dx1*v1 + dn3dx2*v2 + dn3dx3*v3;
    }
    else {
      dv1 = dn1dx1*v1 + dn2dx1*v2 + dn3dx1*v3;
      dv2 = dn1dx2*v1 + dn2dx2*v2 + dn3dx2*v3;
      dv3 = dn1dx3*v1 + dn2dx3*v2 + dn3dx3*v3;
    }
  }

  void dNda(Real &dv1, Real &dv2, Real &dv3,
            const Real v1, const Real v2, const Real v3,
            const Real x1, const Real x2, const Real x3,
            const Real a1, const Real a2, const Real a3,
            const bool trans = false) const {
    Real one(1);
    Real dZdx = dSdx(x1,x2), dZdy = dSdy(x1,x2);

    Real dn3da1 = m_*(-dZdx)/(dZdx*dZdx+dZdy*dZdy+one);
    Real dn3da2 = m_*(-dZdy)/(dZdx*dZdx+dZdy*dZdy+one);
    Real dn3da3 = m_*one/(dZdx*dZdx+dZdy*dZdy+one);

    Real dn1da1 = -dZdx*dn3da1;
    Real dn1da2 = -dZdx*dn3da2;
    Real dn1da3 = -dZdx*dn3da3;

    Real dn2da1 = -dZdy*dn3da1;
    Real dn2da2 = -dZdy*dn3da2;
    Real dn2da3 = -dZdy*dn3da3;

    if (!trans) {
      dv1 = dn1da1*v1 + dn1da2*v2 + dn1da3*v3;
      dv2 = dn2da1*v1 + dn2da2*v2 + dn2da3*v3;
      dv3 = dn3da1*v1 + dn3da2*v2 + dn3da3*v3;
    }
    else {
      dv1 = dn1da1*v1 + dn2da1*v2 + dn3da1*v3;
      dv2 = dn1da2*v1 + dn2da2*v2 + dn3da2*v3;
      dv3 = dn1da3*v1 + dn2da3*v2 + dn3da3*v3;
    }
  }

/*
  void N(Real &n1, Real &n2, Real &n3,
         const Real x1, const Real x2, const Real x3,
         const Real a1, const Real a2, const Real a3) const {
    Real dZdx = dSdx(x1,x2), dZdy = dSdy(x1,x2), one(1);
    n3 = m_*g_/(dZdx*dZdx+dZdy*dZdy+one);
    n1 = -dZdx*n3;
    n2 = -dZdy*n3;
  }

  void dNdx(Real &dv1, Real &dv2, Real &dv3,
            const Real v1, const Real v2, const Real v3,
            const Real x1, const Real x2, const Real x3,
            const Real a1, const Real a2, const Real a3) const {
    Real one(1), two(2);
    Real dZdx = dSdx(x1,x2), dZdy = dSdy(x1,x2);
    Real d2Zdx2 = d2Sdx2(x1,x2), d2Zdxdy = d2Sdxdy(x1,x2), d2Zdy2 = d2Sdy2(x1,x2);
    Real n1(0), n2(0), n3(0);
    N(n1,n2,n3,x1,x2,x3,a1,a2,a3);

    Real dn3dx1 = -m_*g_*two*(dZdx*d2Zdx2+dZdy*d2Zdxdy)/std::pow((dZdx*dZdx+dZdy*dZdy+one),2);
    Real dn3dx2 = -m_*g_*two*(dZdx*d2Zdxdy+dZdy*d2Zdy2)/std::pow((dZdx*dZdx+dZdy*dZdy+one),2);
    Real dn3dx3 = static_cast<Real>(0);

    Real dn1dx1 = -d2Zdx2*n3 - dZdx*dn3dx1;
    Real dn1dx2 = -d2Zdxdy*n3 - dZdx*dn3dx2;
    Real dn1dx3 = static_cast<Real>(0);

    Real dn2dx1 = -d2Zdxdy*n3 - dZdy*dn3dx1;
    Real dn2dx2 = -d2Zdy2*n3 - dZdy*dn3dx2;
    Real dn2dx3 = static_cast<Real>(0);

    dv1 = dn1dx1*v1 + dn1dx2*v2 + dn1dx3*v3;
    dv2 = dn2dx1*v1 + dn2dx2*v2 + dn2dx3*v3;
    dv3 = dn3dx1*v1 + dn3dx2*v2 + dn3dx3*v3;
  }

  void dNda(Real &dv1, Real &dv2, Real &dv3,
            const Real v1, const Real v2, const Real v3,
            const Real x1, const Real x2, const Real x3,
            const Real a1, const Real a2, const Real a3) const {
    dv1 = static_cast<Real>(0);
    dv2 = static_cast<Real>(0);
    dv3 = static_cast<Real>(0);
  }
*/

  void F(Real &f1, Real &f2, Real &f3,
         const Real x1, const Real x2, const Real x3,
         const Real v1, const Real v2, const Real v3,
         const Real a1, const Real a2, const Real a3) const {
    Real n1(0), n2(0), n3(0);
    N(n1,n2,n3,x1,x2,x3,a1,a2,a3);
    Real nmag  = std::sqrt(n1*n1+n2*n2+n3*n3);
    Real speed = std::sqrt(v1*v1+v2*v2+v3*v3);
    f1 = -mu_*nmag*v1/speed;
    f2 = -mu_*nmag*v2/speed;
    f3 = -mu_*nmag*v3/speed;
  }

  void dFdx(Real &dv1, Real &dv2, Real &dv3,
            const Real w1, const Real w2, const Real w3,
            const Real x1, const Real x2, const Real x3,
            const Real v1, const Real v2, const Real v3,
            const Real a1, const Real a2, const Real a3,
            const bool trans = false) const {
    Real n1(0), n2(0), n3(0);
    N(n1,n2,n3,x1,x2,x3,a1,a2,a3);
    Real nmag  = std::sqrt(n1*n1+n2*n2+n3*n3);
    Real speed = std::sqrt(v1*v1+v2*v2+v3*v3);
    Real d1v(0), d2v(0), d3v(0);
    if (!trans) {
      dNdx(d1v,d2v,d3v,w1,w2,w3,x1,x2,x3,a1,a2,a3,trans);
      Real ndv = n1*d1v + n2*d2v + n3*d3v;
      dv1 = -mu_*(ndv/nmag)*(v1/speed);
      dv2 = -mu_*(ndv/nmag)*(v2/speed);
      dv3 = -mu_*(ndv/nmag)*(v3/speed);
    }
    else {
      dNdx(d1v,d2v,d3v,n1,n2,n3,x1,x2,x3,a1,a2,a3,trans);
      Real vw = v1*w1 + v2*w2 + v3*w3;
      dv1 = -mu_*(d1v/nmag)*(vw/speed);
      dv2 = -mu_*(d2v/nmag)*(vw/speed);
      dv3 = -mu_*(d3v/nmag)*(vw/speed);
    }
  }

  void dFda(Real &dv1, Real &dv2, Real &dv3,
            const Real w1, const Real w2, const Real w3,
            const Real x1, const Real x2, const Real x3,
            const Real v1, const Real v2, const Real v3,
            const Real a1, const Real a2, const Real a3,
            const bool trans = false) const {
    Real n1(0), n2(0), n3(0);
    N(n1,n2,n3,x1,x2,x3,a1,a2,a3);
    Real nmag  = std::sqrt(n1*n1+n2*n2+n3*n3);
    Real speed = std::sqrt(v1*v1+v2*v2+v3*v3);
    Real d1v(0), d2v(0), d3v(0);
    if (!trans) {
      dNda(d1v,d2v,d3v,w1,w2,w3,x1,x2,x3,a1,a2,a3,trans);
      Real ndv = n1*d1v + n2*d2v + n3*d3v;
      dv1 = -mu_*(ndv/nmag)*(v1/speed);
      dv2 = -mu_*(ndv/nmag)*(v2/speed);
      dv3 = -mu_*(ndv/nmag)*(v3/speed);
    }
    else {
      dNda(d1v,d2v,d3v,n1,n2,n3,x1,x2,x3,a1,a2,a3,trans);
      Real vw = v1*w1 + v2*w2 + v3*w3;
      dv1 = -mu_*(d1v/nmag)*(vw/speed);
      dv2 = -mu_*(d2v/nmag)*(vw/speed);
      dv3 = -mu_*(d3v/nmag)*(vw/speed);
    }
  }

  void dFdv(Real &dv1, Real &dv2, Real &dv3,
            const Real w1, const Real w2, const Real w3,
            const Real x1, const Real x2, const Real x3,
            const Real v1, const Real v2, const Real v3,
            const Real a1, const Real a2, const Real a3,
            const bool trans = false) const {
    // Real two(2);
    Real n1(0), n2(0), n3(0);
    N(n1,n2,n3,x1,x2,x3,a1,a2,a3);
    Real nmag  = std::sqrt(n1*n1+n2*n2+n3*n3);
    Real speed = std::sqrt(v1*v1+v2*v2+v3*v3);
    Real vw    = v1*w1 + v2*w2 + v3*w3;
    dv1 = -mu_*nmag*(w1/speed - v1*vw/(speed*speed*speed));
    dv2 = -mu_*nmag*(w2/speed - v2*vw/(speed*speed*speed));
    dv3 = -mu_*nmag*(w3/speed - v3*vw/(speed*speed*speed));
  }

  class Jacobian : public ROL::LinearOperator<Real> {
  private:
    const ROL::Ptr<ROL::Constraint_SimOpt<Real> > con_;
    const ROL::Ptr<const ROL::Vector<Real> > u_, z_;
  public:
    Jacobian(const ROL::Ptr<ROL::Constraint_SimOpt<Real> > &con,
             const ROL::Ptr<const ROL::Vector<Real> > &u,
             const ROL::Ptr<const ROL::Vector<Real> > &z) : con_(con), u_(u), z_(z) {}
    void apply(ROL::Vector<Real> &Jv, const ROL::Vector<Real> &v, Real &tol) const {
      con_->applyJacobian_1(Jv,v,*u_,*z_,tol);
    }
  };

  class AdjointJacobian : public ROL::LinearOperator<Real> {
  private:
    const ROL::Ptr<ROL::Constraint_SimOpt<Real> > con_;
    const ROL::Ptr<const ROL::Vector<Real> > u_, z_;
  public:
    AdjointJacobian(const ROL::Ptr<ROL::Constraint_SimOpt<Real> > &con,
                    const ROL::Ptr<const ROL::Vector<Real> > &u,
                    const ROL::Ptr<const ROL::Vector<Real> > &z) : con_(con), u_(u), z_(z) {}
    void apply(ROL::Vector<Real> &Jv, const ROL::Vector<Real> &v, Real &tol) const {
      con_->applyAdjointJacobian_1(Jv,v,*u_,*z_,tol);
    }
  };

  class Precond : public ROL::LinearOperator<Real> {
  public:
    Precond() {}
    void apply(ROL::Vector<Real> &Jv, const ROL::Vector<Real> &v, Real &tol) const {
      Jv.set(v);
    }
    void applyInverse(ROL::Vector<Real> &Jv, const ROL::Vector<Real> &v, Real &tol) const {
      Jv.set(v);
    }
  };

public:

  PuttingConstraint(const Real g = 9.8, const Real m = 0.01,
                    const Real x0 = 1.0, const Real y0 = 2.0,
                    const Real xn = 1.0, const Real yn = -2.0,
                    const Real mu = 0.07 )
    : g_(g), m_(m), x0_(x0), y0_(y0), xn_(xn), yn_(yn), mu_(mu) {}

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, 
                  const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> >       cp = getVector(c);
    ROL::Ptr<const std::vector<Real> > up = getConstVector(u);
    ROL::Ptr<const std::vector<Real> > zp = getConstVector(z);
//    const Real one(1);
    // Get number of time steps
    int n = getSize(*up);
    // Parse state vector
    Real T(0);
    std::vector<Real> x1(n+1), x2(n+1), x3(n+1); // Position
    std::vector<Real> v1(n+1), v2(n+1), v3(n+1); // Velocity
    std::vector<Real> a1(n+1), a2(n+1), a3(n+1); // Accelleration
    parseState(x1,x2,x3,v1,v2,v3,a1,a2,a3,T,*up);
    // Compute time step
    Real dt = T/static_cast<Real>(n);
    // Initial forces
    Real n1(0), n2(0), n3(0), f1(0), f2(0), f3(0);
    N(n1,n2,n3,x1[0],x2[0],x3[0],a1[0],a2[0],a3[0]);                   // Initial normal force
    F(f1,f2,f3,x1[0],x2[0],x3[0],v1[0],v2[0],v3[0],a1[0],a2[0],a3[0]); // Initial friction force
    // Build residual
    (*cp)[0]     = x1[0] - x0_;                                               // Initial x1-position
    (*cp)[n+1]   = x2[0] - y0_;                                               // Initial x2-position
    (*cp)[2*n+2] = x3[0] - S(x0_,y0_);                                        // Initial x2-position
    (*cp)[3*n+3] = v1[0] - (*zp)[0];                                          // Initial x1-velocity
    (*cp)[4*n+4] = v2[0] - (*zp)[1];                                          // Initial x2-velocity
    (*cp)[5*n+5] = v3[0] - (dSdx(x0_,y0_)*(*zp)[0] + dSdy(x0_,y0_)*(*zp)[1]); // Initial x3-velocity
    (*cp)[6*n+6] = m_*a1[0] - n1 - f1;
    (*cp)[7*n+7] = m_*a2[0] - n2 - f2;
    for (int i = 1; i < n+1; ++i) {
      // Compute forces
      N(n1,n2,n3,x1[i],x2[i],x3[i],a1[i],a2[i],a3[i]);                   // Normal force
      F(f1,f2,f3,x1[i],x2[i],x3[i],v1[i],v2[i],v3[i],a1[i],a2[i],a3[i]); // Friction force
      // Build residual
      (*cp)[i]       = (x1[i] - x1[i-1]) - dt * static_cast<Real>(0.5) * (v1[i] + v1[i-1]); // Definition of v1
      (*cp)[n+1+i]   = (x2[i] - x2[i-1]) - dt * static_cast<Real>(0.5) * (v2[i] + v2[i-1]); // Definition of v2
      (*cp)[2*n+2+i] = (x3[i] - x3[i-1]) - dt * static_cast<Real>(0.5) * (v3[i] + v3[i-1]); // Definition of v3
      (*cp)[3*n+3+i] = (v1[i] - v1[i-1]) - dt * static_cast<Real>(0.5) * (a1[i] + a1[i-1]); // Definition of a1
      (*cp)[4*n+4+i] = (v2[i] - v2[i-1]) - dt * static_cast<Real>(0.5) * (a2[i] + a2[i-1]); // Definition of a2
      (*cp)[5*n+5+i] = (v3[i] - v3[i-1]) - dt * static_cast<Real>(0.5) * (a3[i] + a3[i-1]); // Definition of a3
      (*cp)[6*n+6+i] = m_*a1[i] - n1 - f1; // Newton's 2nd law
      (*cp)[7*n+7+i] = m_*a2[i] - n2 - f2; // Newton's 2nd law
      (*cp)[8*n+7+i] = x3[i] - S(x1[i],x2[i]); // Topography constraint
    }
    (*cp)[9*n+8] = x1[n] - xn_;        // Final x1-position
    (*cp)[9*n+9] = x2[n] - yn_;        // Final x2-position
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                       const ROL::Vector<Real> &z, Real &tol) {
    //ROL::Constraint_SimOpt<Real>::applyJacobian_1(jv,v,u,z,tol);
    ROL::Ptr<std::vector<Real> >      jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real> > vp = getConstVector(v);
    ROL::Ptr<const std::vector<Real> > up = getConstVector(u);
    const Real half(0.5); // one(1), two(2);
    // Get number of time steps
    int n = getSize(*up);
    // Parse state vector
    Real T(0);
    std::vector<Real> x1(n+1), x2(n+1), x3(n+1); // Position
    std::vector<Real> v1(n+1), v2(n+1), v3(n+1); // Velocity
    std::vector<Real> a1(n+1), a2(n+1), a3(n+1); // Accelleration
    parseState(x1,x2,x3,v1,v2,v3,a1,a2,a3,T,*up);
    // Parse direction vector
    Real vT(0);
    std::vector<Real> vx1(n+1), vx2(n+1), vx3(n+1); // Position
    std::vector<Real> vv1(n+1), vv2(n+1), vv3(n+1); // Velocity
    std::vector<Real> va1(n+1), va2(n+1), va3(n+1); // Accelleration
    parseState(vx1,vx2,vx3,vv1,vv2,vv3,va1,va2,va3,vT,*vp);
    // Compute time step
    Real dt = T/static_cast<Real>(n);
    Real vdt = vT/static_cast<Real>(n);
    // Initial force jacobians
    Real DNDXv1(0), DNDXv2(0), DNDXv3(0);
    Real DNDAv1(0), DNDAv2(0), DNDAv3(0);
    Real DFDXv1(0), DFDXv2(0), DFDXv3(0);
    Real DFDVv1(0), DFDVv2(0), DFDVv3(0);
    Real DFDAv1(0), DFDAv2(0), DFDAv3(0);
    dNdx(DNDXv1,DNDXv2,DNDXv3,vx1[0],vx2[0],vx3[0],x1[0],x2[0],x3[0],a1[0],a2[0],a3[0]);
    dNda(DNDAv1,DNDAv2,DNDAv3,va1[0],va2[0],va3[0],x1[0],x2[0],x3[0],a1[0],a2[0],a3[0]);
    dFdx(DFDXv1,DFDXv2,DFDXv3,vx1[0],vx2[0],vx3[0],x1[0],x2[0],x3[0],v1[0],v2[0],v3[0],a1[0],a2[0],a3[0]);
    dFdv(DFDVv1,DFDVv2,DFDVv3,vv1[0],vv2[0],vv3[0],x1[0],x2[0],x3[0],v1[0],v2[0],v3[0],a1[0],a2[0],a3[0]);
    dFda(DFDAv1,DFDAv2,DFDAv3,va1[0],va2[0],va3[0],x1[0],x2[0],x3[0],v1[0],v2[0],v3[0],a1[0],a2[0],a3[0]);
    // Build jacobian
    (*jvp)[0]       = vx1[0]; // Initial x1-position
    (*jvp)[n+1]     = vx2[0]; // Initial x2-position
    (*jvp)[2*n+2]   = vx3[0]; // Initial x2-position
    (*jvp)[3*n+3]   = vv1[0]; // Initial x1-velocity
    (*jvp)[4*n+4]   = vv2[0]; // Initial x2-velocity
    (*jvp)[5*n+5]   = vv3[0]; // Initial x3-velocity
    (*jvp)[6*n+6]   = m_*va1[0] - (DNDXv1 + DNDAv1) - (DFDXv1 + DFDVv1 + DFDAv1);
    (*jvp)[7*n+7]   = m_*va2[0] - (DNDXv2 + DNDAv2) - (DFDXv2 + DFDVv2 + DFDAv2);
    for (int i = 1; i < n+1; ++i) {
      // Compute force jacobians
      dNdx(DNDXv1,DNDXv2,DNDXv3,vx1[i],vx2[i],vx3[i],x1[i],x2[i],x3[i],a1[i],a2[i],a3[i]);
      dNda(DNDAv1,DNDAv2,DNDAv3,va1[i],va2[i],va3[i],x1[i],x2[i],x3[i],a1[i],a2[i],a3[i]);
      dFdx(DFDXv1,DFDXv2,DFDXv3,vx1[i],vx2[i],vx3[i],x1[i],x2[i],x3[i],v1[i],v2[i],v3[i],a1[i],a2[i],a3[i]);
      dFdv(DFDVv1,DFDVv2,DFDVv3,vv1[i],vv2[i],vv3[i],x1[i],x2[i],x3[i],v1[i],v2[i],v3[i],a1[i],a2[i],a3[i]);
      dFda(DFDAv1,DFDAv2,DFDAv3,va1[i],va2[i],va3[i],x1[i],x2[i],x3[i],v1[i],v2[i],v3[i],a1[i],a2[i],a3[i]);
      // Build jacobian
      (*jvp)[i]       = (vx1[i] - vx1[i-1]) - dt * half * (vv1[i] + vv1[i-1]) - vdt * half * (v1[i] + v1[i-1]); // Derivative of x1
      (*jvp)[n+1+i]   = (vx2[i] - vx2[i-1]) - dt * half * (vv2[i] + vv2[i-1]) - vdt * half * (v2[i] + v2[i-1]); // Derivative of x2
      (*jvp)[2*n+2+i] = (vx3[i] - vx3[i-1]) - dt * half * (vv3[i] + vv3[i-1]) - vdt * half * (v3[i] + v3[i-1]); // Derivative of x3
      (*jvp)[3*n+3+i] = (vv1[i] - vv1[i-1]) - dt * half * (va1[i] + va1[i-1]) - vdt * half * (a1[i] + a1[i-1]); // Derivative of v1
      (*jvp)[4*n+4+i] = (vv2[i] - vv2[i-1]) - dt * half * (va2[i] + va2[i-1]) - vdt * half * (a2[i] + a2[i-1]); // Derivative of v2
      (*jvp)[5*n+5+i] = (vv3[i] - vv3[i-1]) - dt * half * (va3[i] + va3[i-1]) - vdt * half * (a3[i] + a3[i-1]); // Derivative of v3
      (*jvp)[6*n+6+i] = m_*va1[i] - (DNDXv1 + DNDAv1) - (DFDXv1 + DFDVv1 + DFDAv1); // Newton's 2nd law
      (*jvp)[7*n+7+i] = m_*va2[i] - (DNDXv2 + DNDAv2) - (DFDXv2 + DFDVv2 + DFDAv2); // Newton's 2nd law
      (*jvp)[8*n+7+i] = vx3[i] - dSdx(x1[i],x2[i])*vx1[i] - dSdy(x1[i],x2[i])*vx2[i];
    }
    (*jvp)[9*n+8] = vx1[n]; // Final x1-position
    (*jvp)[9*n+9] = vx2[n]; // Final x2-position
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    jv.zero();
    ROL::Ptr<std::vector<Real> >      jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real> > vp = getConstVector(v);
    ROL::Ptr<const std::vector<Real> > up = getConstVector(u);
    // Get number of time steps
    int n = getSize(*up);
    // Build jacobian
    (*jvp)[3*n+3] = -(*vp)[0];
    (*jvp)[4*n+4] = -(*vp)[1];
    (*jvp)[5*n+5] = -(dSdx(x0_,y0_)*(*vp)[0] + dSdy(x0_,y0_)*(*vp)[1]);
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                              const ROL::Vector<Real> &z, Real &tol) {
    //ROL::Constraint_SimOpt<Real>::applyAdjointJacobian_1(ajv,v,u,z,tol);
    ROL::Ptr<std::vector<Real> >      jvp = getVector(ajv);
    ROL::Ptr<const std::vector<Real> > vp = getConstVector(v);
    ROL::Ptr<const std::vector<Real> > up = getConstVector(u);
    const Real half(0.5), one(1); // two(2);
    // Get number of time steps
    int n = getSize(*up);
    // Parse state vector
    Real T(0);
    std::vector<Real> x1(n+1), x2(n+1), x3(n+1); // Position
    std::vector<Real> v1(n+1), v2(n+1), v3(n+1); // Velocity
    std::vector<Real> a1(n+1), a2(n+1), a3(n+1); // Accelleration
    parseState(x1,x2,x3,v1,v2,v3,a1,a2,a3,T,*up);
    // Compute time step
    Real dt = T/static_cast<Real>(n);
    Real rN  = one/static_cast<Real>(n);
    // Parse direction vector
    std::vector<Real> cv1(n+1), cv2(n+1), cv3(n+1);
    std::vector<Real> ca1(n+1), ca2(n+1), ca3(n+1);
    std::vector<Real> cn1(n+1), cn2(n+1), cx3(n+1);
    Real cx1(0), cx2(0);
    parseResidual(cv1,cv2,cv3,ca1,ca2,ca3,cn1,cn2,cx3,cx1,cx2,*vp);
    // Initial force jacobians
    Real DNDXv1(0), DNDXv2(0), DNDXv3(0);
    Real DNDAv1(0), DNDAv2(0), DNDAv3(0);
    Real DFDXv1(0), DFDXv2(0), DFDXv3(0);
    Real DFDVv1(0), DFDVv2(0), DFDVv3(0);
    Real DFDAv1(0), DFDAv2(0), DFDAv3(0);
    dNdx(DNDXv1,DNDXv2,DNDXv3,cn1[0],cn2[0],0,x1[0],x2[0],x3[0],a1[0],a2[0],a3[0],true);
    dNda(DNDAv1,DNDAv2,DNDAv3,cn1[0],cn2[0],0,x1[0],x2[0],x3[0],a1[0],a2[0],a3[0],true);
    dFdx(DFDXv1,DFDXv2,DFDXv3,cn1[0],cn2[0],0,x1[0],x2[0],x3[0],v1[0],v2[0],v3[0],a1[0],a2[0],a3[0],true);
    dFdv(DFDVv1,DFDVv2,DFDVv3,cn1[0],cn2[0],0,x1[0],x2[0],x3[0],v1[0],v2[0],v3[0],a1[0],a2[0],a3[0],true);
    dFda(DFDAv1,DFDAv2,DFDAv3,cn1[0],cn2[0],0,x1[0],x2[0],x3[0],v1[0],v2[0],v3[0],a1[0],a2[0],a3[0],true);
    // Build jacobian
    (*jvp)[0]     = cv1[0] - cv1[1] - DNDXv1 - DFDXv1;// - dSdx(x1[0],x2[0])*cx3[0];
    (*jvp)[n+1]   = cv2[0] - cv2[1] - DNDXv2 - DFDXv2;// - dSdy(x1[0],x2[0])*cx3[0];
    (*jvp)[2*n+2] = cv3[0] - cv3[1] - DNDXv3 - DFDXv3;// + cx3[0];
    (*jvp)[3*n+3] = ca1[0] - ca1[1] - dt*half*cv1[1] - DFDVv1;
    (*jvp)[4*n+4] = ca2[0] - ca2[1] - dt*half*cv2[1] - DFDVv2;
    (*jvp)[5*n+5] = ca3[0] - ca3[1] - dt*half*cv3[1] - DFDVv3;
    (*jvp)[6*n+6] = m_*cn1[0] - DNDAv1 - DFDAv1 - dt*half*ca1[1];
    (*jvp)[7*n+7] = m_*cn2[0] - DNDAv2 - DFDAv2 - dt*half*ca2[1];
    (*jvp)[8*n+8] =           - DNDAv3 - DFDAv3 - dt*half*ca3[1];
    (*jvp)[9*n+9] = static_cast<Real>(0);
    for (int i = 1; i < n; ++i) {
      // Compute force jacobians
      dNdx(DNDXv1,DNDXv2,DNDXv3,cn1[i],cn2[i],0,x1[i],x2[i],x3[i],a1[i],a2[i],a3[i],true);
      dNda(DNDAv1,DNDAv2,DNDAv3,cn1[i],cn2[i],0,x1[i],x2[i],x3[i],a1[i],a2[i],a3[i],true);
      dFdx(DFDXv1,DFDXv2,DFDXv3,cn1[i],cn2[i],0,x1[i],x2[i],x3[i],v1[i],v2[i],v3[i],a1[i],a2[i],a3[i],true);
      dFdv(DFDVv1,DFDVv2,DFDVv3,cn1[i],cn2[i],0,x1[i],x2[i],x3[i],v1[i],v2[i],v3[i],a1[i],a2[i],a3[i],true);
      dFda(DFDAv1,DFDAv2,DFDAv3,cn1[i],cn2[i],0,x1[i],x2[i],x3[i],v1[i],v2[i],v3[i],a1[i],a2[i],a3[i],true);
      // Build jacobian
      (*jvp)[i]       = cv1[i] - cv1[i+1] - DNDXv1 - DFDXv1 - dSdx(x1[i],x2[i])*cx3[i-1];
      (*jvp)[n+1+i]   = cv2[i] - cv2[i+1] - DNDXv2 - DFDXv2 - dSdy(x1[i],x2[i])*cx3[i-1];
      (*jvp)[2*n+2+i] = cv3[i] - cv3[i+1] - DNDXv3 - DFDXv3 + cx3[i-1];
      (*jvp)[3*n+3+i] = ca1[i] - ca1[i+1] - dt*half*(cv1[i] + cv1[i+1]) - DFDVv1;
      (*jvp)[4*n+4+i] = ca2[i] - ca2[i+1] - dt*half*(cv2[i] + cv2[i+1]) - DFDVv2;
      (*jvp)[5*n+5+i] = ca3[i] - ca3[i+1] - dt*half*(cv3[i] + cv3[i+1]) - DFDVv3;
      (*jvp)[6*n+6+i] = m_*cn1[i] - DNDAv1 - DFDAv1 - dt*half*(ca1[i]+ca1[i+1]);
      (*jvp)[7*n+7+i] = m_*cn2[i] - DNDAv2 - DFDAv2 - dt*half*(ca2[i]+ca2[i+1]);
      (*jvp)[8*n+8+i] =           - DNDAv3 - DFDAv3 - dt*half*(ca3[i]+ca3[i+1]);
      (*jvp)[9*n+9]  -= rN*half * (cv1[i]*(v1[i]+v1[i-1]) + cv2[i]*(v2[i]+v2[i-1]) + cv3[i]*(v3[i]+v3[i-1])
                                +  ca1[i]*(a1[i]+a1[i-1]) + ca2[i]*(a2[i]+a2[i-1]) + ca3[i]*(a3[i]+a3[i-1]));
    }
    // Compute force jacobians
    dNdx(DNDXv1,DNDXv2,DNDXv3,cn1[n],cn2[n],0,x1[n],x2[n],x3[n],a1[n],a2[n],a3[n],true);
    dNda(DNDAv1,DNDAv2,DNDAv3,cn1[n],cn2[n],0,x1[n],x2[n],x3[n],a1[n],a2[n],a3[n],true);
    dFdx(DFDXv1,DFDXv2,DFDXv3,cn1[n],cn2[n],0,x1[n],x2[n],x3[n],v1[n],v2[n],v3[n],a1[n],a2[n],a3[n],true);
    dFdv(DFDVv1,DFDVv2,DFDVv3,cn1[n],cn2[n],0,x1[n],x2[n],x3[n],v1[n],v2[n],v3[n],a1[n],a2[n],a3[n],true);
    dFda(DFDAv1,DFDAv2,DFDAv3,cn1[n],cn2[n],0,x1[n],x2[n],x3[n],v1[n],v2[n],v3[n],a1[n],a2[n],a3[n],true);
    // Build jacobian
    (*jvp)[n]      = cv1[n] + cx1 - DNDXv1 - DFDXv1 - dSdx(x1[n],x2[n])*cx3[n-1];
    (*jvp)[2*n+1]  = cv2[n] + cx2 - DNDXv2 - DFDXv2 - dSdy(x1[n],x2[n])*cx3[n-1];
    (*jvp)[3*n+2]  = cv3[n] - DNDXv3 - DFDXv3 + cx3[n-1];
    (*jvp)[4*n+3]  = ca1[n] - dt*half*cv1[n] - DFDVv1;
    (*jvp)[5*n+4]  = ca2[n] - dt*half*cv2[n] - DFDVv2;
    (*jvp)[6*n+5]  = ca3[n] - dt*half*cv3[n] - DFDVv3;
    (*jvp)[7*n+6]  = m_*cn1[n] - DNDAv1 - DFDAv1 - dt*half*ca1[n];
    (*jvp)[8*n+7]  = m_*cn2[n] - DNDAv2 - DFDAv2 - dt*half*ca2[n];
    (*jvp)[9*n+8]  =           - DNDAv3 - DFDAv3 - dt*half*ca3[n];
    (*jvp)[9*n+9] -= rN*half * (cv1[n]*(v1[n]+v1[n-1]) + cv2[n]*(v2[n]+v2[n-1]) + cv3[n]*(v3[n]+v3[n-1])
                             +  ca1[n]*(a1[n]+a1[n-1]) + ca2[n]*(a2[n]+a2[n-1]) + ca3[n]*(a3[n]+a3[n-1]));
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> >      jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real> > vp = getConstVector(v);
    ROL::Ptr<const std::vector<Real> > up = getConstVector(u);
    // Get number of time steps
    int n = getSize(*up);
    // Build adjoint jacobian
    (*jvp)[0] = -(*vp)[3*n+3]-dSdx(x0_,y0_)*(*vp)[5*n+5];
    (*jvp)[1] = -(*vp)[4*n+4]-dSdy(x0_,y0_)*(*vp)[5*n+5];
  }

  void applyInverseJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    int iter(0), flag(0);
    ROL::Ptr<const ROL::Vector<Real> > u_ptr = ROL::makePtrFromRef(u);
    ROL::Ptr<const ROL::Vector<Real> > z_ptr = ROL::makePtrFromRef(z);
    ROL::Ptr<ROL::Constraint_SimOpt<Real> > con = ROL::makePtrFromRef(*this);
    ROL::Ptr<ROL::LinearOperator<Real> > jac
      = ROL::makePtr<Jacobian>(con,u_ptr,z_ptr);
    ROL::Ptr<ROL::LinearOperator<Real> > precond
      = ROL::makePtr<Precond>();

    ROL::ParameterList parlist;
    parlist.sublist("General").sublist("Krylov").set("Type","GMRES");
    parlist.sublist("General").sublist("Krylov").set("Absolute Tolerance", 1e-8);
    parlist.sublist("General").sublist("Krylov").set("Relative Tolerance", 1e-4);
    parlist.sublist("General").sublist("Krylov").set("Iteration Limit", 600);
    ROL::Ptr<ROL::Krylov<Real> > krylov = ROL::KrylovFactory<Real>(parlist);

    krylov->run(ijv, *jac, v, *precond, iter, flag );
  }

  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv, const ROL::Vector<Real> &v,
                                     const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    int iter(0), flag(0);
    ROL::Ptr<const ROL::Vector<Real> > u_ptr = ROL::makePtrFromRef(u);
    ROL::Ptr<const ROL::Vector<Real> > z_ptr = ROL::makePtrFromRef(z);
    ROL::Ptr<ROL::Constraint_SimOpt<Real> > con = ROL::makePtrFromRef(*this);
    ROL::Ptr<ROL::LinearOperator<Real> > jac
      = ROL::makePtr<AdjointJacobian>(con,u_ptr,z_ptr);
    ROL::Ptr<ROL::LinearOperator<Real> > precond
      = ROL::makePtr<Precond>();

    ROL::ParameterList parlist;
    parlist.sublist("General").sublist("Krylov").set("Type","GMRES");
    parlist.sublist("General").sublist("Krylov").set("Absolute Tolerance", 1e-8);
    parlist.sublist("General").sublist("Krylov").set("Relative Tolerance", 1e-4);
    parlist.sublist("General").sublist("Krylov").set("Iteration Limit", 600);
    ROL::Ptr<ROL::Krylov<Real> > krylov = ROL::KrylovFactory<Real>(parlist);

    krylov->run(iajv, *jac, v, *precond, iter, flag );
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Constraint_SimOpt<Real>::applyAdjointHessian_11(ahwv,w,v,u,z,tol);
  }
  
  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }
  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }
  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }
};

template<class Real>
class PuttingObjective : public ROL::Objective_SimOpt<Real> {
/******************************************************************************/
/* DATA LAYOUT FOR SIMULATION VARIABLE:                                       */
/*   0:n            X-position                                                */
/*   n+1:2n+1       Y-position                                                */
/*   2n+2:3n+2      Z-position                                                */
/*   3n+3:4n+3      X-velocity                                                */
/*   4n+4:5n+4      Y-velocity                                                */
/*   5n+5:6n+5      Z-velocity                                                */
/*   6n+6:7n+6      X-acceleration                                            */
/*   7n+7:8n+7      Y-acceleration                                            */
/*   8n+8:9n+8      Z-acceleration                                            */
/*   9n+9           Final time                                                */
/******************************************************************************/
private:
  const Real target_;

  int getSize(const std::vector<Real> &u) const {
    return static_cast<int>(u.size()-10)/9;
  }

  void parseState(std::vector<Real> &x1, std::vector<Real> &x2, std::vector<Real> &x3,
                  std::vector<Real> &v1, std::vector<Real> &v2, std::vector<Real> &v3,
                  std::vector<Real> &a1, std::vector<Real> &a2, std::vector<Real> &a3,
                  Real &T, const std::vector<Real> &u) const {
    int n = getSize(u);
    T = u[9*n+9];
    for (int i = 0; i < n+1; ++i) {
      x1[i] = u[i];       x2[i] = u[n+1+i];   x3[i] = u[2*n+2+i];
      v1[i] = u[3*n+3+i]; v2[i] = u[4*n+4+i]; v3[i] = u[5*n+5+i];
      a1[i] = u[6*n+6+i]; a2[i] = u[7*n+7+i]; a3[i] = u[8*n+8+i];
    }
  }

  ROL::Ptr<const std::vector<Real> > getConstVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<std::vector<Real> > getVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
  }

public:
  PuttingObjective(const Real target = 1e-2) : target_(target) {}

  Real value( const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    // Unwrap u
    ROL::Ptr<const std::vector<Real> > up = getConstVector(u);
    // Get number of time steps
    int n = getSize(*up);
    // Compute final x-velocity
    Real vx = (*up)[4*n+3];
    // Compute final y-velocity
    Real vy = (*up)[5*n+4];
    // Return final speed
    return static_cast<Real>(0.25)*std::pow(vx*vx + vy*vy - target_,2);
    //return static_cast<Real>(0.5)*(vx*vx + vy*vy);
  }

  void gradient_1( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> >       gp = getVector(g);
    ROL::Ptr<const std::vector<Real> > up = getConstVector(u);
    // Get number of time steps
    int n = getSize(*up);
    // Compute final x-velocity
    Real vx = (*up)[4*n+3];
    // Compute final y-velocity
    Real vy = (*up)[5*n+4];
    // Return derivative of final speed
    (*gp)[4*n+3] = vx*(vx*vx + vy*vy - target_);
    (*gp)[5*n+4] = vy*(vx*vx + vy*vy - target_);
    //(*gp)[4*n+3] = vx;
    //(*gp)[5*n+4] = vy;
  }

  void gradient_2( ROL::Vector<Real> &g, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    g.zero();
  }

  void hessVec_11( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<std::vector<Real> >      hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real> > vp = getConstVector(v);
    ROL::Ptr<const std::vector<Real> > up = getConstVector(u);
    // Get number of time steps
    int n = getSize(*vp);
    // Return derivative of final speed
    Real two(2), three(3);
    Real vx = (*up)[4*n+3];
    Real vy = (*up)[5*n+4];
    (*hvp)[4*n+3] = (three*vx*vx + vy*vy - target_)*(*vp)[4*n+3] + two*vx*vy*(*vp)[5*n+4];
    (*hvp)[5*n+4] = two*vx*vy*(*vp)[4*n+3] + (vx*vx + three*vy*vy - target_)*(*vp)[5*n+4];
    // Return derivative of final speed
    //(*hvp)[4*n+3] = (*vp)[4*n+3];
    //(*hvp)[5*n+4] = (*vp)[5*n+4];
  }

  void hessVec_12( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_21( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_22( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, 
                   const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    hv.zero();
  }
};

template<class Real>
class GreenConstraint : public ROL::Constraint_SimOpt<Real> {
/******************************************************************************/
/* DATA LAYOUT FOR SIMULATION VARIABLE:                                       */
/*   0:n            X-position                                                */
/*   n+1:2n+1       Y-position                                                */
/*   2n+2:3n+2      Z-position                                                */
/*   3n+3:4n+3      X-velocity                                                */
/*   4n+4:5n+4      Y-velocity                                                */
/*   5n+5:6n+5      Z-velocity                                                */
/*   6n+6:7n+6      X-acceleration                                            */
/*   7n+7:8n+7      Y-acceleration                                            */
/*   8n+8:9n+8      Z-acceleration                                            */
/*   9n+9           Final time                                                */
/******************************************************************************/
private:
  int getSize(const std::vector<Real> &u) const {
    return static_cast<int>(u.size()-10)/9;
  }

  void parseState(std::vector<Real> &x1, std::vector<Real> &x2, std::vector<Real> &x3,
                  std::vector<Real> &v1, std::vector<Real> &v2, std::vector<Real> &v3,
                  std::vector<Real> &a1, std::vector<Real> &a2, std::vector<Real> &a3,
                  Real &T, const std::vector<Real> &u) const {
    int n = getSize(u);
    T = u[9*n+9];
    for (int i = 0; i < n+1; ++i) {
      x1[i] = u[i];       x2[i] = u[n+1+i];   x3[i] = u[2*n+2+i];
      v1[i] = u[3*n+3+i]; v2[i] = u[4*n+4+i]; v3[i] = u[5*n+5+i];
      a1[i] = u[6*n+6+i]; a2[i] = u[7*n+7+i]; a3[i] = u[8*n+8+i];
    }
  }

  ROL::Ptr<const std::vector<Real> > getConstVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<std::vector<Real> > getVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
  }

public:

  GreenConstraint(void) {}

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, 
                  const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> >       cp = getVector(c);
    ROL::Ptr<const std::vector<Real> > up = getConstVector(u);
    // Get number of time steps
    int n = getSize(*up);
    // Parse state vector
    Real T(0);
    std::vector<Real> x1(n+1), x2(n+1), x3(n+1);
    std::vector<Real> v1(n+1), v2(n+1), v3(n+1);
    std::vector<Real> a1(n+1), a2(n+1), a3(n+1);
    parseState(x1,x2,x3,v1,v2,v3,a1,a2,a3,T,*up);
    // Build residual
    for (int i = 0; i < n+1; ++i) {
      (*cp)[i] = x1[i]*x1[i] + x2[i]*x2[i];
    }
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                       const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> >      jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real> > vp = getConstVector(v);
    ROL::Ptr<const std::vector<Real> > up = getConstVector(u);
    // Get number of time steps
    int n = getSize(*up);
    // Parse state vector
    Real T(0);
    std::vector<Real> x1(n+1), x2(n+1), x3(n+1);
    std::vector<Real> v1(n+1), v2(n+1), v3(n+1);
    std::vector<Real> a1(n+1), a2(n+1), a3(n+1);
    parseState(x1,x2,x3,v1,v2,v3,a1,a2,a3,T,*up);
    // Parse direction vector
    Real vT(0);
    std::vector<Real> vx1(n+1), vx2(n+1), vx3(n+1);
    std::vector<Real> vv1(n+1), vv2(n+1), vv3(n+1);
    std::vector<Real> va1(n+1), va2(n+1), va3(n+1);
    parseState(vx1,vx2,vx3,vv1,vv2,vv3,va1,va2,va3,vT,*vp);
    // Build jacobian
    for (int i = 0; i < n+1; ++i) {
      (*jvp)[i] = static_cast<Real>(2)*(x1[i]*vx1[i] + x2[i]*vx2[i]);
    }
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                       const ROL::Vector<Real> &z, Real &tol) {
    jv.zero();
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, 
                              const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> >      jvp = getVector(ajv);
    ROL::Ptr<const std::vector<Real> > vp = getConstVector(v);
    ROL::Ptr<const std::vector<Real> > up = getConstVector(u);
    // Get number of time steps
    int n = getSize(*up);
    // Parse state vector
    Real T(0);
    std::vector<Real> x1(n+1), x2(n+1), x3(n+1);
    std::vector<Real> v1(n+1), v2(n+1), v3(n+1);
    std::vector<Real> a1(n+1), a2(n+1), a3(n+1);
    parseState(x1,x2,x3,v1,v2,v3,a1,a2,a3,T,*up);
    // Build adjoint jacobian
    for (int i = 0; i < n+1; ++i) {
      (*jvp)[i]       = static_cast<Real>(2)*x1[i]*(*vp)[i];
      (*jvp)[n+1+i]   = static_cast<Real>(2)*x2[i]*(*vp)[i];
    }
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u,
                              const ROL::Vector<Real> &z, Real &tol) {
    jv.zero();
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> >    ahwvp = getVector(ahwv);
    ROL::Ptr<const std::vector<Real> > wp = getConstVector(w);
    ROL::Ptr<const std::vector<Real> > vp = getConstVector(v);
    // Get number of time steps
    int n = getSize(*vp);
    // Parse direction vector
    Real vT(0);
    std::vector<Real> vx1(n+1), vx2(n+1), vx3(n+1);
    std::vector<Real> vv1(n+1), vv2(n+1), vv3(n+1);
    std::vector<Real> va1(n+1), va2(n+1), va3(n+1);
    parseState(vx1,vx2,vx3,vv1,vv2,vv3,va1,va2,va3,vT,*vp);
    // Build hessian
    for (int i = 0; i < n+1; ++i) {
      (*ahwvp)[i]     = static_cast<Real>(2)*(*wp)[i]*vx1[i];
      (*ahwvp)[n+1+i] = static_cast<Real>(2)*(*wp)[i]*vx2[i];
    }
  }
  
  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }
  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }
  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ahwv.zero();
  }
};
