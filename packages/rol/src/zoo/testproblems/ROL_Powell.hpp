// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for Powell's badly scaled function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_POWELL_HPP
#define ROL_POWELL_HPP

#include "ROL_ScaledStdVector.hpp"
#include "ROL_TestProblem.hpp"

namespace ROL {
namespace ZOO {

/** \brief Powell's badly scaled function.
 */
template<class Real>
class Objective_Powell : public Objective<Real> {

public:
  Objective_Powell() {}

  Real value( const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > xp
      = dynamic_cast<const PrimalScaledStdVector<Real>&>(x).getVector();

    Real f1 = 1.e4*(*xp)[0]*(*xp)[1] - 1.0;
    Real f2 = std::exp(-(*xp)[0]) + std::exp(-(*xp)[1]) - 1.0001;

    return f1*f1+f2*f2;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp
      = dynamic_cast<DualScaledStdVector<Real>&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > xp
      = dynamic_cast<const PrimalScaledStdVector<Real>&>(x).getVector();

    Real f1 = 1.e4*(*xp)[0]*(*xp)[1] - 1.0;
    Real f2 = std::exp(-(*xp)[0]) + std::exp(-(*xp)[1]) - 1.0001;

    Real f11 = 1.e4*(*xp)[1];
    Real f12 = 1.e4*(*xp)[0];
    Real f21 = -std::exp(-(*xp)[0]);
    Real f22 = -std::exp(-(*xp)[1]);

    (*gp)[0] = 2.0*(f11*f1 + f21*f2);
    (*gp)[1] = 2.0*(f12*f1 + f22*f2);
  }
#if USE_HESSVEC
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<std::vector<Real> > hvp
      = dynamic_cast<DualScaledStdVector<Real>&>(hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const PrimalScaledStdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > xp
      = dynamic_cast<const PrimalScaledStdVector<Real>&>(x).getVector();

    Real f1 = 1.e4*(*xp)[0]*(*xp)[1] - 1.0;
    Real f2 = std::exp(-(*xp)[0]) + std::exp(-(*xp)[1]) - 1.0001;

    Real f11 = 1.e4*(*xp)[1];
    Real f12 = 1.e4*(*xp)[0];
    Real f21 = -std::exp(-(*xp)[0]);
    Real f22 = -std::exp(-(*xp)[1]);

    Real f111 = 0.0;
    Real f112 = 1.e4;
    Real f121 = 1.e4;
    Real f122 = 0.0;
    Real f211 = std::exp(-(*xp)[0]);
    Real f212 = 0.0;
    Real f221 = 0.0;
    Real f222 = std::exp(-(*xp)[1]);

    Real h11 = 2.0*(f111*f1 + f11*f11) + 2.0*(f211*f2 + f21*f21);
    Real h12 = 2.0*(f112*f1 + f11*f12) + 2.0*(f212*f2 + f21*f22);
    Real h21 = 2.0*(f121*f1 + f12*f11) + 2.0*(f221*f2 + f22*f21);
    Real h22 = 2.0*(f122*f1 + f12*f12) + 2.0*(f222*f2 + f22*f22);

    (*hvp)[0] = h11*(*vp)[0] + h12*(*vp)[1];
    (*hvp)[1] = h21*(*vp)[0] + h22*(*vp)[1];
  }
#endif
  void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<std::vector<Real> > hvp
      = dynamic_cast<PrimalScaledStdVector<Real>&>(hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const DualScaledStdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > xp
      = dynamic_cast<const PrimalScaledStdVector<Real>&>(x).getVector();

    Real f1 = 1.e4*(*xp)[0]*(*xp)[1] - 1.0;
    Real f2 = std::exp(-(*xp)[0]) + std::exp(-(*xp)[1]) - 1.0001;

    Real f11 = 1.e4*(*xp)[1];
    Real f12 = 1.e4*(*xp)[0];
    Real f21 = -std::exp(-(*xp)[0]);
    Real f22 = -std::exp(-(*xp)[1]);

    Real f111 = 0.0;
    Real f112 = 1.e4;
    Real f121 = 1.e4;
    Real f122 = 0.0;
    Real f211 = std::exp(-(*xp)[0]);
    Real f212 = 0.0;
    Real f221 = 0.0;
    Real f222 = std::exp(-(*xp)[1]);

    Real h11 = 2.0*(f111*f1 + f11*f11) + 2.0*(f211*f2 + f21*f21);
    Real h12 = 2.0*(f112*f1 + f11*f12) + 2.0*(f212*f2 + f21*f22);
    Real h21 = 2.0*(f121*f1 + f21*f11) + 2.0*(f221*f2 + f22*f21);
    Real h22 = 2.0*(f122*f1 + f12*f12) + 2.0*(f222*f2 + f22*f22);

    (*hvp)[0] = (1.0/(h11*h22-h12*h21))*( h22*(*vp)[0] - h21*(*vp)[1]);
    (*hvp)[1] = (1.0/(h11*h22-h12*h21))*(-h12*(*vp)[0] + h11*(*vp)[1]);
  }
};

template<class Real>
class getPowell : public TestProblem<Real> {
private:
  ROL::Ptr<std::vector<Real> > scale_;

public:
  getPowell(void) {
    // Get problem scaling
    scale_ = ROL::makePtr<std::vector<Real>>(2,0.0);
    (*scale_)[0] = 1.e10; (*scale_)[1] = 1.e-2;
  }

  Ptr<Objective<Real>> getObjective(void) const {
    // Instantiate Objective Function
    return ROL::makePtr<Objective_Powell<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem dimension
    int n = 2;
    // Get Initial Guess
    ROL::Ptr<std::vector<Real> > x0p = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*x0p)[0] = 0.0; (*x0p)[1] = 1.0;
    return ROL::makePtr<PrimalScaledStdVector<Real>>(x0p,scale_);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem dimension
    int n = 2;
    // Get Solution: (*xp)[0] = 1.0981770261368074e-05;
    ROL::Ptr<std::vector<Real> > xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*xp)[0] = 1.098e-05; (*xp)[1] = 9.106;
    return ROL::makePtr<PrimalScaledStdVector<Real>>(xp,scale_);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
