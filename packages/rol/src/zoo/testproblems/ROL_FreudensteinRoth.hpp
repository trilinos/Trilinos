// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for Freudenstein and Roth's function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_FREUDENSTEINROTH_HPP
#define ROL_FREUDENSTEINROTH_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"

namespace ROL {
namespace ZOO {

/** \brief Freudenstein and Roth's function.
 */
template<class Real>
class Objective_FreudensteinRoth : public Objective<Real> {
public:
  Objective_FreudensteinRoth() {}

  Real value( const Vector<Real> &x, Real &tol ) {
    Ptr<const std::vector<Real> > ex
      = dynamic_cast<const StdVector<Real>&>(x).getVector();

    Real f1 = -13.0 + (*ex)[0] + ((5.0-(*ex)[1])*(*ex)[1] - 2.0)*(*ex)[1];
    Real f2 = -29.0 + (*ex)[0] + (((*ex)[1]+1.0)*(*ex)[1] - 14.0)*(*ex)[1];

    return f1*f1+f2*f2;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    Ptr<std::vector<Real> > eg
      = dynamic_cast<StdVector<Real>&>(g).getVector();
    Ptr<const std::vector<Real> > ex
      = dynamic_cast<const StdVector<Real>&>(x).getVector();

    Real f1 = -13.0 + (*ex)[0] + ((5.0-(*ex)[1])*(*ex)[1] - 2.0)*(*ex)[1];
    Real f2 = -29.0 + (*ex)[0] + (((*ex)[1]+1.0)*(*ex)[1] - 14.0)*(*ex)[1];

    Real f11 = 1.0;
    Real f12 = 10.0*(*ex)[1] - 3.0*(*ex)[1]*(*ex)[1] - 2.0;
    Real f21 = 1.0;
    Real f22 = 3.0*(*ex)[1]*(*ex)[1] + 2.0*(*ex)[1] - 14.0;

    (*eg)[0] = 2.0*(f11*f1 + f21*f2);
    (*eg)[1] = 2.0*(f12*f1 + f22*f2);
  }
#if USE_HESSVEC
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Ptr<std::vector<Real> > ehv
      = dynamic_cast<StdVector<Real>&>(hv).getVector();
    Ptr<const std::vector<Real> > ev
      = dynamic_cast<const StdVector<Real>&>(v).getVector();
    Ptr<const std::vector<Real> > ex
      = dynamic_cast<const StdVector<Real>&>(x).getVector();

    Real f1 = -13.0 + (*ex)[0] + ((5.0-(*ex)[1])*(*ex)[1] - 2.0)*(*ex)[1];
    Real f2 = -29.0 + (*ex)[0] + (((*ex)[1]+1.0)*(*ex)[1] - 14.0)*(*ex)[1];

    Real f11 = 1.0;
    Real f12 = 10.0*(*ex)[1] - 3.0*(*ex)[1]*(*ex)[1] - 2.0;
    Real f21 = 1.0;
    Real f22 = 3.0*(*ex)[1]*(*ex)[1] + 2.0*(*ex)[1] - 14.0;

    Real f122 = 10.0 - 6.0*(*ex)[1];
    Real f222 = 6.0*(*ex)[1] + 2.0;

    Real h11 = 2.0*(f11*f11) + 2.0*(f21*f21);
    Real h12 = 2.0*(f12*f11) + 2.0*(f22*f21);
    Real h22 = 2.0*(f122*f1 + f12*f12) + 2.0*(f222*f2 + f22*f22);

    (*ehv)[0] = h11*(*ev)[0] + h12*(*ev)[1];
    (*ehv)[1] = h12*(*ev)[0] + h22*(*ev)[1];
  }
#endif
  void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Ptr<std::vector<Real> > ehv
      = dynamic_cast<StdVector<Real>&>(hv).getVector();
    Ptr<const std::vector<Real> > ev
      = dynamic_cast<const StdVector<Real>&>(v).getVector();
    Ptr<const std::vector<Real> > ex
      = dynamic_cast<const StdVector<Real>&>(x).getVector();

    Real f1 = -13.0 + (*ex)[0] + ((5.0-(*ex)[1])*(*ex)[1] - 2.0)*(*ex)[1];
    Real f2 = -29.0 + (*ex)[0] + (((*ex)[1]+1.0)*(*ex)[1] - 14.0)*(*ex)[1];

    Real f11 = 1.0;
    Real f12 = 10.0*(*ex)[1] - 3.0*(*ex)[1]*(*ex)[1] - 2.0;
    Real f21 = 1.0;
    Real f22 = 3.0*(*ex)[1]*(*ex)[1] + 2.0*(*ex)[1] - 14.0;

    Real f122 = 10.0 - 6.0*(*ex)[1];
    Real f222 = 6.0*(*ex)[1] + 2.0;

    Real h11 = 2.0*(f11*f11) + 2.0*(f21*f21);
    Real h12 = 2.0*(f12*f11) + 2.0*(f22*f21);
    Real h22 = 2.0*(f122*f1 + f12*f12) + 2.0*(f222*f2 + f22*f22);

    (*ehv)[0] = (1.0/(h11*h22-h12*h12))*( h22*(*ev)[0] - h12*(*ev)[1]);
    (*ehv)[1] = (1.0/(h11*h22-h12*h12))*(-h12*(*ev)[0] + h11*(*ev)[1]);
  }
};

template<class Real>
class getFreudensteinRoth : public TestProblem<Real> {
public:
  getFreudensteinRoth(void){}

  Ptr<Objective<Real>> getObjective(void) const {
    // Instantiate Objective Function
    return makePtr<Objective_FreudensteinRoth<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem dimension
    int n = 2;
    // Get Initial Guess
    Ptr<std::vector<Real> > x0p = makePtr<std::vector<Real>>(n,0.0);
    (*x0p)[0] = 0.5; (*x0p)[1] = -2.0;
    return makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem dimension
    int n = 2;
    // Get Solution
    Ptr<std::vector<Real> > xp = makePtr<std::vector<Real>>(n,0.0);
    if (i == 0) {
      (*xp)[0] = 5.0; (*xp)[1] = 4.0;
    }
    else if (i == 1) {
      (*xp)[0] = 11.412779; (*xp)[1] = -0.896805;
    }
    else {
      throw Exception::NotImplemented(">>> ROL::FreudensteinRoth : The index i must be between 0 and 1!");
    }
    return makePtr<StdVector<Real>>(xp);
  }

  int getNumSolutions(void) const {
    return 2;
  }
};


} // End ZOO Namespace
} // End ROL Namespace

#endif
