// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for Beale's function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_BEALE_HPP
#define ROL_BEALE_HPP

#include "ROL_ScaledStdVector.hpp"
#include "ROL_TestProblem.hpp"

namespace ROL {
namespace ZOO {

  /** \brief Beale's function.
   */
  template<class Real>
  class Objective_Beale : public Objective<Real> {
  private:
    std::vector<Real> y_;

  public:
    Objective_Beale() {
      y_.clear();
      y_.push_back(static_cast<Real>(1.5));
      y_.push_back(static_cast<Real>(2.25));
      y_.push_back(static_cast<Real>(2.625));
    }

    Real value( const Vector<Real> &x, Real &tol ) {
      Ptr<const std::vector<Real> > ex
        = dynamic_cast<const StdVector<Real>&>(x).getVector();

      Real f1 = static_cast<Real>(1.5)-(*ex)[0]*(static_cast<Real>(1)-(*ex)[1]);
      Real f2 = static_cast<Real>(2.25)-(*ex)[0]*(static_cast<Real>(1)-pow((*ex)[1],2));
      Real f3 = static_cast<Real>(2.625)-(*ex)[0]*(static_cast<Real>(1)-pow((*ex)[1],3));

      return pow(f1,2)+pow(f2,2)+pow(f3,2);
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      Ptr<std::vector<Real> > eg
        = dynamic_cast<StdVector<Real>&>(g).getVector();
      Ptr<const std::vector<Real> > ex
        = dynamic_cast<const StdVector<Real>&>(x).getVector();

      Real f1 = static_cast<Real>(1.5)-(*ex)[0]*(static_cast<Real>(1)-(*ex)[1]);
      Real f2 = static_cast<Real>(2.25)-(*ex)[0]*(static_cast<Real>(1)-pow((*ex)[1],2));
      Real f3 = static_cast<Real>(2.625)-(*ex)[0]*(static_cast<Real>(1)-pow((*ex)[1],3));
      Real df1dx = -(static_cast<Real>(1)-(*ex)[1]);
      Real df1dy = (*ex)[0];
      Real df2dx = -(static_cast<Real>(1)-pow((*ex)[1],2));
      Real df2dy = static_cast<Real>(2)*(*ex)[0]*(*ex)[1];
      Real df3dx = -(static_cast<Real>(1)-pow((*ex)[1],3));
      Real df3dy = static_cast<Real>(3)*(*ex)[0]*pow((*ex)[1],2);

      (*eg)[0] = static_cast<Real>(2)*df1dx*f1+static_cast<Real>(2)*df2dx*f2+static_cast<Real>(2)*df3dx*f3;
      (*eg)[1] = static_cast<Real>(2)*df1dy*f1+static_cast<Real>(2)*df2dy*f2+static_cast<Real>(2)*df3dy*f3;
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Ptr<std::vector<Real> > ehv
        = dynamic_cast<StdVector<Real>&>(hv).getVector();
      Ptr<const std::vector<Real> > ev
        = dynamic_cast<const StdVector<Real>&>(v).getVector();
      Ptr<const std::vector<Real> > ex
        = dynamic_cast<const StdVector<Real>&>(x).getVector();

      Real f1 = static_cast<Real>(1.5)-(*ex)[0]*(static_cast<Real>(1)-(*ex)[1]);
      Real f2 = static_cast<Real>(2.25)-(*ex)[0]*(static_cast<Real>(1)-pow((*ex)[1],2));
      Real f3 = static_cast<Real>(2.625)-(*ex)[0]*(static_cast<Real>(1)-pow((*ex)[1],3));
      Real df1dx = -(static_cast<Real>(1)-(*ex)[1]);
      Real df1dy = (*ex)[0];
      Real df2dx = -(static_cast<Real>(1)-pow((*ex)[1],2));
      Real df2dy = static_cast<Real>(2)*(*ex)[0]*(*ex)[1];
      Real df3dx = -(static_cast<Real>(1)-pow((*ex)[1],3));
      Real df3dy = static_cast<Real>(3)*(*ex)[0]*pow((*ex)[1],2);
      Real d2f1dx2 = static_cast<Real>(0);
      Real d2f1dy2 = static_cast<Real>(0);
      Real d2f1dxdy = static_cast<Real>(1);
      Real d2f2dx2 = static_cast<Real>(0);
      Real d2f2dy2 = static_cast<Real>(2)*(*ex)[0];
      Real d2f2dxdy = static_cast<Real>(2)*(*ex)[1];
      Real d2f3dx2 = static_cast<Real>(0);
      Real d2f3dy2 = static_cast<Real>(6)*(*ex)[0]*(*ex)[1];
      Real d2f3dxdy = static_cast<Real>(3)*pow((*ex)[1],2);

      Real H11 = static_cast<Real>(2)*(d2f1dx2*f1+df1dx*df1dx)+static_cast<Real>(2)*(d2f2dx2*f2+df2dx*df2dx)
                  +static_cast<Real>(2)*(d2f3dx2*f3+df3dx*df3dx);
      Real H22 = static_cast<Real>(2)*(d2f1dy2*f1+df1dy*df1dy)+static_cast<Real>(2)*(d2f2dy2*f2+df2dy*df2dy)
                  +static_cast<Real>(2)*(d2f3dy2*f3+df3dy*df3dy);
      Real H12 = static_cast<Real>(2)*(d2f1dxdy*f1 + df1dx*df1dy)+static_cast<Real>(2)*(d2f2dxdy*f2 + df2dx*df2dy)
                  +static_cast<Real>(2)*(d2f3dxdy*f3 + df3dx*df3dy);

      (*ehv)[0] = H11*(*ev)[0]+H12*(*ev)[1];
      (*ehv)[1] = H12*(*ev)[0]+H22*(*ev)[1];
    }
#endif
    void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Ptr<std::vector<Real> > ehv
        = dynamic_cast<StdVector<Real>&>(hv).getVector();
      Ptr<const std::vector<Real> > ev
        = dynamic_cast<const StdVector<Real>&>(v).getVector();
      Ptr<const std::vector<Real> > ex
        = dynamic_cast<const StdVector<Real>&>(x).getVector();

      Real f1 = static_cast<Real>(1.5)-(*ex)[0]*(static_cast<Real>(1)-(*ex)[1]);
      Real f2 = static_cast<Real>(2.25)-(*ex)[0]*(static_cast<Real>(1)-pow((*ex)[1],2));
      Real f3 = static_cast<Real>(2.625)-(*ex)[0]*(static_cast<Real>(1)-pow((*ex)[1],3));
      Real df1dx = -(static_cast<Real>(1)-(*ex)[1]);
      Real df1dy = (*ex)[0];
      Real df2dx = -(static_cast<Real>(1)-pow((*ex)[1],2));
      Real df2dy = static_cast<Real>(2)*(*ex)[0]*(*ex)[1];
      Real df3dx = -(static_cast<Real>(1)-pow((*ex)[1],3));
      Real df3dy = static_cast<Real>(3)*(*ex)[0]*pow((*ex)[1],2);
      Real d2f1dx2 = static_cast<Real>(0);
      Real d2f1dy2 = static_cast<Real>(0);
      Real d2f1dxdy = static_cast<Real>(1);
      Real d2f2dx2 = static_cast<Real>(0);
      Real d2f2dy2 = static_cast<Real>(2)*(*ex)[0];
      Real d2f2dxdy = static_cast<Real>(2)*(*ex)[1];
      Real d2f3dx2 = static_cast<Real>(0);
      Real d2f3dy2 = static_cast<Real>(6)*(*ex)[0]*(*ex)[1];
      Real d2f3dxdy = static_cast<Real>(3)*pow((*ex)[1],2);

      Real H11 = static_cast<Real>(2)*(d2f1dx2*f1+df1dx*df1dx)+static_cast<Real>(2)*(d2f2dx2*f2+df2dx*df2dx)
                  +static_cast<Real>(2)*(d2f3dx2*f3+df3dx*df3dx);
      Real H22 = static_cast<Real>(2)*(d2f1dy2*f1+df1dy*df1dy)+static_cast<Real>(2)*(d2f2dy2*f2+df2dy*df2dy)
                  +static_cast<Real>(2)*(d2f3dy2*f3+df3dy*df3dy);
      Real H12 = static_cast<Real>(2)*(d2f1dxdy*f1 + df1dx*df1dy)+static_cast<Real>(2)*(d2f2dxdy*f2 + df2dx*df2dy)
                  +static_cast<Real>(2)*(d2f3dxdy*f3 + df3dx*df3dy);

      (*ehv)[0] = (static_cast<Real>(1)/(H11*H22-H12*H12))*( H22*(*ev)[0] - H12*(*ev)[1]);
      (*ehv)[1] = (static_cast<Real>(1)/(H11*H22-H12*H12))*(-H12*(*ev)[0] + H11*(*ev)[1]);
    }
  };

  template<class Real>
  class getBeale : public TestProblem<Real> {
  public:
    getBeale(void) {}

    Ptr<Objective<Real>> getObjective(void) const {
      return makePtr<Objective_Beale<Real>>();
    }

    Ptr<Vector<Real>> getInitialGuess(void) const {
      int n = 2;
      // Build scale
      Ptr<std::vector<Real> > scale = makePtr<std::vector<Real>>(n,static_cast<Real>(0));
      (*scale)[0] = static_cast<Real>(1.e-1);
      (*scale)[1] = static_cast<Real>(1.e1);

      // Get Initial Guess
      Ptr<std::vector<Real> > x0p = makePtr<std::vector<Real>>(n,static_cast<Real>(0)); 
      (*x0p)[0] = static_cast<Real>(1);
      (*x0p)[1] = static_cast<Real>(1);

      return makePtr<PrimalScaledStdVector<Real>>(x0p,scale);
    }

    Ptr<Vector<Real>> getSolution(const int i = 0) const {
      int n = 2;
      // Build scale
      Ptr<std::vector<Real> > scale = makePtr<std::vector<Real>>(n,static_cast<Real>(0));
      (*scale)[0] = static_cast<Real>(1.e-1);
      (*scale)[1] = static_cast<Real>(1.e1);

      // Get Solution
      Ptr<std::vector<Real> > xp  = makePtr<std::vector<Real>>(n,static_cast<Real>(0));
      (*xp)[0] = static_cast<Real>(3);
      (*xp)[1] = static_cast<Real>(0.5);

      return makePtr<PrimalScaledStdVector<Real>>(xp,scale);
    }
  };

}// End ZOO Namespace
}// End ROL Namespace

#endif
