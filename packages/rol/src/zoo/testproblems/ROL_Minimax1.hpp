// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MINIMAX1_HPP
#define ROL_MINIMAX1_HPP

#include "ROL_TestProblem.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Ptr.hpp"

namespace ROL {
namespace ZOO {

template<class Real>
class Minimax1 : public Objective<Real> {

  typedef std::vector<Real>  vector;
  typedef Vector<Real>       V;
  typedef StdVector<Real>    SV;  

private:

  Ptr<const vector> getVector( const V& x ) {
    return dynamic_cast<const SV&>(x).getVector();
  }

  Ptr<vector> getVector( V& x ) {
    return dynamic_cast<SV&>(x).getVector();
  }

public:
  Minimax1(void) {}

  Real value(const Vector<Real> &x, Real &tol) {
    Ptr<const vector> xp = getVector(x);
    Real f1 = std::pow((*xp)[0],2.0) + std::pow((*xp)[1],4.0);
    Real f2 = std::pow(2.0-(*xp)[0],2.0) + std::pow(2.0-(*xp)[1],2.0);
    Real f3 = 2.0*std::exp(-(*xp)[0] + (*xp)[1]);
    return std::max(f1,std::max(f2,f3));
  }

  void gradient(Vector<Real> &g, const Vector<Real> &x, Real &tol) {
    Ptr<const vector> xp = getVector(x);
    Ptr<vector> gp = getVector(g);
    Real f1 = std::pow((*xp)[0],2.0) + std::pow((*xp)[1],4.0);
    Real f2 = std::pow(2.0-(*xp)[0],2.0) + std::pow(2.0-(*xp)[1],2.0);
    Real f3 = 2.0*std::exp(-(*xp)[0] + (*xp)[1]);

    //int ind = ((f1 >= f2) ? ((f1 < f3) ? 3 : 1) : ((f2 < f3) ? 3 : 2));
    if( f1 >= std::max(f2,f3) ) {
      (*gp)[0] = 2.0*(*xp)[0];
      (*gp)[1] = 4.0*std::pow((*xp)[1],3.0);
    }
    else if ( f2 >= std::max(f1,f3) ) {
      (*gp)[0] = 2.0*(*xp)[0]-4.0;
      (*gp)[1] = 2.0*(*xp)[1]-4.0;
    }
    else if ( f3 >= std::max(f1,f2) ) {
      (*gp)[0] = -2.0*std::exp(-(*xp)[0]+(*xp)[1]);
      (*gp)[1] = 2.0*std::exp(-(*xp)[0]+(*xp)[1]);
    } 
  }
}; // class Minimax1

template<class Real>
class getMinimax1 : public TestProblem<Real> {
public:
  getMinimax1(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    return makePtr<Minimax1<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    Ptr<std::vector<Real> > x_ptr = makePtr<std::vector<Real>>(2, 0.0);
    (*x_ptr)[0] = 1.0; (*x_ptr)[1] = -0.1;
    return makePtr<StdVector<Real>>(x_ptr);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    Ptr<std::vector<Real> > z_ptr = makePtr<std::vector<Real>>(2, 0.0);
    (*z_ptr)[0] = 1.13904; (*z_ptr)[1] = 0.89956;
    return makePtr<StdVector<Real>>(z_ptr);
  }
};

} // namespace ZOO
} // namespace ROL

#endif
