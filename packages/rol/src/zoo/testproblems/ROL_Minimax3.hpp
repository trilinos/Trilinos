// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MINIMAX3_HPP
#define ROL_MINIMAX3_HPP

#include "ROL_TestProblem.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Ptr.hpp"

namespace ROL {
namespace ZOO {

template<class Real>
class Minimax3 : public Objective<Real> {

  typedef std::vector<Real>  vector;
  typedef Vector<Real>       V;
  typedef StdVector<Real>    SV;  

private:

  ROL::Ptr<const vector> getVector( const V& x ) {
    
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::Ptr<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

public:
  Minimax3(void) {}

  Real value(const Vector<Real> &x, Real &tol) {
   
    
    ROL::Ptr<const vector> xp = getVector(x);

    Real F  = std::pow((*xp)[0],2.0) + std::pow((*xp)[1],2.0) + 2.0*std::pow((*xp)[2],2.0) 
            + std::pow((*xp)[3],2.0) - 5.0*(*xp)[0] - 5.0*(*xp)[1] - 21.0*(*xp)[2] + 7.0*(*xp)[3];
    Real g2 = -std::pow((*xp)[0],2.0)-std::pow((*xp)[1],2.0)-std::pow((*xp)[2],2.0)-std::pow((*xp)[3],2.0)
              -(*xp)[0]+(*xp)[1]-(*xp)[2]+(*xp)[3]+8.0;
    Real g3 = -std::pow((*xp)[0],2.0)-2.0*std::pow((*xp)[1],2.0)-std::pow((*xp)[2],2.0)
              -2.0*std::pow((*xp)[3],2.0)+(*xp)[0]+(*xp)[3]+10.0;
    Real g4 = -std::pow((*xp)[0],2.0)-std::pow((*xp)[1],2.0)-std::pow((*xp)[2],2.0)
              -2.0*(*xp)[0]+(*xp)[1]+(*xp)[3]+5.0;
    Real a2 = 10.0, a3 = 10.0, a4 = 10.0;
    return std::max(F,std::max(F-a2*g2,std::max(F-a3*g3,F-a4*g4)));
  }

  void gradient(Vector<Real> &g, const Vector<Real> &x, Real &tol) {

    
    ROL::Ptr<const vector> xp = getVector(x);
    ROL::Ptr<vector> gp = getVector(g); 

    Real F  = std::pow((*xp)[0],2.0) + std::pow((*xp)[1],2.0) + 2.0*std::pow((*xp)[2],2.0) 
            + std::pow((*xp)[3],2.0) - 5.0*(*xp)[0] - 5.0*(*xp)[1] - 21.0*(*xp)[2] + 7.0*(*xp)[3];
    Real g2 = -std::pow((*xp)[0],2.0)-std::pow((*xp)[1],2.0)-std::pow((*xp)[2],2.0)-std::pow((*xp)[3],2.0)
              -(*xp)[0]+(*xp)[1]-(*xp)[2]+(*xp)[3]+8.0;
    Real g3 = -std::pow((*xp)[0],2.0)-2.0*std::pow((*xp)[1],2.0)-std::pow((*xp)[2],2.0)
              -2.0*std::pow((*xp)[3],2.0)+(*xp)[0]+(*xp)[3]+10.0;
    Real g4 = -std::pow((*xp)[0],2.0)-std::pow((*xp)[1],2.0)-std::pow((*xp)[2],2.0)
              -2.0*(*xp)[0]+(*xp)[1]+(*xp)[3]+5.0;
    Real a2 = 10.0, a3 = 10.0, a4 = 10.0;

    (*gp)[0] = 2.0*(*xp)[0] - 5.0;
    (*gp)[1] = 2.0*(*xp)[1] - 5.0;
    (*gp)[2] = 4.0*(*xp)[2] - 21.0;
    (*gp)[3] = 2.0*(*xp)[3] + 7.0;
    if ( F-a2*g2 >= std::max(F,std::max(F-a3*g3,F-a4*g4)) ) {
      (*gp)[0] += a2*(2.0*(*xp)[0] + 1.0);
      (*gp)[1] += a2*(2.0*(*xp)[1] - 1.0);
      (*gp)[2] += a2*(2.0*(*xp)[2] + 1.0);
      (*gp)[3] += a2*(2.0*(*xp)[3] - 1.0);
    }
    else if ( F-a3*g3 >= std::max(F,std::max(F-a2*g2,F-a4*g4)) ) {
      (*gp)[0] += a2*(2.0*(*xp)[0] - 1.0);
      (*gp)[1] += a2*(4.0*(*xp)[1]);
      (*gp)[2] += a2*(2.0*(*xp)[2]);
      (*gp)[3] += a2*(4.0*(*xp)[3] - 1.0);
    } 
    else if ( F-a4*g4 >= std::max(F,std::max(F-a2*g2,F-a3*g3)) ) {
      (*gp)[0] += a2*(2.0*(*xp)[0] + 2.0);
      (*gp)[1] += a2*(2.0*(*xp)[1] - 1.0);
      (*gp)[2] += a2*(2.0*(*xp)[2]);
      (*gp)[3] += a2*(-1.0);
    } 
  }
}; // class Minimax3

template<class Real>
class getMinimax3 : public TestProblem<Real> {
public:
  getMinimax3(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    return makePtr<Minimax3<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    ROL::Ptr<std::vector<Real> > x_ptr = ROL::makePtr<std::vector<Real>>(4, 0.0);
    return makePtr<StdVector<Real>>(x_ptr);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    ROL::Ptr<std::vector<Real> > z_ptr = ROL::makePtr<std::vector<Real>>(4, 0.0);
    (*z_ptr)[0] = 0.0; (*z_ptr)[1] = 1.0;
    (*z_ptr)[2] = 2.0; (*z_ptr)[3] = -1.0;
    return makePtr<StdVector<Real>>(z_ptr);
  }
};

} // namespace ZOO
} // namespace ROL

#endif
