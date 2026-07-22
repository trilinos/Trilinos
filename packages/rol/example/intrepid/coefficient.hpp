// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef COEFFICIENT_HPP
#define COEFFICIENT_HPP

#include "discretization.hpp"
#include <vector>

/* Abstract BVP Coefficient class for problems of the form
 *
 * -[(diff(u,z,x)u'] + advec(u,z,x)u' + react(u,z,x) = 0
 *
 *  with homogeneous Neumann conditions
 *
 */
template<class Real, class ScalarT>
class Coefficient {

  typedef std::vector<Real>                 vec;
  typedef Intrepid::FieldContainer<Real>    FC;
  typedef Intrepid::FieldContainer<ScalarT> FCT;

  public:
 
    virtual ~Coefficient() {}

    virtual void reaction(FCT &result, const FC &x, const FCT &u, const FCT &z, const vec &param)  = 0;
 
    virtual void advection(FCT &result, const FC &x, const FCT &u, const FCT &z, const vec &param) = 0;

    virtual void diffusion(FCT &result, const FC &x, const FCT &u, const FCT &z, const vec &param) = 0;

};


// A particular instance of a coefficient class where
//
// diff(u,z,x) = 1
// advec(u,z,x) = 0
// react(u,z,x) = u^2*z^2-1

template<class Real, class ScalarT>
class ExampleCoefficient : public Coefficient<Real,ScalarT> {

 
  typedef std::vector<Real>                   vec;
  typedef Intrepid::FieldContainer<Real>      FC;
  typedef Intrepid::FieldContainer<ScalarT>   FCT;

  public:

    ExampleCoefficient() {}
    ~ExampleCoefficient() {}

    void reaction(FCT &result, const FC &x, const FCT &u, const FCT &z, const vec &param) {

      int numCells  = x.dimension(0);
      int numCubPts = x.dimension(1);

      for(int c=0; c<numCells; ++c) {
        for(int p=0; p<numCubPts; ++p) {
          result(c,p) = z(c,p)*z(c,p)*u(c,p)*u(c,p)-1.0;
        }
      }
    }

    void advection(FCT &result, const FC &x, const FCT &u, const FCT &z, const vec &param) {

      int numCells  = x.dimension(0);
      int numCubPts = x.dimension(1);
      int spaceDim  = x.dimension(2);

      for(int c=0; c<numCells; ++c) {
        for(int p=0; p<numCubPts; ++p) {
          for(int d=0; d<spaceDim; ++d) {
            result(c,p,d) = ScalarT(0.0);
          }
        }
      }
    }

    void diffusion(FCT &result, const FC &x, const FCT &u, const FCT &z, const vec &param) {

      int numCells  = x.dimension(0);
      int numCubPts = x.dimension(1);

      for(int c=0; c<numCells; ++c) {
        for(int p=0; p<numCubPts; ++p) {
          result(c,p) = ScalarT(1.0);
        }
      }
    }

};


#endif // COEFFICIENT

