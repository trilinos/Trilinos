// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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

