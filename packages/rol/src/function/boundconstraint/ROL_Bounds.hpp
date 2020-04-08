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

#ifndef ROL_BOUNDS_H
#define ROL_BOUNDS_H

#include "ROL_BoundConstraint.hpp"

/** @ingroup func_group
    \class ROL::Bounds
    \brief Provides the elementwise interface to apply upper and lower bound
           constraints.

*/

namespace ROL {

template<typename Real>
class Bounds : public BoundConstraint<Real> {
private:
  const Real scale_;
  const Real feasTol_;

  using BoundConstraint<Real>::lower_;
  using BoundConstraint<Real>::upper_;

  Ptr<Vector<Real>> mask_;

  Real min_diff_;

  Elementwise::ReductionMin<Real> minimum_;

  class Active : public Elementwise::BinaryFunction<Real> {
    public:
    Active(Real offset) : offset_(offset) {}
    Real apply( const Real &x, const Real &y ) const {
      return ((y <= offset_) ? 0 : x);
    }
    private:
    Real offset_;
  };

  class UpperBinding : public Elementwise::BinaryFunction<Real> {
    public:
    UpperBinding(Real xeps, Real geps) : xeps_(xeps), geps_(geps) {}
    Real apply( const Real &x, const Real &y ) const {
      return ((y < -geps_ && x <= xeps_) ? 0 : 1);
    }
    private:
    Real xeps_, geps_;
  };

  class LowerBinding : public Elementwise::BinaryFunction<Real> {
    public:
    LowerBinding(Real xeps, Real geps) : xeps_(xeps), geps_(geps) {}
    Real apply( const Real &x, const Real &y ) const {
      return ((y > geps_ && x <= xeps_) ? 0 : 1);
    }
    private:
    Real xeps_, geps_;
  };

  class PruneBinding : public Elementwise::BinaryFunction<Real> {
    public:
      Real apply( const Real &x, const Real &y ) const {
        return ((y == 1) ? x : 0);
      }
  } prune_;

public:

  Bounds(const Vector<Real> &x,
         bool isLower = true,
         Real scale = 1,
         Real feasTol = 1e-2);

  Bounds(const Ptr<Vector<Real>> &x_lo,
         const Ptr<Vector<Real>> &x_up,
         const Real scale = 1,
         const Real feasTol = 1e-2);

  void project( Vector<Real> &x ) override;

  void projectInterior( Vector<Real> &x ) override;

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) override;

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) override;

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) override;

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) override;

  bool isFeasible( const Vector<Real> &v ) override;

}; // class Bounds

} // namespace ROL

#include "ROL_Bounds_Def.hpp"

#endif
