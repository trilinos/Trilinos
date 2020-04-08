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

/** \file
    \brief  Contains definitions for std::vector bound constraints.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_STDBOUNDCONSTRAINT_HPP
#define ROL_STDBOUNDCONSTRAINT_HPP

#include "ROL_StdVector.hpp"
#include "ROL_BoundConstraint.hpp"

namespace ROL {

template<class Real>
class StdBoundConstraint : public BoundConstraint<Real> {
private:
  int dim_;
  std::vector<Real> x_lo_;
  std::vector<Real> x_up_;
  Real min_diff_;
  Real scale_;

  using BoundConstraint<Real>::lower_;
  using BoundConstraint<Real>::upper_;

public:
  StdBoundConstraint(std::vector<Real> &x, bool isLower = false, Real scale = Real(1));

  StdBoundConstraint(std::vector<Real> &l, std::vector<Real> &u, Real scale = Real(1));

  bool isFeasible( const Vector<Real> &x ) override;

  void project( Vector<Real> &x ) override;

  void projectInterior( Vector<Real> &x ) override;

  void pruneLowerActive(Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0)) override;

  void pruneUpperActive(Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0)) override;

  void pruneLowerActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) override;

  void pruneUpperActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) override;
};

}// End ROL Namespace

#include "ROL_StdBoundConstraint_Def.hpp"

#endif
