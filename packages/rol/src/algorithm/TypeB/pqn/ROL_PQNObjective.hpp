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

#ifndef ROL_PQNOBJECTIVE_H
#define ROL_PQNOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Secant.hpp"

/** @ingroup func_group
    \class ROL::PQNObjective
    \brief Provides the interface to evaluate the quadratic quasi-Newton objective.

    This class implements the PQN quasi-Newton objective for use with
    ROL::TypeB::QuasiNewtonAlgorithm.  Given a function
    \f$f:\mathcal{X}\to\mathbb{R}\f$ and a Hessian approximation
    \f$B_k:\mathcal{X}}\to\mathcal{X}^*\f$, the functional is
    \f[
       q_k(x) = \frac{1}{2}\langle B_k (x-x_k),(x-x_k)\rangle_{\mathcal{X}^*,\mathcal{X}}
                + \langle f'(x_k), (x-x_k)\rangle_{\mathcal{X}^*,\mathcal{X}}.
    \f]

    ---
*/


namespace ROL {

template<typename Real>
class PQNObjective : public Objective<Real> {
private:
  const Ptr<Secant<Real>> secant_;
  const Ptr<Vector<Real>> x_, g_, pwa_, dwa_;

public:
  PQNObjective(const Ptr<Secant<Real>> &secant,
               const Vector<Real> &x,
               const Vector<Real> &g);

  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

  void setAnchor(const Vector<Real> &x, const Vector<Real> &g);
}; // class PQNObjective

} // namespace ROL

#include "ROL_PQNObjective_Def.hpp"

#endif
