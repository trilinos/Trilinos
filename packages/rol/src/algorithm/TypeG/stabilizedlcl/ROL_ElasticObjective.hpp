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

#ifndef ROL_ELASTICOBJECTIVE_H
#define ROL_ELASTICOBJECTIVE_H

#include "ROL_AugmentedLagrangianObjective.hpp"
#include "ROL_PartitionedVector.hpp"

/** @ingroup func_group
    \class ROL::ElasticObjective
    \brief Provides the interface to evaluate the elastic augmented Lagrangian.

    This class implements the elastic augmented Lagrangian functional for use with
    ROL::StablizedLCLAlgorithm.  Given a function
    \f$f:\mathcal{X}\to\mathbb{R}\f$ and an equality constraint
    \f$c:\mathcal{X}\to\mathcal{C}\f$, the augmented Lagrangian functional is
    \f[
       L_A(x,\lambda,\mu) = f(x) +
           \langle \lambda, c(x)\rangle_{\mathcal{C}^*,\mathcal{C}} +
           \frac{\mu}{2} \langle \mathfrak{R}c(x),c(x)\rangle_{\mathcal{C}^*,\mathcal{C}}
           + \sigma\langle \mathfrak{R} e, u-v\ranlge_{\mathcal{C}^*,\mathcal{C}}
    \f]
    where \f$\lambda\in\mathcal{C}^*\f$ denotes the Lagrange multiplier estimate,
    \f$\mu > 0\f$ and \f$\sigma>0\f$ are penalty parameters,
    \f$e\in\mathcal{C}\f$ is the constant one vector, and
    \f$\mathfrak{R}\in\mathcal{L}(\mathcal{C},\mathcal{C}^*)\f$ is the Riesz operator
    on the constraint space.

    This implementation permits the scaling of \f$L_A\f$ by \f$\mu^{-1}\f$ and also
    permits the Hessian approximation
    \f[
        \nabla^2_x L_A(x,\lambda,\mu)v \approx \nabla^2 f(x) v + \mu c'(x)^*\mathfrak{R} c'(x)v.
    \f]

    ---
*/


namespace ROL {

template<typename Real>
class ElasticObjective : public Objective<Real> {
private:
  // Required for Augmented Lagrangian definition
  Ptr<AugmentedLagrangianObjective<Real>> alobj_;
  Ptr<Vector<Real>> e_, tmp_;
  Real sigma_, cscale_;

public:
  ElasticObjective(const Ptr<Objective<Real>> &obj,
                   const Ptr<Constraint<Real>> &con,
                   const Real penaltyParameter,
                   const Real sigma,
                   const Vector<Real> &dualOptVec,
                   const Vector<Real> &primConVec,
                   const Vector<Real> &dualConVec,
                   ParameterList &parlist);

  ElasticObjective(const Ptr<Objective<Real>> &obj,
                   const Ptr<Constraint<Real>> &con,
                   const Real penaltyParameter,
                   const Real sigma,
                   const Vector<Real> &dualOptVec,
                   const Vector<Real> &primConVec,
                   const Vector<Real> &dualConVec,
                   const bool scaleLagrangian,
                   const int HessianApprox);

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

  // Set problem data scalings
  void setScaling(const Real fscale = 1.0, const Real cscale = 1.0);
  // Return objective function value
  Real getObjectiveValue(const Vector<Real> &x, Real &tol);
  // Compute objective function gradient
  const Ptr<const Vector<Real>> getObjectiveGradient(const Vector<Real> &x, Real &tol);
  // Return constraint value
  const Ptr<const Vector<Real>> getConstraintVec(const Vector<Real> &x, Real &tol);
  // Return total number of constraint evaluations
  int getNumberConstraintEvaluations(void) const;
  // Return total number of objective evaluations
  int getNumberFunctionEvaluations(void) const;
  // Return total number of gradient evaluations
  int getNumberGradientEvaluations(void) const;
  // Reset with upated penalty parameter
  void reset(const Vector<Real> &multiplier, Real penaltyParameter, Real sigma);
  // Return augmented Lagrangian
  const Ptr<AugmentedLagrangianObjective<Real>> getAugmentedLagrangian(void) const;
}; // class ElasticObjective

} // namespace ROL

#include "ROL_ElasticObjective_Def.hpp"

#endif
