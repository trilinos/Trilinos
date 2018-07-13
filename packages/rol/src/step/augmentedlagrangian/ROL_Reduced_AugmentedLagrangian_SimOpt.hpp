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

#ifndef ROL_REDUCED_AUGMENTEDLAGRANGIAN_SIMOPT_H
#define ROL_REDUCED_AUGMENTEDLAGRANGIAN_SIMOPT_H

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_AugmentedLagrangian_SimOpt.hpp"
#include "ROL_AugmentedLagrangian.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Ptr.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::Reduced_AugmentedLagrangian_SimOpt
    \brief Provides the interface to evaluate the reduced SimOpt augmented Lagrangian.

    This class implements the reduced SimOpt augmented Lagrangian functional for use with
    ROL::AugmentedLagrangianStep.  Given a function
    \f$f:\mathcal{U}\times\mathcal{Z}\to\mathbb{R}\f$, a reducible equality constaint
    \f$c_r:\mathcal{U}\times\mathcal{Z}\to\mathcal{C}_r\f$ and another equality constraint
    \f$c_a:\mathcal{U}\times\mathcal{Z}\to\mathcal{C}_a\f$, the (partially)
    augmented Lagrangian functional is
    \f[
       L_A(u,z,\lambda,\mu) = f(u,z) +
           \langle \lambda, c_a(u,z)\rangle_{\mathcal{C}_a^*,\mathcal{C}_a} +
           \frac{\mu}{2} \langle \mathfrak{R}c_a(u,z),c_a(u,z)\rangle_{\mathcal{C}_a^*,\mathcal{C}_a}
    \f]
    where \f$\lambda\in\mathcal{C}_a^*\f$ denotes the Lagrange multiplier estimate,
    \f$\mu > 0\f$ is the penalty parameter and
    \f$\mathfrak{R}\in\mathcal{L}(\mathcal{C}_a,\mathcal{C}_a^*)\f$ is the Riesz operator
    on the constraint space \f$\mathcal{C}_a\f$.  Since \f$c_r\f$ is reducible, there exists
    a solution operator \f$S:\mathcal{Z}\to\mathcal{U}\f$ such that
    \f[
       c_r(S(z),z) = 0 \quad\forall\, z\in\mathcal{Z}.
    \f]
    Substituting \f$S(z)\f$ into \f$L_A\f$ yields the reduced augmented Lagrangian
    \f$\bar{L}_A(z,\lambda,\mu) = L_A(S(z),z,\lambda,\mu)\f$.

    This implementation permits the scaling of \f$L_A\f$ by \f$\mu^{-1}\f$ and also
    permits the Hessian approximation
    \f[
        \nabla^2_z \bar{L}_A(z,\lambda,\mu)v \approx
           \nabla^2 \bar{f}(z) v + \mu \bar{c}_a'(z)^*\mathfrak{R} \bar{c}_a'(z)v
    \f]
    where \f$\bar{f}(z) = f(S(z),z)\f$ and \f$\bar{c}_a(z) = c_a(S(z),z)\f$.

    ---
*/


namespace ROL {

template <class Real>
class Reduced_AugmentedLagrangian_SimOpt : public AugmentedLagrangian<Real> {
private:
  ROL::Ptr<AugmentedLagrangian_SimOpt<Real> > augLagSimOpt_;
  ROL::Ptr<Reduced_Objective_SimOpt<Real> > rAugLagSimOpt_;
  ROL::Ptr<Vector<Real> > state_;

  // Evaluation counters
  int ngval_;

public:
  Reduced_AugmentedLagrangian_SimOpt(const ROL::Ptr<Objective_SimOpt<Real> > &obj,
                                     const ROL::Ptr<Constraint_SimOpt<Real> > &redCon,
                                     const ROL::Ptr<Constraint_SimOpt<Real> > &augCon,
                                     const ROL::Ptr<Vector<Real> > &state,
                                     const ROL::Ptr<Vector<Real> > &control,
                                     const ROL::Ptr<Vector<Real> > &adjoint,
                                     const ROL::Ptr<Vector<Real> > &augConVec,
                                     const ROL::Ptr<Vector<Real> > &multiplier,
                                     const Real penaltyParameter,
                                     ROL::ParameterList &parlist) : state_(state),
                                     ngval_(0) {

    augLagSimOpt_ = ROL::makePtr<AugmentedLagrangian_SimOpt<Real>>(obj,
                                                                      augCon,
                                                                      *multiplier,
                                                                      penaltyParameter,
                                                                      *state,
                                                                      *control,
                                                                      *augConVec,
                                                                      parlist);
    rAugLagSimOpt_ = ROL::makePtr<Reduced_Objective_SimOpt<Real>>(augLagSimOpt_,redCon,state,control,adjoint);
    rAugLagSimOpt_->update(*control);
    Real tol = 1e-8;
    rAugLagSimOpt_->value(*control,tol);
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    rAugLagSimOpt_->update(x,flag,iter);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    return rAugLagSimOpt_->value(x,tol);
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    ngval_++;
    rAugLagSimOpt_->gradient(g,x,tol);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    rAugLagSimOpt_->hessVec(hv,v,x,tol);
  }

  // Return objective function value
  Real getObjectiveValue(const Vector<Real> &x) {
    return augLagSimOpt_->getObjectiveValue(*state_,x);
  }

  // Return constraint value
  void getConstraintVec(Vector<Real> &c, const Vector<Real> &x) {
    augLagSimOpt_->getConstraintVec(c,*state_,x);
  }

  // Return total number of constraint evaluations
  int getNumberConstraintEvaluations(void) const {
    return augLagSimOpt_->getNumberConstraintEvaluations();
  }

  // Return total number of objective evaluations
  int getNumberFunctionEvaluations(void) const {
    return augLagSimOpt_->getNumberFunctionEvaluations();
  }

  // Return total number of gradient evaluations
  int getNumberGradientEvaluations(void) const {
    return ngval_;
    //return augLagSimOpt_->getNumberGradientEvaluations();
  }

  // Reset with upated penalty parameter
  void reset(const Vector<Real> &multiplier, const Real penaltyParameter) {
    ngval_ = 0;
    augLagSimOpt_->reset(multiplier,penaltyParameter);
  }

// For parametrized (stochastic) objective functions and constraints
public:
  void setParameter(const std::vector<Real> &param) {
    AugmentedLagrangian<Real>::setParameter(param);
    rAugLagSimOpt_->setParameter(param);
  }
}; // class AugmentedLagrangian

} // namespace ROL

#endif
