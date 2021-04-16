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


#ifndef ROL_FLETCHEROBJECTVEBASE_H
#define ROL_FLETCHEROBJECTVEBASE_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_Types.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_KrylovFactory.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_ScalarController.hpp"
#include <iostream>

namespace ROL {

template<typename Real>
class FletcherObjectiveBase : public Objective<Real> {
protected:
  const Ptr<Objective<Real>> obj_;
  const Ptr<Constraint<Real>> con_;

  Real sigma_;                // penalty parameter
  Real delta_;                // regularization parameter
  Real quadPenaltyParameter_; // augmented Lagrangian penalty parameter
  bool useInexact_;
  int HessianApprox_;

  // Evaluation counters
  int nfval_, ngval_, ncval_;

  Ptr<ScalarController<Real,int>> fPhi_; // value of penalty function
  Ptr<VectorController<Real,int>> gPhi_; // gradient of penalty function
  Ptr<VectorController<Real,int>> y_;    // multiplier estimate
  Ptr<ScalarController<Real,int>> fval_; // value of objective function
  Ptr<VectorController<Real,int>> g_;    // gradient of objective value
  Ptr<VectorController<Real,int>> c_;    // constraint value

  Ptr<Vector<Real>> scaledc_;  // penaltyParameter_ * c_
  Ptr<Vector<Real>> gL_;       // gradient of Lagrangian (g - A*y)
  Ptr<Vector<Real>> gLdual_;   // dual gradient of Lagrangian (g - A*y)
  Ptr<Vector<Real>> xprim_, xdual_, cprim_, cdual_;

  Real cnorm_;                 // norm of constraint violation

  Real multSolverError_;        // Error from augmented system solve in value()
  Real gradSolveError_;         // Error from augmented system solve in gradient()

  // For Augmented system solves
  Ptr<Krylov<Real>> krylov_;
  int iterKrylov_, flagKrylov_;
  Ptr<Vector<Real>> v1_, v2_, b1_, b2_, w1_, w2_;
  Ptr<PartitionedVector<Real>> vv_, bb_, ww_;

public:
  FletcherObjectiveBase(const Ptr<Objective<Real>> &obj,
                        const Ptr<Constraint<Real>> &con,
                        const Vector<Real> &xprim,
                        const Vector<Real> &xdual,
                        const Vector<Real> &cprim,
                        const Vector<Real> &cdual,
                        ParameterList &parlist);

  virtual void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;

  // Accessors
  Ptr<const Vector<Real>> getLagrangianGradient(const Vector<Real>& x);
  Ptr<const Vector<Real>> getConstraintVec(const Vector<Real>& x);
  Ptr<const Vector<Real>> getMultiplierVec(const Vector<Real>& x);
  Ptr<const Vector<Real>> getGradient(const Vector<Real>& x);
  Real getObjectiveValue(const Vector<Real>& x);
  int getNumberFunctionEvaluations() const;
  int getNumberGradientEvaluations() const;
  int getNumberConstraintEvaluations() const;
  void reset(Real sigma, Real delta);

protected:
  Real objValue(const Vector<Real>& x, Real &tol);
  void objGrad(Vector<Real> &g, const Vector<Real>& x, Real &tol);
  void conValue(Vector<Real> &c, const Vector<Real>&x, Real &tol);
  void computeMultipliers(Vector<Real> &y, Vector<Real> &gL, const Vector<Real> &x, Vector<Real> &g, Vector<Real> &c, Real tol);
  virtual void solveAugmentedSystem(Vector<Real> &v1, Vector<Real> &v2, const Vector<Real> &b1, const Vector<Real> &b2, const Vector<Real> &x, Real &multSolverError_, bool refine) = 0;

}; // class FletcherObjectiveBase

} // namespace ROL

#include "ROL_FletcherObjectiveBase_Def.hpp"

#endif
