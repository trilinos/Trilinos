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


#ifndef ROL_FLETCHERBASE_H
#define ROL_FLETCHERBASE_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_KrylovFactory.hpp"
#include "ROL_PartitionedVector.hpp"
#include <iostream>

namespace ROL {

template <class Real>
class FletcherBase : public Objective<Real> {

protected:
  const Ptr<Objective<Real> > obj_;
  const Ptr<Constraint<Real> > con_;

  Real penaltyParameter_;
  Real quadPenaltyParameter_;

  // Evaluation counters
  int nfval_;
  int ngval_;
  int ncval_;

  Real fPhi_;                   // value of penalty function
  Ptr<Vector<Real> > gPhi_;     // gradient of penalty function

  Ptr<Vector<Real> > y_;        // multiplier estimate

  Real fval_;                   // value of objective function
  Ptr<Vector<Real> > g_;        // gradient of objective value
  Ptr<Vector<Real> > c_;        // constraint value
  Ptr<Vector<Real> > scaledc_;  // penaltyParameter_ * c_
  Ptr<Vector<Real> > gL_;       // gradient of Lagrangian (g - A*y)

  Real cnorm_;                  // norm of constraint violation

  bool isValueComputed_;
  bool isGradientComputed_;
  bool isMultiplierComputed_;
  bool isObjValueComputed_;
  bool isObjGradComputed_;
  bool isConValueComputed_;

  Real multSolverError_;         // Error from augmented system solve in value()
  Real gradSolveError_;          // Error from augmented system solve in gradient()

  Real delta_;                  // regularization parameter

  bool useInexact_;

  // For Augmented system solves
  Ptr<Krylov<Real> > krylov_;
  int iterKrylov_;
  int flagKrylov_;
  Ptr<Vector<Real> > v1_;
  Ptr<Vector<Real> > v2_;
  Ptr<PartitionedVector<Real> > vv_;
  Ptr<Vector<Real> > b1_;
  Ptr<Vector<Real> > b2_;
  Ptr<PartitionedVector<Real> > bb_;
  Ptr<Vector<Real> > w1_;
  Ptr<Vector<Real> > w2_;
  Ptr<PartitionedVector<Real> > ww_;

  void objValue(const Vector<Real>& x, Real &tol) {
    if( !isObjValueComputed_ ) {
      fval_ = obj_->value(x,tol); nfval_++;
      isObjValueComputed_ = true;
    }
  }

  void objGrad(const Vector<Real>& x, Real &tol) {
    if( !isObjGradComputed_ ) {
      obj_->gradient(*g_, x, tol); ngval_++;
      isObjGradComputed_ = true;
    }    
  }

  void conValue(const Vector<Real>&x, Real &tol) {
    if( !isConValueComputed_ ) {
      con_->value(*c_,x,tol); ncval_++;
      scaledc_->set(*c_);
      scaledc_->scale(penaltyParameter_);
      isConValueComputed_ = true;
    }
  }

  virtual void computeMultipliers(const Vector<Real>& x, Real tol) {}

public:
  FletcherBase(const ROL::Ptr<Objective<Real> > &obj,
               const ROL::Ptr<Constraint<Real> > &con)
  : obj_(obj), con_(con), nfval_(0), ngval_(0), ncval_(0), fPhi_(0), 
    isValueComputed_(false), isGradientComputed_(false),
    isMultiplierComputed_(false), isObjValueComputed_(false),
    isObjGradComputed_(false), isConValueComputed_(false),
    multSolverError_(0), gradSolveError_(0),
    iterKrylov_(0), flagKrylov_(0) {}

  // Accessors
  const Ptr<Vector<Real>> getLagrangianGradient(const Vector<Real>& x) {
    // TODO: Figure out reasonable tolerance
    if( !isMultiplierComputed_ ) {
      Real tol = static_cast<Real>(1e-12);
      computeMultipliers(x, tol);
    }
    return gL_;
  }

  const Ptr<Vector<Real>> getConstraintVec(const Vector<Real>& x) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    conValue(x, tol);  
    return c_;
  }

  const Ptr<Vector<Real>> getMultiplierVec(const Vector<Real>& x) {
    // TODO: Figure out reasonable tolerance
    if( !isMultiplierComputed_ ) {
      Real tol = static_cast<Real>(1e-12);
      computeMultipliers(x, tol);
    }
    return y_;
  }

  const Ptr<Vector<Real>> getGradient(const Vector<Real>& x) {
    if( !isGradientComputed_ ) {
      // TODO: Figure out reasonable tolerance
      Real tol = static_cast<Real>(1e-12);
      this->gradient(*gPhi_, x, tol);
    }
    return gPhi_;
  }

  Real getObjectiveValue(const Vector<Real>& x) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    objValue(x, tol);

    return fval_;
  }

  int getNumberFunctionEvaluations() const {
    return nfval_;
  } 

  int getNumberGradientEvaluations() const {
    return ngval_;
  } 

  int getNumberConstraintEvaluations() const {
    return ncval_;
  }

  void setDelta(Real delta) {
    delta_ = delta;
    isValueComputed_ = false;
    isGradientComputed_ = false;
  }

  void setPenaltyParameter( Real sigma ) {
    penaltyParameter_ = sigma;
    isValueComputed_ = false;
    isGradientComputed_ = false;
  }

}; // class Fletcher

} // namespace ROL

#include "ROL_Fletcher.hpp"
#include "ROL_BoundFletcher.hpp"

#endif
