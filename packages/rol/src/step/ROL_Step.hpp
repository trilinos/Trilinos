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

#ifndef ROL_STEP_H
#define ROL_STEP_H

#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_ParameterList.hpp"

/** @ingroup step_group
    \class ROL::Step
    \brief Provides the interface to compute optimization steps.
*/


namespace ROL {

template <class Real>
class Step {
private:
  Teuchos::RCP<StepState<Real> > state_;

protected:
  Teuchos::RCP<StepState<Real> > getState(void) {
    return this->state_;
  }

public:

  virtual ~Step() {}

  Step(void) { 
    state_ = Teuchos::rcp( new StepState<Real> );
  }


  /** \brief Initialize step with bound constraint.
  */
  virtual void initialize( Vector<Real> &x, Objective<Real> &obj, BoundConstraint<Real> &con, 
                           AlgorithmState<Real> &algo_state ) {
    Real tol = std::sqrt(ROL_EPSILON);
    // Initialize state descent direction and gradient storage
    state_->descentVec   = x.clone();
    state_->gradientVec  = x.clone();
    state_->searchSize = 0.0;
    // Project x onto constraint set
    if ( con.isActivated() ) {
      con.project(x);
    }
    // Update objective function, get value, and get gradient
    obj.update(x,true,algo_state.iter);
    algo_state.value = obj.value(x,tol);
    algo_state.nfval++;
    obj.gradient(*(state_->gradientVec),x,tol);
    algo_state.ngrad++;
    if ( con.isActivated() ) {
      Teuchos::RCP<Vector<Real> > xnew = x.clone();
      xnew->set(x);
      xnew->axpy(-1.0,*(Step<Real>::state_->gradientVec));
      con.project(*xnew);
      xnew->axpy(-1.0,x);
      algo_state.gnorm = xnew->norm();
    }
    else {
      algo_state.gnorm = (state_->gradientVec)->norm();
    }
  }

  /** \brief Initialize step with equality constraint.
  */
  virtual void initialize( Vector<Real> &x, Vector<Real> &g, Vector<Real> &l, Vector<Real> &c,
                           Objective<Real> &obj, EqualityConstraint<Real> &con, 
                           AlgorithmState<Real> &algo_state ) {
  }

  /** \brief Compute step.
  */
  virtual void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, 
                        BoundConstraint<Real> &con, 
                        AlgorithmState<Real> &algo_state ) = 0;

  /** \brief Update step, if successful.
  */
  virtual void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, 
                       BoundConstraint<Real> &con,
                       AlgorithmState<Real> &algo_state ) = 0;

  /** \brief Compute step (equality constraints).
  */
  virtual void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                        Objective<Real> &obj, EqualityConstraint<Real> &con, 
                        AlgorithmState<Real> &algo_state ) {}

  /** \brief Update step, if successful (equality constraints).
  */
  virtual void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
                       Objective<Real> &obj, EqualityConstraint<Real> &con,
                       AlgorithmState<Real> &algo_state ) {}

  /** \brief Print iterate header.
  */
  virtual std::string printHeader( void ) const = 0;

  /** \brief Print step name.
  */
  virtual std::string printName( void ) const = 0;

  /** \brief Print iterate status.
  */
  virtual std::string print( AlgorithmState<Real> &algo_state, bool printHeader = false ) const = 0;

  /** \brief Get state for step object.
  */
  virtual Teuchos::RCP<const StepState<Real> > getStepState(void) const {
    return this->state_;
  }

  // struct StepState (scalars, vectors) map?

  // getState

  // setState

}; // class Step

} // namespace ROL

#endif
