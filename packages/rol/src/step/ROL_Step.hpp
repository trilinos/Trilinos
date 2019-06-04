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
#include "ROL_Constraint.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_Types.hpp"

#include "ROL_ParameterList.hpp"

/** @ingroup step_group
    \class ROL::Step
    \brief Provides the interface to compute optimization steps.
*/


namespace ROL {

// We need a forward declaration here, because some steps are algorithms.
template<class Real>
class Algorithm;

template <class Real>
class Step {
private:
  ROL::Ptr<StepState<Real> > state_;

protected:
  ROL::Ptr<StepState<Real> > getState(void) {
    return state_;
  }

public:

  virtual ~Step() {}

  Step(void) { 
    state_ = ROL::makePtr<StepState<Real>>();
  }


  /** \brief Initialize step with bound constraint.
  */
  virtual void initialize( Vector<Real> &x, const Vector<Real> &g,
                           Objective<Real> &obj, BoundConstraint<Real> &con,
                           AlgorithmState<Real> &algo_state ) {
    initialize(x,x,g,obj,con,algo_state);
  }

  /** \brief Initialize step with bound constraint.
  */
  virtual void initialize( Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g,
                           Objective<Real> &obj, BoundConstraint<Real> &con,
                           AlgorithmState<Real> &algo_state ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>()), one(1), zero(0);
    // Initialize state descent direction and gradient storage
    state_->descentVec  = s.clone();
    state_->gradientVec = g.clone();
    state_->searchSize  = zero;
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
      ROL::Ptr<Vector<Real> > xnew = x.clone();
      xnew->set(x);
      xnew->axpy(-one,(Step<Real>::state_->gradientVec)->dual());
      con.project(*xnew);
      xnew->axpy(-one,x);
      algo_state.gnorm = xnew->norm();
    }
    else {
      algo_state.gnorm = (state_->gradientVec)->norm();
    }
  }

  /** \brief Initialize step with equality constraint.
  */
  virtual void initialize( Vector<Real> &x, const Vector<Real> &g, Vector<Real> &l, const Vector<Real> &c,
                           Objective<Real> &obj, Constraint<Real> &con,
                           AlgorithmState<Real> &algo_state ) {
  }

  /** \brief Initialize step with equality constraint.
  */
  virtual void initialize( Vector<Real> &x, const Vector<Real> &g, Vector<Real> &l, const Vector<Real> &c,
                           Objective<Real> &obj, Constraint<Real> &con, BoundConstraint<Real> &bnd,
                           AlgorithmState<Real> &algo_state ) {
  }

  /** \brief Compute step.
  */
  virtual void compute( Vector<Real> &s, const Vector<Real> &x,
                        Objective<Real> &obj, BoundConstraint<Real> &bnd,
                        AlgorithmState<Real> &algo_state ) {
    throw Exception::NotImplemented(">>> ROL::Step::compute(s,x,obj,bnd,algo_state) is not implemented!");
  }

  /** \brief Update step, if successful.
  */
  virtual void update( Vector<Real> &x, const Vector<Real> &s,
                       Objective<Real> &obj, BoundConstraint<Real> &bnd,
                       AlgorithmState<Real> &algo_state ) {
    throw Exception::NotImplemented(">>> ROL::Step::update(x,s,obj,bnd,algo_state) is not implemented!");
  }

  /** \brief Compute step (equality constraints).
  */
  virtual void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                        Objective<Real> &obj, Constraint<Real> &con,
                        AlgorithmState<Real> &algo_state ) {
    throw Exception::NotImplemented(">>> ROL::Step::compute(s,x,l,obj,con,algo_state) is not implemented!");
  }

  /** \brief Update step, if successful (equality constraints).
  */
  virtual void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
                       Objective<Real> &obj, Constraint<Real> &con,
                       AlgorithmState<Real> &algo_state ) {
    throw Exception::NotImplemented(">>> ROL::Step::update(x,s,l,obj,bnd,con,algo_state) is not implemented!");
  }

  /** \brief Compute step (equality constraints).
  */
  virtual void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                        Objective<Real> &obj, Constraint<Real> &con,
                        BoundConstraint<Real> &bnd,
                        AlgorithmState<Real> &algo_state ) {
    throw Exception::NotImplemented(">>> ROL::Step::compute(s,x,l,obj,bnd,con,algo_state) is not implemented!");
  }

  /** \brief Update step, if successful (equality constraints).
  */
  virtual void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
                       Objective<Real> &obj, Constraint<Real> &con,
                       BoundConstraint<Real> &bnd,
                       AlgorithmState<Real> &algo_state ) {
    throw Exception::NotImplemented(">>> ROL::Step::update(x,s,l,obj,bnd,con,algo_state) is not implemented!");
  }



  // Methods using an Optimization problem

  void initialize( OptimizationProblem<Real> &opt, AlgorithmState<Real> &algo_state ) {

    

    ROL::Ptr<Objective<Real> >          obj = opt.getObjective();
    ROL::Ptr<Vector<Real> >             x   = opt.getSolutionVector();
    ROL::Ptr<BoundConstraint<Real> >    bnd = opt.getBoundConstraint();
    ROL::Ptr<Constraint<Real> >         con = opt.getConstraint();
    ROL::Ptr<Vector<Real> >             l   = opt.getMultiplierVector();

    if( con == ROL::nullPtr ) { // has no equality constraint
      if( bnd == ROL::nullPtr ) { // has no bound constraint or inactive
        bnd = ROL::makePtr<BoundConstraint<Real>>();
        bnd->deactivate();
      }
      initialize(*x, x->dual(), *obj, *bnd, algo_state);
    }
    else { // has equality constraint 

      if( bnd == ROL::nullPtr ) {
        initialize(*x,x->dual(),*l,l->dual(),*obj,*con,algo_state );
      }
      initialize(*x,x->dual(),*l,l->dual(),*obj,*con,*bnd,algo_state);
    }
  }

  void compute( Vector<Real> &s, OptimizationProblem<Real> &opt, AlgorithmState<Real> &algo_state ) {
    

    ROL::Ptr<Objective<Real> >          obj = opt.getObjective();
    ROL::Ptr<Vector<Real> >             x   = opt.getSolutionVector();
    ROL::Ptr<BoundConstraint<Real> >    bnd = opt.getBoundConstraint();
    ROL::Ptr<Constraint<Real> >         con = opt.getConstraint();
    ROL::Ptr<Vector<Real> >             l   = opt.getMultiplierVector();

    if( con == ROL::nullPtr ) { // has no equality constraint
      if( bnd == ROL::nullPtr ) { // has no bound constraint
        bnd = ROL::makePtr<BoundConstraint<Real>>();
        bnd->deactivate();
      }
      compute(s,*x, *obj, *bnd, algo_state);
    }
    else { // has equality constraint 
      if( bnd == ROL::nullPtr ) {
        compute(s,*x,*l,*obj,*con,algo_state);
      }
      compute(s,*x,*l,*obj,*con,*bnd,algo_state);
    }
 
  }

  void update( OptimizationProblem<Real> &opt, const Vector<Real> &s, AlgorithmState<Real> &algo_state ) {
    

    ROL::Ptr<Objective<Real> >          obj = opt.getObjective();
    ROL::Ptr<Vector<Real> >             x   = opt.getSolutionVector();
    ROL::Ptr<BoundConstraint<Real> >    bnd = opt.getBoundConstraint();
    ROL::Ptr<Constraint<Real> >         con = opt.getConstraint();
    ROL::Ptr<Vector<Real> >             l   = opt.getMultiplierVector();

    if( con == ROL::nullPtr ) { // has no equality constraint
      if( bnd == ROL::nullPtr ) { // has no bound constraint
        bnd = ROL::makePtr<BoundConstraint<Real>>();
        bnd->deactivate();
      }
      update(*x, s, *obj, *bnd, algo_state);
    }
    else { // has equality constraint
      if( bnd == ROL::nullPtr ) {
        update(*x,*l,s,*obj,*con,algo_state);
      }
      update(*x,*l,s,*obj,*con,*bnd,algo_state);
    }
 
  }


  /** \brief Print iterate header.
  */
  virtual std::string printHeader( void ) const {
    throw Exception::NotImplemented(">>> ROL::Step::printHeader() is not implemented!");
  }

  /** \brief Print step name.
  */
  virtual std::string printName( void ) const {
    throw Exception::NotImplemented(">>> ROL::Step::printName() is not implemented!");
  }

  /** \brief Print iterate status.
  */
  virtual std::string print( AlgorithmState<Real> &algo_state, bool printHeader = false ) const {
    throw Exception::NotImplemented(">>> ROL::Step::print(algo_state,printHeader) is not implemented!");
  }

  /** \brief Get state for step object.
  */
  const ROL::Ptr<const StepState<Real> > getStepState(void) const {
    return state_;
  }

  /** \brief Get state for step object.
  */
  void reset(const Real searchSize = 1.0) {
    state_->reset(searchSize);
  }

  // struct StepState (scalars, vectors) map?

  // getState

  // setState

}; // class Step

} // namespace ROL

#endif
