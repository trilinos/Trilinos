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

#ifndef ROL_ALGORITHM_H
#define ROL_ALGORITHM_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_StepFactory.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_StatusTestFactory.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_EqualityConstraint.hpp"

#include "ROL_OptimizationProblem.hpp"

/** \class ROL::Algorithm
    \brief Provides an interface to run optimization algorithms.
*/


namespace ROL {

template<class Real>
class StepFactory;

template<class Real>
class StatusTestFactory;

template <class Real>
class Algorithm {
private:
  Teuchos::RCP<Step<Real> >           step_;
  Teuchos::RCP<StatusTest<Real> >     status_;
  Teuchos::RCP<AlgorithmState<Real> > state_;

  bool printHeader_;

public:

  virtual ~Algorithm() {}

  /** \brief Constructor, given a step and a status test.
  */
  Algorithm( const Teuchos::RCP<Step<Real> > & step,
             const Teuchos::RCP<StatusTest<Real> > & status,
             bool printHeader = false ) {
    step_ = step;
    status_ = status;
    state_ = Teuchos::rcp(new AlgorithmState<Real>);
    printHeader_ = printHeader;
  }

  /** \brief Constructor, given a step, a status test, and a
             previously defined algorithm state.
  */
  Algorithm( const Teuchos::RCP<Step<Real> > & step,
             const Teuchos::RCP<StatusTest<Real> > & status,
             const Teuchos::RCP<AlgorithmState<Real> > & state,
             bool printHeader = false ) {
    step_ = step;
    status_ = status;
    state_ = state;
    printHeader_ = printHeader;
  }

  /** \brief Constructor, given a string, for the step, and a
             parameter list of various options.  The status
             test is determined based on the step string.
  */
  Algorithm( const std::string &stepname,
             Teuchos::ParameterList &parlist,
             bool printHeader = false) {
    EStep els = StringToEStep(stepname);
    TEUCHOS_TEST_FOR_EXCEPTION( !(isValidStep(els)),
                                std::invalid_argument,
                                "Invalid step name in algorithm constructor!");
    StepFactory<Real> stepFactory;
    StatusTestFactory<Real> statusTestFactory;
    step_   = stepFactory.getStep(stepname,parlist);
    status_ = statusTestFactory.getStatusTest(stepname,parlist);
    state_  = Teuchos::rcp(new AlgorithmState<Real>);
    printHeader_ = printHeader;
  }

  /** \brief Run algorithm on unconstrained problems (Type-U).
             This is the primary Type-U interface.
  */
  virtual std::vector<std::string> run( Vector<Real>      &x,
                                        Objective<Real>   &obj,
                                        bool              print = false,
                                        std::ostream      &outStream = std::cout ) {
    BoundConstraint<Real> con;
    con.deactivate();
    return run(x,x.dual(),obj,con,print,outStream);
  }

  /** \brief Run algorithm on unconstrained problems (Type-U).
             This general interface supports the use of dual optimization vector spaces,
             where the user does not define the dual() method.
  */
  virtual std::vector<std::string> run( Vector<Real>       &x,
                                        const Vector<Real> &g, 
                                        Objective<Real>    &obj,
                                        bool               print = false,
                                        std::ostream       &outStream = std::cout ) {
    BoundConstraint<Real> con;
    con.deactivate();
    return run(x,g,obj,con,print,outStream);
  }

  /** \brief Run algorithm on bound constrained problems (Type-B).
             This is the primary Type-B interface.
  */
  virtual std::vector<std::string> run( Vector<Real>          &x, 
                                        Objective<Real>       &obj,
                                        BoundConstraint<Real> &con,
                                        bool                  print = false,
                                        std::ostream          &outStream = std::cout ) {
    return run(x,x.dual(),obj,con,print,outStream);
  }

  /** \brief Run algorithm on bound constrained problems (Type-B).
             This general interface supports the use of dual optimization vector spaces,
             where the user does not define the dual() method.
  */
  virtual std::vector<std::string> run( Vector<Real>          &x, 
                                        const Vector<Real>    &g, 
                                        Objective<Real>       &obj,
                                        BoundConstraint<Real> &con,
                                        bool                  print = false,
                                        std::ostream          &outStream = std::cout ) {
    std::vector<std::string> output;

    // Initialize Current Iterate Container 
    if ( state_->iterateVec == Teuchos::null ) {
      state_->iterateVec = x.clone();
    }
    state_->iterateVec->set(x);

    // Initialize Step Container
    Teuchos::RCP<Vector<Real> > s = x.clone();

    // Initialize Step
    step_->initialize(x, g, obj, con, *state_);
    output.push_back(step_->print(*state_,true));
    if ( print ) {
      outStream << step_->print(*state_,true);
    }

    // Initialize Minimum Value and Vector
    if ( state_->minIterVec == Teuchos::null ) {
      state_->minIterVec = x.clone();
    }
    state_->minIterVec->set(x);
    state_->minIter = state_->iter;
    state_->minValue = state_->value;

    // Run Algorithm
    while (status_->check(*state_)) {
      step_->compute(*s, x, obj, con, *state_);
      step_->update(x, *s, obj, con, *state_);
      // Store Minimal Value and Vector
      if ( state_->minValue > state_->value ) {
        state_->minIterVec->set(*(state_->iterateVec));
        state_->minValue = state_->value;
        state_->minIter = state_->iter;
      }
      // Update Output
      output.push_back(step_->print(*state_,printHeader_));
      if ( print ) {
        outStream << step_->print(*state_,printHeader_);
      }
    }
    return output;
  }


  /** \brief Run algorithm on equality constrained problems (Type-E).
             This is the primary Type-E interface.
  */
  virtual std::vector<std::string> run( Vector<Real>             &x,
                                        Vector<Real>             &l, 
                                        Objective<Real>          &obj,
                                        EqualityConstraint<Real> &con,
                                        bool                     print = false,
                                        std::ostream             &outStream = std::cout ) {

    return run(x, x.dual(), l, l.dual(), obj, con, print, outStream);

  }


  /** \brief Run algorithm on equality constrained problems (Type-E).
             This general interface supports the use of dual optimization and
             constraint vector spaces, where the user does not define the dual() method.
  */
  virtual std::vector<std::string> run( Vector<Real>             &x,
                                        const Vector<Real>       &g, 
                                        Vector<Real>             &l, 
                                        const Vector<Real>       &c, 
                                        Objective<Real>          &obj,
                                        EqualityConstraint<Real> &con,
                                        bool                     print = false,
                                        std::ostream             &outStream = std::cout ) {
    std::vector<std::string> output;

    // Initialize Current Iterate Container 
    if ( state_->iterateVec == Teuchos::null ) {
      state_->iterateVec = x.clone();
    }
    state_->iterateVec->set(x);

    // Initialize Current Lagrange Multiplier Container 
    if ( state_->lagmultVec == Teuchos::null ) {
      state_->lagmultVec = l.clone();
    }
    state_->lagmultVec->set(l);

    // Initialize Step Container
    Teuchos::RCP<Vector<Real> > s = x.clone();

    // Initialize Step
    step_->initialize(x, g, l, c, obj, con, *state_);
    output.push_back(step_->print(*state_,true));
    if ( print ) {
      outStream << step_->print(*state_,true);
    }

    // Initialize Minimum Value and Vector
    if ( state_->minIterVec == Teuchos::null ) {
      state_->minIterVec = x.clone();
    }
    state_->minIterVec->set(x);
    state_->minIter = state_->iter;
    state_->minValue = state_->value;

    // Run Algorithm
    while (status_->check(*state_)) {
      step_->compute(*s, x, l, obj, con, *state_);
      step_->update(x, l, *s, obj, con, *state_);
      output.push_back(step_->print(*state_,printHeader_));
      if ( print ) {
        outStream << step_->print(*state_,printHeader_);
      }
    }
    return output;
  }

  /** \brief Run algorithm on equality and bound constrained problems (Type-EB).
             This is the primary Type-EB interface.
  */
  virtual std::vector<std::string> run( Vector<Real>             &x,
                                        Vector<Real>             &l, 
                                        Objective<Real>          &obj,
                                        EqualityConstraint<Real> &con,
                                        BoundConstraint<Real>    &bnd,
                                        bool                     print = false,
                                        std::ostream             &outStream = std::cout ) {
    return run(x,x.dual(),l,l.dual(),obj,con,bnd,print,outStream);
  }

  /** \brief Run algorithm on equality and bound constrained problems (Type-EB).
             This general interface supports the use of dual optimization and
             constraint vector spaces, where the user does not define the dual() method.
  */
  virtual std::vector<std::string> run( Vector<Real>             &x,
                                        const Vector<Real>       &g, 
                                        Vector<Real>             &l, 
                                        const Vector<Real>       &c, 
                                        Objective<Real>          &obj,
                                        EqualityConstraint<Real> &con,
                                        BoundConstraint<Real>    &bnd,
                                        bool                     print = false,
                                        std::ostream             &outStream = std::cout ) {
    std::vector<std::string> output;

    // Initialize Current Iterate Container 
    if ( state_->iterateVec == Teuchos::null ) {
      state_->iterateVec = x.clone();
    }
    state_->iterateVec->set(x);

    // Initialize Current Lagrange Multiplier Container 
    if ( state_->lagmultVec == Teuchos::null ) {
      state_->lagmultVec = l.clone();
    }
    state_->lagmultVec->set(l);

    // Initialize Step Container
    Teuchos::RCP<Vector<Real> > s = x.clone();

    // Initialize Step
    step_->initialize(x, g, l, c, obj, con, bnd, *state_);
    output.push_back(step_->print(*state_,true));
    if ( print ) {
      outStream << step_->print(*state_,true);
    }

    // Initialize Minimum Value and Vector
    if ( state_->minIterVec == Teuchos::null ) {
      state_->minIterVec = x.clone();
    }
    state_->minIterVec->set(x);
    state_->minIter = state_->iter;
    state_->minValue = state_->value;

    // Run Algorithm
    while (status_->check(*state_)) {
      step_->compute(*s, x, l, obj, con, bnd, *state_);
      step_->update(x, l, *s, obj, con, bnd, *state_);
      output.push_back(step_->print(*state_,printHeader_));
      if ( print ) {
        outStream << step_->print(*state_,printHeader_);
      }
    }
    return output;
  }

  /** \brief Run algorithm using a ROL::OptimizationProblem.
  */
  virtual std::vector<std::string> run( OptimizationProblem<Real> &opt,
                                        bool                     print = false,
                                        std::ostream             &outStream = std::cout ) {
    // Get components of optimization problem
    Teuchos::RCP<Objective<Real> >          obj = opt.getObjective();
    Teuchos::RCP<Vector<Real> >             x   = opt.getSolutionVector();
    Teuchos::RCP<BoundConstraint<Real> >    bnd = opt.getBoundConstraint();
    Teuchos::RCP<EqualityConstraint<Real> > con = opt.getEqualityConstraint();
    Teuchos::RCP<Vector<Real> >             l   = opt.getMultiplierVector();
    // Call appropriate run function
    if ( con == Teuchos::null ) {
      if ( bnd == Teuchos::null ) {
        return run(*x,*obj,print,outStream);
      }
      else {
        return run(*x,*obj,*bnd,print,outStream);
      }
    }
    else {
      if ( bnd == Teuchos::null ) {
        return run(*x,*l,*obj,*con,print,outStream);
      }
      else {
        return run(*x,*l,*obj,*con,*bnd,print,outStream);
      }
    }
  }

  std::string getIterHeader(void) {
    return step_->printHeader();
  }

  std::string getIterInfo(bool withHeader = false) {
    return step_->print(*state_,withHeader);
  }

  Teuchos::RCP<const AlgorithmState<Real> > getState(void) const {
    return state_;
  }

  void reset(void) {
    state_  = Teuchos::rcp(new AlgorithmState<Real>);
  }

}; // class Algorithm


} // namespace ROL

#endif
