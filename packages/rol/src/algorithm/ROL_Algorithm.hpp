// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ALGORITHM_H
#define ROL_ALGORITHM_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_ValidParameters.hpp"

/** \class ROL::Algorithm
    \brief Provides an interface to run optimization algorithms.
*/


namespace ROL {

template <class Real>
class Algorithm {
private:
  ROL::Ptr<Step<Real> >           step_;
  ROL::Ptr<StatusTest<Real> >     status_;
  ROL::Ptr<AlgorithmState<Real> > state_;

  bool printHeader_;

public:

  virtual ~Algorithm() {}

  /** \brief Constructor, given a step and a status test.
  */
  Algorithm( const ROL::Ptr<Step<Real> > & step,
             const ROL::Ptr<StatusTest<Real> > & status,
             bool printHeader = false ) {
    step_ = step;
    status_ = status;
    state_ = ROL::makePtr<AlgorithmState<Real>>();
    printHeader_ = printHeader;
  }

  /** \brief Constructor, given a step, a status test, and a
             previously defined algorithm state.
  */
  Algorithm( const ROL::Ptr<Step<Real> > & step,
             const ROL::Ptr<StatusTest<Real> > & status,
             const ROL::Ptr<AlgorithmState<Real> > & state,
             bool printHeader = false ) {
    step_ = step;
    status_ = status;
    state_ = state;
    printHeader_ = printHeader;
  }

  /** \brief Run algorithm on unconstrained problems (Type-U).
             This is the primary Type-U interface.
  */
  virtual std::vector<std::string> run( Vector<Real>      &x,
                                        Objective<Real>   &obj,
                                        bool              print = false,
                                        std::ostream      &outStream = std::cout,
                                        bool              printVectors = false,
                                        std::ostream      &vectorStream = std::cout ) {
    BoundConstraint<Real> bnd;
    bnd.deactivate();
    return run(x,x.dual(),obj,bnd,print,outStream,printVectors,vectorStream);
  }

  /** \brief Run algorithm on unconstrained problems (Type-U).
             This general interface supports the use of dual optimization vector spaces,
             where the user does not define the dual() method.
  */
  virtual std::vector<std::string> run( Vector<Real>       &x,
                                        const Vector<Real> &g, 
                                        Objective<Real>    &obj,
                                        bool               print = false,
                                        std::ostream       &outStream = std::cout,
                                        bool               printVectors = false,
                                        std::ostream       &vectorStream = std::cout ) {
    BoundConstraint<Real> bnd;
    bnd.deactivate();
    return run(x,g,obj,bnd,print,outStream,printVectors,vectorStream);
  }

  /** \brief Run algorithm on bound constrained problems (Type-B).
             This is the primary Type-B interface.
  */
  virtual std::vector<std::string> run( Vector<Real>          &x, 
                                        Objective<Real>       &obj,
                                        BoundConstraint<Real> &bnd,
                                        bool                  print = false,
                                        std::ostream          &outStream = std::cout,
                                        bool                  printVectors = false,
                                        std::ostream          &vectorStream = std::cout ) {
    return run(x,x.dual(),obj,bnd,print,outStream,printVectors,vectorStream);
  }

  /** \brief Run algorithm on bound constrained problems (Type-B).
             This general interface supports the use of dual optimization vector spaces,
             where the user does not define the dual() method.
  */
  virtual std::vector<std::string> run( Vector<Real>          &x, 
                                        const Vector<Real>    &g, 
                                        Objective<Real>       &obj,
                                        BoundConstraint<Real> &bnd,
                                        bool                  print = false,
                                        std::ostream          &outStream = std::cout,
                                        bool                  printVectors = false,
                                        std::ostream          &vectorStream = std::cout ) {
    if(printVectors) {
      x.print(vectorStream);
    }

    std::vector<std::string> output;

    // Initialize Current Iterate Container 
    if ( state_->iterateVec == ROL::nullPtr ) {
      state_->iterateVec = x.clone();
    }
    state_->iterateVec->set(x);

    // Initialize Step Container
    ROL::Ptr<Vector<Real> > s = x.clone();

    // Initialize Step
    step_->initialize(x, g, obj, bnd, *state_);
    output.push_back(step_->print(*state_,true));
    if ( print ) {
      outStream << step_->print(*state_,true);
    }

    // Initialize Minimum Value and Vector
    if ( state_->minIterVec == ROL::nullPtr ) {
      state_->minIterVec = x.clone();
    }
    state_->minIterVec->set(x);
    state_->minIter = state_->iter;
    state_->minValue = state_->value;

    // Run Algorithm
    while (status_->check(*state_)) {
      step_->compute(*s, x, obj, bnd, *state_);
      step_->update(x, *s, obj, bnd, *state_);

      if( printVectors ) {
        x.print(vectorStream);
      }

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
    std::stringstream hist;
    hist << "Optimization Terminated with Status: ";
    hist << EExitStatusToString(state_->statusFlag);
    hist << "\n";
    output.push_back(hist.str());
    if ( print ) {
      outStream << hist.str();
    }
    return output;
  }


  /** \brief Run algorithm on equality constrained problems (Type-E).
             This is the primary Type-E interface.
  */
  virtual std::vector<std::string> run( Vector<Real>             &x,
                                        Vector<Real>             &l, 
                                        Objective<Real>          &obj,
                                        Constraint<Real>         &con,
                                        bool                     print = false,
                                        std::ostream             &outStream = std::cout,
                                        bool                     printVectors = false,
                                        std::ostream             &vectorStream = std::cout ) {

    return run(x, x.dual(), l, l.dual(), obj, con, print, outStream, printVectors, vectorStream);

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
                                        Constraint<Real>         &con,
                                        bool                     print = false,
                                        std::ostream             &outStream = std::cout,
                                        bool                     printVectors = false,
                                        std::ostream             &vectorStream = std::cout ) {
    if( printVectors ) {
      x.print(vectorStream);
    } 

    std::vector<std::string> output;

    // Initialize Current Iterate Container 
    if ( state_->iterateVec == ROL::nullPtr ) {
      state_->iterateVec = x.clone();
    }
    state_->iterateVec->set(x);

    // Initialize Current Lagrange Multiplier Container 
    if ( state_->lagmultVec == ROL::nullPtr ) {
      state_->lagmultVec = l.clone();
    }
    state_->lagmultVec->set(l);

    // Initialize Step Container
    ROL::Ptr<Vector<Real> > s = x.clone();

    // Initialize Step
    step_->initialize(x, g, l, c, obj, con, *state_);
    output.push_back(step_->print(*state_,true));
    if ( print ) {
      outStream << step_->print(*state_,true);
    }

    // Initialize Minimum Value and Vector
    if ( state_->minIterVec == ROL::nullPtr ) {
      state_->minIterVec = x.clone();
    }
    state_->minIterVec->set(x);
    state_->minIter = state_->iter;
    state_->minValue = state_->value;

    // Run Algorithm
    while (status_->check(*state_)) {
      step_->compute(*s, x, l, obj, con, *state_);
      step_->update(x, l, *s, obj, con, *state_);

      if( printVectors ) { 
        x.print(vectorStream);
      } 

      output.push_back(step_->print(*state_,printHeader_));
      if ( print ) {
        outStream << step_->print(*state_,printHeader_);
      }
    }
    std::stringstream hist;
    hist << "Optimization Terminated with Status: ";
    hist << EExitStatusToString(state_->statusFlag);
    hist << "\n";
    output.push_back(hist.str());
    if ( print ) {
      outStream << hist.str();
    }
    return output;
  }

  /** \brief Run algorithm on equality and bound constrained problems (Type-EB).
             This is the primary Type-EB interface.
  */
  virtual std::vector<std::string> run( Vector<Real>             &x,
                                        Vector<Real>             &l, 
                                        Objective<Real>          &obj,
                                        Constraint<Real>         &con,
                                        BoundConstraint<Real>    &bnd,
                                        bool                     print = false,
                                        std::ostream             &outStream = std::cout,
                                        bool                     printVectors = false,
                                        std::ostream             &vectorStream = std::cout) {
    return run(x,x.dual(),l,l.dual(),obj,con,bnd,print,outStream,printVectors,vectorStream);
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
                                        Constraint<Real>         &con,
                                        BoundConstraint<Real>    &bnd,
                                        bool                     print = false,
                                        std::ostream             &outStream = std::cout,
                                        bool                     printVectors = false,
                                        std::ostream             &vectorStream = std::cout ) {
    if(printVectors) {
      x.print(vectorStream); 
    } 

    std::vector<std::string> output;

    // Initialize Current Iterate Container 
    if ( state_->iterateVec == ROL::nullPtr ) {
      state_->iterateVec = x.clone();
    }
    state_->iterateVec->set(x);

    // Initialize Current Lagrange Multiplier Container 
    if ( state_->lagmultVec == ROL::nullPtr ) {
      state_->lagmultVec = l.clone();
    }
    state_->lagmultVec->set(l);

    // Initialize Step Container
    ROL::Ptr<Vector<Real> > s = x.clone();

    // Initialize Step
    step_->initialize(x, g, l, c, obj, con, bnd, *state_);
    output.push_back(step_->print(*state_,true));
    if ( print ) {
      outStream << step_->print(*state_,true);
    }

    // Initialize Minimum Value and Vector
    if ( state_->minIterVec == ROL::nullPtr ) {
      state_->minIterVec = x.clone();
    }
    state_->minIterVec->set(x);
    state_->minIter = state_->iter;
    state_->minValue = state_->value;

    // Run Algorithm
    while (status_->check(*state_)) {
      step_->compute(*s, x, l, obj, con, bnd, *state_);
      step_->update(x, l, *s, obj, con, bnd, *state_);
      if( printVectors ) {
        x.print(vectorStream);
      }
      output.push_back(step_->print(*state_,printHeader_));
      if ( print ) {
        outStream << step_->print(*state_,printHeader_);
      }
    }
    std::stringstream hist;
    hist << "Optimization Terminated with Status: ";
    hist << EExitStatusToString(state_->statusFlag);
    hist << "\n";
    output.push_back(hist.str());
    if ( print ) {
      outStream << hist.str();
    }
    return output;
  }

  std::string getIterHeader(void) {
    return step_->printHeader();
  }

  std::string getIterInfo(bool withHeader = false) {
    return step_->print(*state_,withHeader);
  }

  ROL::Ptr<const AlgorithmState<Real> > getState(void) const {
    return state_;
  }

  void reset(void) {
    state_->reset();
  }






}; // class Algorithm


} // namespace ROL

#endif
