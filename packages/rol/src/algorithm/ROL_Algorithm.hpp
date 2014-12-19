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
#include "ROL_StatusTest.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_EqualityConstraint.hpp"

/** \class ROL::Algorithm
    \brief Provides an interface to run optimization algorithms.
*/


namespace ROL {

template <class Real>
class DefaultAlgorithm {
private:
  Teuchos::RCP<Step<Real> >           step_;
  Teuchos::RCP<StatusTest<Real> >     status_;
  Teuchos::RCP<AlgorithmState<Real> > state_;

  bool printHeader_;

public:

  virtual ~DefaultAlgorithm() {}

  DefaultAlgorithm(Step<Real> & step, StatusTest<Real> & status, bool printHeader = false ) {
    step_   = Teuchos::rcp(&step,   false);
    status_ = Teuchos::rcp(&status, false);
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
      state_->iterateVec->set(x);
    }

    // Initialize Step Container
    Teuchos::RCP<Vector<Real> > s = x.clone();

    // Initialize Step
    step_->initialize(x, g, obj, con, *state_);
    output.push_back(step_->print(*state_,true));
    if ( print ) {
      outStream << step_->print(*state_,true);
    }

    // Run Algorithm
    while (status_->check(*state_)) {
      step_->compute(*s, x, obj, con, *state_);
      step_->update(x, *s, obj, con, *state_);
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
      state_->iterateVec->set(x);
    }

    // Initialize Current Lagrange Multiplier Container 
    if ( state_->lagmultVec == Teuchos::null ) {
      state_->lagmultVec = l.clone();
      state_->lagmultVec->set(l);
    }

    // Initialize Step Container
    Teuchos::RCP<Vector<Real> > s = x.clone();

    // Initialize Step
    step_->initialize(x, g, l, c, obj, con, *state_);
    output.push_back(step_->print(*state_,true));
    if ( print ) {
      outStream << step_->print(*state_,true);
    }

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

  std::string getIterHeader(void) {
    return step_->printHeader();
  }

  std::string getIterInfo(bool withHeader = false) {
    return step_->print(*state_,withHeader);
  }

  Teuchos::RCP<const AlgorithmState<Real> > getState(void) const {
    return state_;
  }
  

}; // class DefaultAlgorithm

} // namespace ROL

#endif
