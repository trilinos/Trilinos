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

#ifndef ROL_OPTIMIZATIONSOLVER_HPP
#define ROL_OPTIMIZATIONSOLVER_HPP

#include "ROL_Algorithm.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_CombinedStatusTest.hpp"

#include "Teuchos_oblackholestream.hpp"

/** \class ROL::OptimizationSolver
    \brief Provides a simplified interface for solving a wide range of
           optimization problems
 */

namespace ROL {

template<class Real>
class OptimizationSolver {
private:

  Teuchos::RCP<Algorithm<Real> >          algo_;
  Teuchos::RCP<Step<Real> >               step_;
  Teuchos::RCP<StatusTest<Real> >         status0_;
  Teuchos::RCP<CombinedStatusTest<Real> > status_;
  Teuchos::RCP<AlgorithmState<Real> >     state_;

  Teuchos::RCP<Vector<Real> > x_;
  Teuchos::RCP<Vector<Real> > g_;
  Teuchos::RCP<Vector<Real> > l_;
  Teuchos::RCP<Vector<Real> > c_;

  Teuchos::RCP<Objective<Real> >       obj_;
  Teuchos::RCP<BoundConstraint<Real> > bnd_;
  Teuchos::RCP<Constraint<Real> >      con_;

  std::vector<std::string>  output_;

  EProblem problemType_;

public:

  OptimizationSolver( OptimizationProblem<Real> &opt,
                      Teuchos::ParameterList &parlist ) {
    // Get step name from parameterlist
    std::string stepname = parlist.sublist("Step").get("Type","Last Type (Dummy)");
    EStep stepType = StringToEStep(stepname);

    // Get optimization problem type: U, E, B, EB
    problemType_ = opt.getProblemType();

    // Set default algorithm if provided step is incompatible with problem type
    if ( !isCompatibleStep(problemType_, stepType) ) {
      switch ( problemType_ ) {
        case TYPE_U:
          stepType = STEP_TRUSTREGION;          break;
        case TYPE_B:
          stepType = STEP_TRUSTREGION;          break;
        case TYPE_E:
          stepType = STEP_COMPOSITESTEP;        break;
        case TYPE_EB:
          stepType = STEP_AUGMENTEDLAGRANGIAN;  break;
        case TYPE_LAST:
        default:
          throw Exception::NotImplemented(">>> ROL::OptimizationSolver: Unknown problem type!");
      }
    }
    stepname = EStepToString(stepType);

    // Build algorithmic components
    StepFactory<Real>        stepFactory;
    StatusTestFactory<Real>  statusTestFactory;
    step_    = stepFactory.getStep(stepname,parlist);
    status0_ = statusTestFactory.getStatusTest(stepname,parlist);
    status_  = Teuchos::rcp( new CombinedStatusTest<Real> );
    state_   = Teuchos::rcp( new AlgorithmState<Real> );

    // Get optimization vector and a vector for the gradient
    x_ = opt.getSolutionVector();
    g_ = x_->dual().clone();

    // If there is an equality constraint, get the multiplier and create a constraint vector
    if( problemType_ == TYPE_E || problemType_ == TYPE_EB ) {
      l_ = opt.getMultiplierVector();
      c_ = l_->dual().clone();
    }
    // Create modified objectives if needed
    if( stepType == STEP_AUGMENTEDLAGRANGIAN ) {
      Teuchos::RCP<Objective<Real> > raw_obj = opt.getObjective();
      con_ = opt.getConstraint();
      // TODO: Provide access to change initial penalty
      obj_ = Teuchos::rcp( new AugmentedLagrangian<Real>(raw_obj,con_,*l_,1.0,*x_,*c_,parlist) );
      bnd_ = opt.getBoundConstraint();
    }
    else if( stepType == STEP_MOREAUYOSIDAPENALTY ) {
      Teuchos::RCP<Objective<Real> > raw_obj = opt.getObjective();
      bnd_ = opt.getBoundConstraint();
      con_ = opt.getConstraint();
      // TODO: Provide access to change initial penalty
      obj_ = Teuchos::rcp( new MoreauYosidaPenalty<Real>(raw_obj,bnd_,*x_,parlist) );
    }
    else if( stepType == STEP_INTERIORPOINT ) {
      Teuchos::RCP<Objective<Real> > raw_obj = opt.getObjective();
      bnd_ = opt.getBoundConstraint();
      con_ = opt.getConstraint();
      // TODO: Provide access to change initial penalty
      obj_ = Teuchos::rcp( new InteriorPoint::PenalizedObjective<Real>(raw_obj,bnd_,*x_,parlist) );
    }
    else {
      obj_   = opt.getObjective();
      bnd_   = opt.getBoundConstraint();
      con_   = opt.getConstraint();
    }
  }

  std::vector<std::string> getOutput(void) const {
    return output_;
  }

  int solve(void) {
    Teuchos::oblackholestream bhs;
    return solve(bhs);
  }

  int solve( std::ostream &outStream,
       const Teuchos::RCP<StatusTest<Real> > &status = Teuchos::null,
       const bool combineStatus = true ) {
    // Build algorithm
    status_->reset();       // Clear previous status test
    status_->add(status0_); // Default StatusTest
    if (status != Teuchos::null) {
      if (!combineStatus) { // Use only user-defined StatusTest
        status_->reset();
      }
      status_->add(status); // Add user-defined StatusTest
    }
    algo_ = Teuchos::rcp( new Algorithm<Real>( step_, status_, state_ ) );

    switch(problemType_) {
      case TYPE_U:
        output_ = algo_->run(*x_,*g_,*obj_,true,outStream);
      break;
      case TYPE_B:
        output_ = algo_->run(*x_,*g_,*obj_,*bnd_,true,outStream);
      break;
      case TYPE_E:
        output_ = algo_->run(*x_,*g_,*l_,*c_,*obj_,*con_,true,outStream);
      break;
      case TYPE_EB:
        output_ = algo_->run(*x_,*g_,*l_,*c_,*obj_,*con_,*bnd_,true,outStream);
      break;
      case TYPE_LAST:
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          "Error in OptimizationSolver::solve() : Unsupported problem type");
      break;
    }

    // TODO: Interrogate AlgorithmState and StatusTest to generate a return code
    //       that indicates why the solver has stopped

    // Return an integer code
    return 0;
  }

  Teuchos::RCP<const AlgorithmState<Real> > getAlgorithmState(void) const {
    return state_;
  }

  void resetAlgorithmState(void) {
    state_ = Teuchos::rcp( new AlgorithmState<Real> );
  }

}; // class OptimizationSolver

} // namespace ROL

#endif // ROL_OPTIMIZATIONSOLVER_HPP


