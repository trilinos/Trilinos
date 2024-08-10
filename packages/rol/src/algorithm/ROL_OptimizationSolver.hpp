// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OPTIMIZATIONSOLVER_HPP
#define ROL_OPTIMIZATIONSOLVER_HPP

#include "ROL_Algorithm.hpp"
#include "ROL_StepFactory.hpp"
#include "ROL_StatusTestFactory.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_CombinedStatusTest.hpp"

#include "ROL_Stream.hpp"

/** \class ROL::OptimizationSolver
    \brief Provides a simplified interface for solving a wide range of
           optimization problems
 */

namespace ROL {

template<class Real>
class OptimizationSolver {
private:

  ROL::Ptr<Algorithm<Real> >          algo_;
  ROL::Ptr<Step<Real> >               step_;
  ROL::Ptr<StatusTest<Real> >         status0_;
  ROL::Ptr<CombinedStatusTest<Real> > status_;
  ROL::Ptr<AlgorithmState<Real> >     state_;

  ROL::Ptr<Vector<Real> > x_;
  ROL::Ptr<Vector<Real> > g_;
  ROL::Ptr<Vector<Real> > l_;
  ROL::Ptr<Vector<Real> > c_;

  ROL::Ptr<Objective<Real> >       obj_;
  ROL::Ptr<BoundConstraint<Real> > bnd_;
  ROL::Ptr<Constraint<Real> >      con_;

  std::vector<std::string>  output_;

  EProblem problemType_;
  EStep stepType_;
  std::string stepname_;

  Real pen_;

public:

  /** \brief Constructor.
  
       @param[in] opt       the OptimizationProblem to be solved
       @param[in] parlist   algorithm and step input parameters

      ---
  */
  OptimizationSolver( OptimizationProblem<Real> &opt,
                      ROL::ParameterList &parlist ) {

    // Get optimization problem type: U, E, B, EB
    problemType_ = opt.getProblemType();

    // Initialize AlgorithmState
    state_ = ROL::makePtr<AlgorithmState<Real>>();

    // Get step name from parameterlist
    stepname_ = parlist.sublist("Step").get("Type","Last Type (Dummy)");
    stepType_ = StringToEStep(stepname_);

    // Set default algorithm if provided step is incompatible with problem type
    if ( !isCompatibleStep(problemType_, stepType_) ) {
      switch ( problemType_ ) {
        case TYPE_U:
          stepType_ = STEP_TRUSTREGION;          break;
        case TYPE_B:
          stepType_ = STEP_TRUSTREGION;          break;
        case TYPE_E:
          stepType_ = STEP_COMPOSITESTEP;        break;
        case TYPE_EB:
          stepType_ = STEP_AUGMENTEDLAGRANGIAN;  break;
        case TYPE_LAST:
        default:
          throw Exception::NotImplemented(">>> ROL::OptimizationSolver: Unknown problem type!");
      }
    }
    stepname_ = EStepToString(stepType_);

    // Build status test
    StatusTestFactory<Real> statusTestFactory;
    status0_ = statusTestFactory.getStatusTest(stepname_,parlist);
    status_  = ROL::makePtr<CombinedStatusTest<Real>>();

    // Get optimization vector and a vector for the gradient
    x_ = opt.getSolutionVector();
    g_ = x_->dual().clone();

    // Initialize Step
    StepFactory<Real> stepFactory;
    step_ = stepFactory.getStep(stepname_,parlist);

    // If there is an equality constraint, get the multiplier and create a constraint vector
    if( problemType_ == TYPE_E || problemType_ == TYPE_EB ) {
      l_ = opt.getMultiplierVector();
      c_ = l_->dual().clone();
    }

    // Create modified objectives if needed
    const Real one(1), ten(10);
    if( stepType_ == STEP_AUGMENTEDLAGRANGIAN ) {
      ROL::Ptr<Objective<Real> > raw_obj = opt.getObjective();
      con_ = opt.getConstraint();
      // TODO: Provide access to change initial penalty
      obj_ = ROL::makePtr<AugmentedLagrangian<Real>>(raw_obj,con_,*l_,1.0,*x_,*c_,parlist);
      bnd_ = opt.getBoundConstraint();
      pen_ = parlist.sublist("Step").sublist("Augmented Lagrangian").get("Initial Penalty Parameter",ten);
    }
    else if( stepType_ == STEP_MOREAUYOSIDAPENALTY ) {
      ROL::Ptr<Objective<Real> > raw_obj = opt.getObjective();
      bnd_ = opt.getBoundConstraint();
      con_ = opt.getConstraint();
      // TODO: Provide access to change initial penalty
      obj_ = ROL::makePtr<MoreauYosidaPenalty<Real>>(raw_obj,bnd_,*x_,parlist);
      pen_ = parlist.sublist("Step").sublist("Moreau-Yosida Penalty").get("Initial Penalty Parameter",ten);
    }
    else if( stepType_ == STEP_INTERIORPOINT ) {
      ROL::Ptr<Objective<Real> > raw_obj = opt.getObjective();
      bnd_ = opt.getBoundConstraint();
      con_ = opt.getConstraint();
      // TODO: Provide access to change initial penalty
      obj_ = ROL::makePtr<InteriorPoint::PenalizedObjective<Real>>(raw_obj,bnd_,*x_,parlist);
      pen_ = parlist.sublist("Step").sublist("Interior Point").get("Initial Barrier Parameter",ten);
    }
    else if( stepType_ == STEP_FLETCHER ) {
      ROL::Ptr<Objective<Real> > raw_obj = opt.getObjective();
      bnd_ = opt.getBoundConstraint();
      con_ = opt.getConstraint();
      if( bnd_->isActivated() ) {
        obj_ = ROL::makePtr<BoundFletcher<Real> >(raw_obj,con_,bnd_,*x_,*c_,parlist);
      }
      else {
        obj_ = ROL::makePtr<Fletcher<Real> >(raw_obj,con_,*x_,*c_,parlist);
      }
      pen_ = parlist.sublist("Step").sublist("Fletcher").get("Penalty Parameter",one);
    }
    else {
      obj_   = opt.getObjective();
      bnd_   = opt.getBoundConstraint();
      con_   = opt.getConstraint();
      if( stepType_ == STEP_TRUSTREGION ) {
        pen_ = parlist.sublist("Step").sublist("Trust Region").get("Initial Radius",ten);
      }
      else if( stepType_ == STEP_BUNDLE ) {
        pen_ = parlist.sublist("Step").sublist("Bundle").get("Initial Trust-Region Parameter",ten);
      }
    }
  }

  /** \brief Returns iteration history as a vector of strings.

      ---
  */
  std::vector<std::string> getOutput(void) const {
    return output_;
  }

  /** \brief Solve optimization problem with no iteration output.

      @param[in] status          is a user-defined StatusTest
      @param[in] combineStatus   if true, the user-defined StatusTest will be combined with the default StatusTest

      ---
  */
  int solve(const ROL::Ptr<StatusTest<Real> > &status = ROL::nullPtr,
            const bool combineStatus = true) {
    ROL::nullstream bhs;
    return solve(bhs,status,combineStatus);
  }

  /** \brief Solve optimization problem.

      @param[in] outStream       is the output stream to collect iteration history
      @param[in] status          is a user-defined StatusTest
      @param[in] combineStatus   if true, the user-defined StatusTest will be combined with the default StatusTest

      ---
  */
  int solve( std::ostream &outStream,
       const ROL::Ptr<StatusTest<Real> > &status = ROL::nullPtr,
       const bool combineStatus = true ) {
    // Build algorithm
    status_->reset();       // Clear previous status test
    status_->add(status0_); // Default StatusTest
    if (status != ROL::nullPtr) {
      if (!combineStatus) { // Use only user-defined StatusTest
        status_->reset();
      }
      status_->add(status); // Add user-defined StatusTest
    }
    algo_ = ROL::makePtr<Algorithm<Real>>( step_, status_, state_ );

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
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          "Error in OptimizationSolver::solve() : Unsupported problem type");
    }

    // TODO: Interrogate AlgorithmState and StatusTest to generate a return code
    //       that indicates why the solver has stopped

    // Return an integer code
    return 0;
  }

  /** \brief Return the AlgorithmState.

      ---
  */
  ROL::Ptr<const AlgorithmState<Real> > getAlgorithmState(void) const {
    return state_;
  }

  /** \brief Reset the AlgorithmState.

      This function does not reset the Step or the StepState.

      ---
  */
  void resetAlgorithmState(void) {
    state_ = ROL::makePtr<AlgorithmState<Real>>();
  }

  /** \brief Reset both Algorithm and Step.

      @param[in] resetAlgo   if true, then AlgorithmState will be reset

      This function will reset the AlgorithmState and reinitialize the
      Step.  This function does not permit changing the Step specified
      upon construction.  To change the Step, reinitialize the
      OptimizationSolver.

      ---
  */
  void reset(const bool resetAlgo = true) {
    // Reset AlgorithmState
    if (resetAlgo) {
      resetAlgorithmState();
    }
    // Reset StepState
    step_->reset(pen_);
    // Reset penalty objectives
    if( stepType_ == STEP_AUGMENTEDLAGRANGIAN ) {
      ROL::dynamicPtrCast<AugmentedLagrangian<Real> >(obj_)->reset(*l_,pen_);
    }
    else if( stepType_ == STEP_MOREAUYOSIDAPENALTY ) {
      ROL::dynamicPtrCast<MoreauYosidaPenalty<Real> >(obj_)->reset(pen_);
    }
    else if( stepType_ == STEP_INTERIORPOINT ) {
      ROL::dynamicPtrCast<InteriorPoint::PenalizedObjective<Real> >(obj_)->updatePenalty(pen_);
    }
  }

  /** \brief Grab step name (after check for consistency).

      @param[out] stepname   Name of step

      This function returns the algorithmic step name.

      ---
  */
  std::string getStepName(void) const {
    return stepname_;
  }

}; // class OptimizationSolver


template<typename Real>
inline Ptr<OptimizationSolver<Real>>
make_OptimizationSolver( OptimizationProblem<Real>& opt, 
                         ParameterList& parlist ) {
  return makePtr<OptimizationSolver<Real>>(opt,parlist);
}

template<typename Real>
inline Ptr<OptimizationSolver<Real>>
make_OptimizationSolver( const Ptr<OptimizationProblem<Real>>& opt, 
                         ParameterList& parlist ) {
  return makePtr<OptimizationSolver<Real>>(*opt,parlist);
}

template<typename Real>
inline Ptr<OptimizationSolver<Real>>
make_OptimizationSolver( OptimizationProblem<Real>& opt, 
                         const Ptr<ParameterList>& parlist ) {
  return makePtr<OptimizationSolver<Real>>(opt,*parlist);
}

template<typename Real>
inline Ptr<OptimizationSolver<Real>>
make_OptimizationSolver( const Ptr<OptimizationProblem<Real>>& opt, 
                         const Ptr<ParameterList>& parlist ) {
  return makePtr<OptimizationSolver<Real>>(*opt,*parlist);
}

} // namespace ROL

#endif // ROL_OPTIMIZATIONSOLVER_HPP


