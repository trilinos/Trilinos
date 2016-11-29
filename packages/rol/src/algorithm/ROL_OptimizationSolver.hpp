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

#define OPTIMIZATION_PROBLEM_REFACTOR

#include "ROL_Algorithm.hpp"
#include "ROL_OptimizationProblemRefactor.hpp"

#include "Teuchos_oblackholestream.hpp"

/** \class ROL::OptimizationSolver 
    \brief Provides a simplified interface for solving a wide range of
           optimization problems
 */

namespace ROL {

template<class Real>
class OptimizationSolver {

  typedef Vector<Real>               V;
  typedef Objective<Real>            OBJ;
  typedef EqualityConstraint<Real>   EQCON;
  typedef BoundConstraint<Real>      BND;

private:

  Teuchos::RCP<Algorithm<Real> >             algo_;
  Teuchos::RCP<Step<Real> >                  step_;
  Teuchos::RCP<StatusTest<Real> >            status_;
  Teuchos::RCP<AlgorithmState<Real> >        state_; 

  Teuchos::RCP<V>         x_;
  Teuchos::RCP<V>         g_;
  Teuchos::RCP<V>         l_;
  Teuchos::RCP<V>         c_;

  Teuchos::RCP<OBJ>       obj_;
  Teuchos::RCP<BND>       bnd_;
  Teuchos::RCP<EQCON>     eqcon_;

  std::vector<std::string>  output_;
  
  EProblem problemType_;

public:

  OptimizationSolver( OptimizationProblem<Real> &opt,
                      Teuchos::ParameterList &parlist ) {

    using Teuchos::RCP; using Teuchos::rcp;

    std::string stepname = parlist.sublist("Step").get("Type","Last Type (Dummy)");

    EStep stepType = StringToEStep(stepname);

    problemType_ = opt.getProblemType();  

    TEUCHOS_TEST_FOR_EXCEPTION( !isValidStep(stepType), std::invalid_argument,
                                "Invalid step name in OptimizationSolver constructor!" );
                      
    TEUCHOS_TEST_FOR_EXCEPTION( !isCompatibleStep(problemType_, stepType), std::logic_error,
      "Error in OptimizationSolver constructor: Step type " << stepname << " does not support "
      << EProblemToString(problemType_) << " problems" << std::endl ); 

    StepFactory<Real>        stepFactory;
    StatusTestFactory<Real>  statusTestFactory;

    step_     = stepFactory.getStep(stepname,parlist);
    status_   = statusTestFactory.getStatusTest(stepname,parlist);
    state_    = rcp( new AlgorithmState<Real> );

    algo_ = rcp( new Algorithm<Real>( step_, status_, state_ ) );

    x_ = opt.getSolutionVector();
    g_ = x_->dual().clone();

    // If there is an equality constraint, get the multiplier and create a constraint vector
    if( problemType_ == TYPE_E || problemType_ == TYPE_EB ) {
      l_ = opt.getMultiplierVector();
      c_ = l_->dual().clone();
    }

    // Create modified objectives if needed
    if( stepType == STEP_AUGMENTEDLAGRANGIAN ) {
      RCP<OBJ> raw_obj = opt.getObjective();
      eqcon_ = opt.getEqualityConstraint();    

      // TODO: Provide access to change initial penalty
      obj_ = rcp( new AugmentedLagrangian<Real>(raw_obj,eqcon_,*l_,1.0,*x_,*c_,parlist) );
      bnd_ = opt.getBoundConstraint();
    }
    else if( stepType == STEP_MOREAUYOSIDAPENALTY ) {
      RCP<OBJ> raw_obj = opt.getObjective();
      bnd_ = opt.getBoundConstraint();        
      eqcon_ = opt.getEqualityConstraint(); 

      // TODO: Provide access to change initial penalty
      obj_ = rcp( new MoreauYosidaPenalty<Real>(raw_obj,bnd_,*x_,1.0) );
    }
    else {
      obj_   = opt.getObjective();   
      bnd_   = opt.getBoundConstraint();
      eqcon_ = opt.getEqualityConstraint();
    }
  }

  virtual std::vector<std::string> getOutput() {
    return output_;
  }

  virtual int solve() {
    Teuchos::oblackholestream bhs;
    return solve(bhs);
  }

  virtual int solve( std::ostream &outStream ) {

    switch(problemType_) {
      case TYPE_U:
        output_ = algo_->run(*x_,*g_,*obj_,true,outStream);    
      break;
      case TYPE_B:
        output_ = algo_->run(*x_,*g_,*obj_,*bnd_,true,outStream);
      break;
      case TYPE_E:
        output_ = algo_->run(*x_,*g_,*l_,*c_,*obj_,*eqcon_,true,outStream);
      break;
      case TYPE_EB:
        output_ = algo_->run(*x_,*g_,*l_,*c_,*obj_,*eqcon_,*bnd_,true,outStream);
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

  Teuchos::RCP<const AlgorithmState<Real> > getAlgorithmState() {
    return state_;
  }

 

}; // class OptimizationSolver

} // namespace ROL

#endif // ROL_OPTIMIZATIONSOLVER_HPP


