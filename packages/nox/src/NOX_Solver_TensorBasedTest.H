// $Id$ 
// $Source$ 

#ifdef WITH_PRERELEASE
//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef NOX_SOLVER_BLOCKTENSORBASED_H
#define NOX_SOLVER_BLOCKTENSORBASED_H

#include "NOX_Solver_Generic.H"	    // base class
#include "NOX_LineSearch_Manager.H" // class data element
#include "NOX_Direction_Manager.H"  // class data element
#include "NOX_Parameter_List.H"	    // class data element
#include "NOX_Utils.H"	            // class data element

#include "NOX_LineSearch_Generic.H" // base class
#include "NOX_LineSearch_Utils_Printing.H"  // class data member
#include "NOX_LineSearch_Utils_Counters.H"  // class data member
#include "NOX_LineSearch_Utils_Slope.H"     // class data member


namespace NOX {
namespace Solver {


class TensorBasedTest : public Generic {

public:

  /*! 
    \brief Constructor

    See reset() for description. 
  */
  TensorBasedTest(NOX::Abstract::Group& xgrp, 
		   NOX::StatusTest::Generic& t, 
		   NOX::Parameter::List& p);

  //! Destructor
  virtual ~TensorBasedTest();

  virtual bool reset(NOX::Abstract::Group& xgrp, 
		     NOX::StatusTest::Generic& t, 
		     NOX::Parameter::List& p);
  virtual bool reset(NOX::Abstract::Group& xgrp, 
		     NOX::StatusTest::Generic& t);
  virtual NOX::StatusTest::StatusType getStatus();
  virtual NOX::StatusTest::StatusType iterate();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const NOX::Parameter::List& getParameterList() const;

protected:
  
  //! Print out initialization information.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

protected:

  //! Subroutine for computing the tensor and Newton directions
  bool NOX::Solver::TensorBasedTest::computeTensorDirection(
			 NOX::Abstract::Group& soln,
			 const NOX::Solver::Generic& solver);  

  //! Subroutine for calculating beta
  double NOX::Solver::TensorBasedTest::calculateBeta(double qa,
						      double qb,
						      double qc,
						      double& qval,
						      double& lambdaBar,
						      double lambda=1.0) const;

  //! Subroutine for computing the curvilinear step
  bool NOX::Solver::TensorBasedTest::computeCurvilinearStep(
				     NOX::Abstract::Vector& dir,
				     const NOX::Abstract::Group& soln,
				     const NOX::Solver::Generic& s,
				     double& lambda);

  //! Subroutine for executing the tensor linesearch
  bool NOX::Solver::TensorBasedTest::implementGlobalStrategy(
				      NOX::Abstract::Group& newGrp,
				      double& step,
				      const NOX::Solver::Generic& s);
    
  
  //! Performs a standard tensor linesearch (tensor or Newton direction)
  bool performLinesearch(NOX::Abstract::Group& newsoln,
			 double& step,
			 const NOX::Abstract::Vector& dir,
			 const NOX::Solver::Generic& s);

  //! Compute the residual norm of the local model
  double NOX::Solver::TensorBasedTest::getNormModelResidual(
				       const NOX::Abstract::Vector& dir,
				       const NOX::Abstract::Group& soln,
				       bool isTensorModel) const;

  //! Print pertinent information about the direction
  void NOX::Solver::TensorBasedTest::printDirectionInfo(char* dirName,
					  const NOX::Abstract::Vector& dir,
					  const NOX::Abstract::Group& soln,
					  bool isTensorModel) const;

  //! Calculate the directional derivative
  double NOX::Solver::TensorBasedTest::getDirectionalDerivative(
				       const NOX::Abstract::Vector& dir,
				       const NOX::Abstract::Group& soln) const;

  //! Select lambda for linesearch (quadratic or halving)
  double selectLambda(double newf, double oldf,
		      double oldfprime, double lambda);

  void NOX::Solver::TensorBasedTest::throwError(const string& functionName,
						 const string& errorMsg) const;
  
protected:
  
  //! Current solution.
  NOX::Abstract::Group* solnPtr;		

  //! Previous solution pointer. 
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Group* oldSolnPtr;	
  //! Previous solution reference.
  NOX::Abstract::Group& oldSoln;	

  //! Current Newton direction.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Vector* newtonVecPtr;
  //! Current Newton direction.reference.
  NOX::Abstract::Vector& newtonVec;

  //! Current tensor direction.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Vector* tensorVecPtr;
  //! Current tensor direction.reference.
  NOX::Abstract::Vector& tensorVec;

  //! Current tensor term vector.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Vector* acVecPtr;
  //! Current tensor term vector
  NOX::Abstract::Vector& acVec;

  //! Vector to previous point.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Vector* scVecPtr;
  //! Vector to previous point
  NOX::Abstract::Vector& scVec;

  //! Working vector.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Vector* tmpVecPtr;
  //! Working vector
  NOX::Abstract::Vector& tmpVec;

  //! Residual vector.pointer.
  NOX::Abstract::Vector* residualVecPtr;

  //! Stopping test.
  NOX::StatusTest::Generic* testPtr;		

  //! Input parameters.
  NOX::Parameter::List* paramsPtr;	

  //! Line Search parameters.
  NOX::Parameter::List& lsParams;	

  //! Direction parameters.
  NOX::Parameter::List& dirParams;	

  //! Linear Solver for the local model parameters.
  NOX::Parameter::List& localParams;	

  //! Printing Utils
  NOX::Utils utils;

  //! Current step.
  double step;			

  //! Value of sc'*dt.
  double beta;

  //! Number of nonlinear iterations.
  int nIter;                    

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

  //! Types of steps 
  enum StepType
    {
      TensorStep,
      NewtonStep
    };

  //! Flag for the base direction to compute after the first iteration 
  StepType requestedBaseStep;

  //! Flag for the direction to be computed this iteration 
  //stepType requestedStep;
  
  //! Enumerated list for each type of line search
  enum LineSearchType
    {
      Curvilinear,
      Standard,
      Dual,
      Newton
    };
  
  //! Choice of line search
  LineSearchType lsType;
  
  //! Algorithms used to determine convergence of the line search
  enum ConvergenceCriteriaType {
    //! Sufficient decrease condition
    ArmijoGoldstein, 
    //! Ared/Pred condition
    AredPred,
    //! Just accept the first step
    None
  };

  //! Choice of convergence criteria   (currently unused)
  ConvergenceCriteriaType convCriteria;

  //! Types of lambda selection 
  enum LambdaSelectionType
    {
      Halving,
      Quadratic
    };

  //! Flag for lambda selection (Halving/Quadratic)
  LambdaSelectionType lambdaSelection;

  //! Flag for Newton direction
  bool isNewtonDirection;

  //! Minimum step length (i.e., when we give up)
  double minStep;
 
  //! Default step
  double defaultStep;

  //! Default step for linesearch failure
  double recoveryStep;

  //! Maximum iterations
  int maxIters;

  //! Scaling factor for the Armijo-Goldstein condition
  double alpha;

  //! Value of sc'*inv(J)*F
  double sctjf;

  //! Value of sc'*inv(J)*a
  double sctja;

  //! Common line search printing utilities.
  NOX::LineSearch::Utils::Printing print;

  //! Common common counters for line searches.
  NOX::LineSearch::Utils::Counters counter;

  //! Common slope calculations for line searches.
  NOX::LineSearch::Utils::Slope slopeObj;

  //! Counter for number of Jacobian-vector products
  mutable int multsJv;

};
} // namespace Solver
} // namespace NOX

#endif

#endif  // WITH_PRERELEASE
