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

#ifndef NOX_SOLVER_TENSORBASED_H
#define NOX_SOLVER_TENSORBASED_H

#include "NOX_Solver_Generic.H"	    // base class
#include "NOX_LineSearch_Manager.H" // class data element
#include "NOX_Direction_Manager.H"  // class data element
#include "NOX_Parameter_List.H"	    // class data element
#include "NOX_Utils.H"	            // class data element

#include "NOX_LineSearch_Generic.H" // base class
#include "NOX_LineSearch_Utils_Printing.H"  // class data member
#include "NOX_LineSearch_Utils_Counters.H"  // class data member
#include "NOX_LineSearch_Utils_Slope.H"     // class data member


//Forward Declarations
namespace NOX {
  namespace Parameter {
    class PrePostOperator;
  }
}

namespace NOX {
namespace Solver {

/*! 
  \brief Nonlinear solver based on a rank-1 tensor method.
  
  Solves \f$F(x)=0\f$ using a rank-1 tensor method and a linesearch
  globalization.
  
  At the kth nonlinear iteration, the solver does the following:

  <ul>

  <li>Computes the tensor direction \f$ d_T \f$ by finding the root or
      smallest magnitude minimizer of the local model \f[ M_T(x_k+d) =
      F_k + J_kd + a_k(s_k^Td)^2, \f] where \f[ a_k = 2(F_{k-1} - F_k
      - J_ks_k) / (s_k^Ts_k)^2 \f] and \f[ s_k = s_{k-1} - s_k. \f]

  <li>Modifies the step according to a global strategy and updates
      \f$x\f$ as \f$x_{k+1} = x_k + d(\lambda) \f$ via a linesearch
      method, where \f$ d(\lambda) \f$ is some function of \f$ \lambda
      \f$.  For instance, the curvilinear step \f$ d_{\lambda T} \f$
      is a function of the linesearch parameter \f$ \lambda \f$ and is
      a parametric step that spans the directions of the tensor step
      and the Newton step.  At \f$ \lambda=1 \f$, the curvilinear step
      equals the full tensor step, and as \f$ \lambda \f$ nears 0, the
      curvilinear step approaches the Newton direction.  This step
      provides a monotonic decrease in the norm of the local tensor
      model as \f$ \lambda \f$ varies from 0 to 1.

  </ul>
  
  The solver iterates until the status tests (see NOX::StatusTest)
  determine either failure or convergence.

  <B>Input Parameters</B>

  To use this solver, set the "Nonlinear Solver" parameter to be
  "Tensor Based".  Then, specify the following sublists with the
  appropriate parameters as indicated below.

  <ul>

  <li> "Direction" - Sublist of the direction parameters, passed to
       the NOX::Direction::Manager constructor. Defaults to an empty
       list.

       <ul>

       <li> "Method" - Name of the direction to be computed in this
       solver.  "Tensor" and "Newton" are the only two valid choices.
       A sublist by that name specifies all of the parameters to
       be passed to the linear solver.  See below.

       <li> "Rescue Bad Newton Solve" (Boolean) - If the linear solve
       does not meet the tolerance specified by the forcing term, then
       use the step anyway. Defaults to true.

       <li> "Linear Solver" - Sublist for the specific linear solver
       parameters that are passed to
       NOX::Abstract::Group::computeNewton() and
       NOX::Abstract::Group::applyJacobianInverse().  "Linear Solver"
       is itself a sublist of the list specified in "Method" above.
       Below is a partial list of standard parameters usually
       available in common linear solvers.  Check with the specific
       linear solver being used for other parameters.

	    <ul>
     
	    <li> "Max Iterations" - Maximum number of Arnoldi
            iterations (also max Krylov space dimension)

            <li> "Tolerance" - Relative tolerance for solving local
            model [default = 1e-4]

	    <li> "Output Frequency" - Print output at every number of
            iterations [default = 20]

	    </ul>

       </ul> 

  <li> "Line Search" - Sublist of the line search parameters. Defaults
       to an empty list.  Be aware that the tensor step is not
       guaranteed to be a descent direction on the function.  Thus,
       not all "basic" line search approaches would be appropriate.
       This solver class approriately handles these considerations.
       The following parameters specify the specific options for this
       line search:

       <ul>

       <li> "Method" - Name of the line search available to tensor
       methods Valid choices are:

            <ul> 

	    <li> "Standard" - Backtrack along tensor direction unless
	    it is not a descent direction, in which case backtrack
	    along Newton direction.

	    <li> "Dual" - Backtrack along both the Newton and tensor
	    directions and choose the better of the two.

            <li> "Curvilinear" - Backtrack along the "curvilinear"
	    path that spans the tensor direction and the Newton
	    direction and that maintains monotonicity on the tensor
	    model. [Default]

            </ul>

       <li> "Lambda selection" - flag for how to calculate next lambda
       (quadratic, halving)

       <li> "Default Step" - starting value of linesearch parameter
       (defaults to 1.0)

       <li> "Minimum Step" - minimum acceptable linesearch parameter
       before the linesearch terminates (defaults to 1.0e-12)

       <li> "Recovery Step" - step parameter to take when the line
       search fails (defaults to value for "Default Step")

       <li> "Max Iters" - maximum number of iterations (i.e., backtracks)

       </ul>

  <li>"Solver Options" - Sublist of general solver options.  
       <ul>
       <li> "User Defined Pre/Post Operator" is supported.  See NOX::Parameter::PrePostOperator for more details.
       </ul>

  </ul>

  <B>Output Parameters</B>

  Every time solve() is called, a sublist for output parameters called
  "Output" will be created and will contain the following parameters:

  <ul>

  <li> "Nonlinear Iterations" - Number of nonlinear iterations
  
  <li> "2-Norm of Residual" - L-2 norm of the final residual \f$ F(x_k) \f$.

  </ul>

  
<b>References</b>

  <ul>

  <li> B. W. Bader, <em>Tensor-Krylov methods for solving large-scale
  systems of nonlinear equations</em>, Ph.D. Thesis, 2003, University
  of Colorado, Boulder, Colorado.

  <li> B. W. Bader, <em>Tensor-Krylov methods for solving large-scale
  systems of nonlinear equations</em>, submitted to SIAM J. Opt.

  <li> B. W. Bader and R. B. Schnabel, <em>Curvilinear linesearch for
  tensor methods</em>, SISC, 25(2):604-622.

  <li> R. B. Schnabel and P. D. Frank, <em>Tensor methods for
  nonlinear equations</em>, SIAM J. Numer. Anal., 21(5):815-843.

  </ul>


  \author Brett Bader (SNL 9233)
  
*/

class TensorBased : public Generic {

public:

  /*! 
    \brief Constructor

    See reset() for description. 
  */
  TensorBased(NOX::Abstract::Group& grp, 
		  NOX::StatusTest::Generic& tests, 
		  NOX::Parameter::List& params);

  //! Destructor
  virtual ~TensorBased();

  virtual bool reset(NOX::Abstract::Group& grp, 
		     NOX::StatusTest::Generic& tests, 
		     NOX::Parameter::List& params);
  virtual bool reset(NOX::Abstract::Group& xgrp, 
		     NOX::StatusTest::Generic& tests);
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
  bool NOX::Solver::TensorBased::computeTensorDirection(
			 NOX::Abstract::Group& soln,
			 const NOX::Solver::Generic& solver);  

  //! Subroutine for calculating beta
  double NOX::Solver::TensorBased::calculateBeta(double qa,
						 double qb,
						 double qc,
						 double& qval,
						 double& lambdaBar,
						 double lambda=1.0) const;

  //! Subroutine for computing the curvilinear step
  bool NOX::Solver::TensorBased::computeCurvilinearStep(
				     NOX::Abstract::Vector& dir,
				     const NOX::Abstract::Group& soln,
				     const NOX::Solver::Generic& s,
				     double& lambda);

  //! Subroutine for executing the tensor linesearch
  bool NOX::Solver::TensorBased::implementGlobalStrategy(
				      NOX::Abstract::Group& newGrp,
				      double& step,
				      const NOX::Solver::Generic& s);
    
  
  //! Performs a standard tensor linesearch (tensor or Newton direction)
  bool performLinesearch(NOX::Abstract::Group& newsoln,
			 double& step,
			 const NOX::Abstract::Vector& dir,
			 const NOX::Solver::Generic& s);

  //! Compute the residual norm of the local model
  double NOX::Solver::TensorBased::getNormModelResidual(
				       const NOX::Abstract::Vector& dir,
				       const NOX::Abstract::Group& soln,
				       bool isTensorModel) const;

  //! Print pertinent information about the direction
  void NOX::Solver::TensorBased::printDirectionInfo(char* dirName,
					  const NOX::Abstract::Vector& dir,
					  const NOX::Abstract::Group& soln,
					  bool isTensorModel) const;

  //! Calculate the directional derivative
  double NOX::Solver::TensorBased::getDirectionalDerivative(
				       const NOX::Abstract::Vector& dir,
				       const NOX::Abstract::Group& soln) const;

  //! Select lambda for linesearch (quadratic or halving)
  double selectLambda(double newf, double oldf,
		      double oldfprime, double lambda);

  //! Throw an error with a method's name and error message
  void NOX::Solver::TensorBased::throwError(const string& functionName,
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
  enum ConvergenceCriteriaType
  {
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

  //! Counter for number of "double" Jacobian-vector products
  mutable int mults2Jv;

  //! Pointer to a user defined NOX::Abstract::PrePostOperator object.
  NOX::Parameter::PrePostOperator* prePostOperatorPtr;

  //! True if a PrePostOperator was registered with the solver.
  bool havePrePostOperator;

};
} // namespace Solver
} // namespace NOX

#endif

#endif  // WITH_PRERELEASE
