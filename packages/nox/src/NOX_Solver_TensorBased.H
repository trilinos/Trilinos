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
#include "NOX_Direction_Tensor.H"   // class data element
#include "NOX_Parameter_List.H"	    // class data element
#include "NOX_Utils.H"              // class data element

#include "NOX_LineSearch_Generic.H" // base class
#include "NOX_LineSearch_Utils_Printing.H"  // class data member
#include "NOX_LineSearch_Utils_Counters.H"  // class data member
#include "NOX_LineSearch_Utils_Slope.H"     // class data member


namespace NOX {
namespace LineSearch { 

class Tensor : public Generic {

public:

  //! Constructor
  Tensor(const NOX::Utils& u, NOX::Parameter::List& params);

  //! Destructor
  ~Tensor();
  
  // derived
  bool reset(NOX::Parameter::List& params);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step, 
	       const NOX::Abstract::Vector& dir,
	       const NOX::Solver::Generic& s);

protected:

  //! Performs a standard tensor linesearch (tensor or Newton direction)
  bool performLinesearch(NOX::Abstract::Group& newsoln,
			 double& step,
			 const NOX::Abstract::Vector& dir,
			 const NOX::Solver::Generic& s,
			 const NOX::Direction::Tensor& direction);

  //! Select lambda for linesearch (quadratic or halving)
  double selectLambda(double newf, double oldf,
		      double oldfprime, double lambda);

protected:

  //! Type of tensor linesearch
  enum LineSearchType {Curvilinear, Standard, Dual, Newton};
  
  //! Algorithms used to determine convergence of the line search
  enum ConvergenceCriteriaType {
    //! Sufficient decrease condition
    ArmijoGoldstein, 
    //! Ared/Pred condition
    AredPred,
    //! Just accept the first step
    None
  };

  //! Types of lambda selection 
  enum LambdaSelectionType {Halving, Quadratic};

  //! Flag for lambda selection (Halving/Quadratic)
  LambdaSelectionType lambdaSelection;

  //! Choice of convergence criteria   (currently unused)
  ConvergenceCriteriaType convCriteria;

  //! Choice of line search
  LineSearchType lsType;
  
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

  //! Pointer to the input parameter list.  We need this to create an "Output" sublist to store output parameters from the counter object.  
  NOX::Parameter::List* paramsPtr;

  //! Common line search printing utilities.
  NOX::LineSearch::Utils::Printing print;

  //! Common common counters for line searches.
  NOX::LineSearch::Utils::Counters counter;

  //! Common slope calculations for line searches.
  NOX::LineSearch::Utils::Slope slopeObj;

  //! Counter for number of Jacobian-vector products
  mutable int multsJv;

};
} // namespace LineSearch
} // namespace NOX



namespace NOX {
namespace Solver {

/*!
  \brief Linesearch-based solver that uses tensor methods.
  
  The following parameters are valid for this solver:
<ul>
<li>   "%Direction" - Sublist of the direction parameters, passed to the
       NOX::Direction::Tensor constructor. Defaults to an empty list.

<li> "Linear %Solver" - Sublist of "Direction" and contains the linear
     solver parameters, passed to
     NOX::Abstract::Group::computeNewton() and to
     NOX::Direction::Tensor::compute().  Valid choices used in tensor
     step calculation are listed below.  More valid choices are
     available for the first step when the Newton step is computed by
     computeNewton().

<ul>
<li>    "Max Iterations" - Maximum number of Arnoldi iterations (also 
        max Krylov space dimension)
<li>    "Output Frequency" - Print output at every number of iterations
        [default = 20]
<li>    "Preconditioning" - Name of preconditioning type.
<li>    "PreconditioningSide" - Side to apply preconditioning 
        (None, Left, Right)
<li>    "Reorthogonalize" - When to reorthogonalize Krylov basis vectors
        (None, As Needed, Always)
<li>    "Tolerance" - Relative tolerance for solving local model 
        [default = 1e-4]
</ul>
<br>

<li>   "Line Search" - Sublist of the line search parameters, passed to
       the NOX::LineSearch::Manager constructor. Defaults to an empty list.

    The following parameters can be specified for this line search:

<ul>
<li> "Method" - Name of the line search available to tensor methods
     Valid choices are
   <ul> 
   <li> "Standard" - Backtrack along tensor direction unless it is not a 
   descent direction, in which case backtrack along Newton direction.
   <li> "Dual" - Backtrack along both the Newton and tensor directions
   and choose the better of the two.
   <li> "Curvilinear" - Backtrack along the "curvilinear" path that 
   spans the tensor direction and the Newton direction and that maintains
   monotonicity on the tensor model. [Default]
   </ul>
<li> "Lambda selection" - flag for how to calculate next lambda
     (quadratic, halving)
<li> "Default Step" - starting step parameter (defaults to 1.0)
<li> "Minimum Step" - minimum acceptable step parameter (defaults to 1.0e-12)
<li> "Recovery Step" - step parameter to take when the line search fails
     (defaults to value for "Default Step")
<li> "Max Iters" - maximum number of iterations (i.e., backtracks)
</ul>
</ul>

    \author Brett W. Bader (CU-Boulder) bader&nbsp;\@&nbsp;cs.colorado.edu 

*/

class TensorBased : public Generic {

public:

  //! %TensorBased-method constructor. 
  TensorBased(NOX::Abstract::Group& grp,
	      NOX::StatusTest::Generic& tests,
	      NOX::Parameter::List& params);

  //! Destructor
  virtual ~TensorBased();

  virtual bool reset(NOX::Abstract::Group& grp,
		     NOX::StatusTest::Generic& tests,
		     NOX::Parameter::List& params);
  virtual bool reset(NOX::Abstract::Group& grp,
		     NOX::StatusTest::Generic& tests);
  virtual NOX::StatusTest::StatusType getStatus();
  virtual NOX::StatusTest::StatusType iterate();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const NOX::Parameter::List& getParameterList() const;

  //!  Return the direction object
  virtual const NOX::Direction::Tensor& getDirection() const;

protected:
  
  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

protected:
  
  //! Current solution.
  NOX::Abstract::Group* solnptr;		

  //! Previous solution pointer. 
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Group* oldsolnptr;	
  //! Previous solution reference.
  NOX::Abstract::Group& oldsoln;	

  //! Current search direction.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Vector* dirptr;
  //! Current search direction.reference.
  NOX::Abstract::Vector& dir;

  //! Stopping test.
  NOX::StatusTest::Generic* testptr;		

  //! Input parameters.
  NOX::Parameter::List* paramsPtr;	

  //! Utils
  NOX::Utils utils;
  
  //! Linesearch
  NOX::LineSearch::Tensor lineSearch;
  
  //! %Search %Direction. 
  NOX::Direction::Tensor direction;  // This is Direction::Manager in LSBased

  //! Current step.
  double step;			

  //! Number of nonlinear iterations.
  int niter;

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;
};
} // namespace Solver
} // namespace NOX

#endif

#endif
