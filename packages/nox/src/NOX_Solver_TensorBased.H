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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_SOLVER_TENSORBASED_H
#define NOX_SOLVER_TENSORBASED_H

#include "NOX_Solver_Generic.H"	    // base class
#include "NOX_LineSearch_Manager.H" // class data element
#include "NOX_Direction_Tensor.H"   // class data element
#include "NOX_Parameter_List.H"	    // class data element
#include "NOX_Abstract_Vector.H"    // for NOX::Abstract::Vector::NormType

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
   <li> "Curvilinear" - Backtrack along the "curvilinear" path that 
   spans the tensor direction and the Newton direction and that maintains
   monotonicity on the tensor model. [Default]
   </ul>
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
  TensorBased(NOX::Abstract::Group& xgrp, NOX::StatusTest::Generic& t, const NOX::Parameter::List& p);

  //! Destructor
  virtual ~TensorBased();

  virtual bool reset(NOX::Abstract::Group& xgrp, NOX::StatusTest::Generic& t, 
		     NOX::Parameter::List& p);
  virtual NOX::StatusTest::StatusType getStatus();
  virtual NOX::StatusTest::StatusType iterate();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const NOX::Parameter::List& getParameterList() const;

protected:
  
  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

  //! Sets the output parameters for the solver in an "Output" sublist
  bool setOutputParameters(); 

  //! Performs a standard tensor linesearch (tensor or Newton direction)
  bool performStandardLinesearch();

  //! Performs a curvilinear linesearch (like a trust region method)
  bool performCurvilinearLinesearch();

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
  NOX::Parameter::List params;	

  //! %Search %Direction. 
  NOX::Direction::Tensor direction; 

  //! Current step.
  double step;			

  //! Number of nonlinear iterations.
  int niter;

  //! Number of linesearch iterations (i.e., backtracks).
  int totalNumLSIterations;

  //! Number of times that a linesearch failed.
  int totalNumFailedLineSearches;

  //! Number of Newton/tensor steps that required a linesearch.
  int totalNumLSSteps;

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

  //! "alpha" parameter in Armijo-Goldstein sufficient decrease condition.
  double alpha;
  

private:  

  //! Minimum step length (i.e., when we give up)
  double minStep;
 
  //! Default step
  double defaultStep;

  //! Default step for linesearch failure
  double recoveryStep;

  //! Maximum number of iterations
  int maxIters;
  
};
} // namespace Solver
} // namespace NOX

#endif

#endif
