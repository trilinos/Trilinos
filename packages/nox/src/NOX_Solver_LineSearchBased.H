// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_SOLVER_LINESEARCHBASED_H
#define NOX_SOLVER_LINESEARCHBASED_H

#include "NOX_Solver_Generic.H"	         // base class
#include "NOX_LineSearch_Manager.H"      // class data element
#include "NOX_Direction_Manager.H"       // class data element
#include "NOX_Solver_PrePostOperator.H"  // class data element
#include "Teuchos_ParameterList.hpp"	         // class data element
#include "NOX_Utils.H"		         // class data element
#include "Teuchos_RefCountPtr.hpp"       // class data element

// Forward declarations
namespace NOX {
  class GlobalData;
}

namespace NOX {
namespace Solver {

/*!
  \brief Nonlinear solver based on a line search (i.e., damping)

  Solves \f$F(x)=0\f$ using an iterative line-search-based method.

  Each iteration, the solver does the following.

  <ul> 
  <li>Compute a search direction \f$d\f$ via a NOX::Direction method

  <li>Compute a step length \f$\lambda\f$ and update \f$x\f$ as \f$x_{\rm new}
      = x_{\rm old} + \lambda d\f$ via a NOX::LineSearch method.

  </ul>
  
  The iterations progress until the status tests (see NOX::StatusTest) determine either
  failure or convergence.

  \note To support several line searches and status tests, this
  version of the solver has a getStepSize() function that returns
  \f$\lambda\f$.

  <B>Input Parameters</B>

  The following parameter list entries are valid for this solver:
 
  - "Line Search" - Sublist of the line search parameters, passed to
    the NOX::LineSearch::Manager constructor. Defaults to an empty list.

  - "Direction" - Sublist of the direction parameters, passed to
    the NOX::Direction::Manager constructor. Defaults to an empty list.

  - "Solver Options" - Sublist of general solver options.  
     <ul>
     <li> "User Defined Pre/Post Operator" is supported.  See NOX::Parameter::PrePostOperator for more details.
     </ul>


  <B>Output Parameters</B>

  Every time solve() is called, a sublist for output parameters called
  "Output" will be created and contain the following parameters.

  "Output":

  - "Nonlinear Iterations" - Number of nonlinear iterations
  
  - "2-Norm of Residual" - Two-norm of final residual

  \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/

class LineSearchBased : public Generic {

public:

  //! Constructor
  /*!
    See reset(NOX::Abstract::Group&, NOX::StatusTest::Generic&, Teuchos::ParameterList&) for description
   */
  LineSearchBased(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
		  const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests, 
		  const Teuchos::RefCountPtr<Teuchos::ParameterList>& params);

  //! Destructor
  virtual ~LineSearchBased();

  virtual bool reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
		     const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests, 
		     const Teuchos::RefCountPtr<Teuchos::ParameterList>& params);
  virtual bool reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
		     const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests);
  virtual bool reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp);
  virtual NOX::StatusTest::StatusType getStatus();
  virtual NOX::StatusTest::StatusType step();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const Teuchos::ParameterList& getList() const;

  //! Return the line search step size from the current iteration
  virtual double getStepSize() const;

protected:
  
  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

protected:
  
  //! Pointer to the global data object.
  Teuchos::RefCountPtr<NOX::GlobalData> globalDataPtr;

  //! Utils
  Teuchos::RefCountPtr<NOX::Utils> utilsPtr;
  
  //! Current solution.
  Teuchos::RefCountPtr<NOX::Abstract::Group> solnPtr;		

  //! Previous solution pointer. 
  Teuchos::RefCountPtr<NOX::Abstract::Group> oldSolnPtr;	

  //! Current search direction pointer.
  Teuchos::RefCountPtr<NOX::Abstract::Vector> dirPtr;

  //! Stopping test.
  Teuchos::RefCountPtr<NOX::StatusTest::Generic> testPtr;		

  //! Input parameters.
  Teuchos::RefCountPtr<Teuchos::ParameterList> paramsPtr;	

  //! Linesearch. 
  NOX::LineSearch::Manager lineSearch; 

  //! %Search %Direction. 
  NOX::Direction::Manager direction; 

  //! Current step.
  double stepSize;			

  //! Number of nonlinear iterations.
  int nIter;                    

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

  //! Type of check to use for status tests.  See NOX::StatusTest for more details.
  NOX::StatusTest::CheckType checkType;

  //! Pointer to a user defined NOX::Abstract::PrePostOperator object.
  NOX::Solver::PrePostOperator prePostOperator;

};
} // namespace Solver
} // namespace NOX

#endif

