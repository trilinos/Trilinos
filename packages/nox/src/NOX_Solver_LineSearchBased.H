// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_SOLVER_LINESEARCHBASED_H
#define NOX_SOLVER_LINESEARCHBASED_H

#include "NOX_Solver_Generic.H"             // base class
#include "Teuchos_ParameterList.hpp"             // class data element
#include "NOX_Utils.H"                 // class data element
#include "Teuchos_RCP.hpp"       // class data element

// Forward declarations
namespace NOX {
  class GlobalData;
  class Observer;
  namespace LineSearch {
    class Generic;
  }
  namespace Direction {
    class Generic;
  }
}

namespace NOX {
namespace Solver {

/*!
  \brief Nonlinear solver based on a line search (i.e., damping)

  Solves \f$ F(x)=0\f$ using an iterative line-search-based method.

  Each iteration, the solver does the following.

  <ul>
  <li>Compute a search direction \f$ d\f$ via a NOX::Direction method

  <li>Compute a step length \f$\lambda\f$ and update \f$ x\f$ as \f$ x_{\rm new}
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
  LineSearchBased(const Teuchos::RCP<NOX::Abstract::Group>& grp,
          const Teuchos::RCP<NOX::StatusTest::Generic>& tests,
          const Teuchos::RCP<Teuchos::ParameterList>& params);

  //! Destructor
  virtual ~LineSearchBased();

  virtual void reset(const NOX::Abstract::Vector& initialGuess,
             const Teuchos::RCP<NOX::StatusTest::Generic>& tests);
  virtual void reset(const NOX::Abstract::Vector& initialGuess);
  virtual void reset();
  virtual NOX::StatusTest::StatusType step();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual NOX::StatusTest::StatusType getStatus() const;
  virtual int getNumIterations() const;
  virtual const Teuchos::ParameterList& getList() const;
  virtual double getStepSize() const;

  inline virtual Teuchos::RCP< const NOX::Abstract::Group > getSolutionGroupPtr() const {return solnPtr;};
  inline virtual Teuchos::RCP< const NOX::Abstract::Group > getPreviousSolutionGroupPtr() const {return oldSolnPtr;};
  inline virtual Teuchos::RCP< const Teuchos::ParameterList > getListPtr() const {return paramsPtr;};
  virtual Teuchos::RCP<const NOX::SolverStats> getSolverStatistics() const;

protected:

  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

protected:

  //! Pointer to the global data object.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Utils
  Teuchos::RCP<NOX::Utils> utilsPtr;

  //! Current solution.
  Teuchos::RCP<NOX::Abstract::Group> solnPtr;

  //! Previous solution pointer.
  Teuchos::RCP<NOX::Abstract::Group> oldSolnPtr;

  //! Current search direction pointer.
  Teuchos::RCP<NOX::Abstract::Vector> dirPtr;

  //! Stopping test.
  Teuchos::RCP<NOX::StatusTest::Generic> testPtr;

  //! Input parameters.
  Teuchos::RCP<Teuchos::ParameterList> paramsPtr;

  //! Linesearch.
  Teuchos::RCP<NOX::LineSearch::Generic> lineSearchPtr;

  //! %Search %Direction.
  Teuchos::RCP<NOX::Direction::Generic> directionPtr;

  //! Current step.
  double stepSize;

  //! Number of nonlinear iterations.
  int nIter;

  //! Determine if we want to catch and ignore exceptions thrown from the linear solver stack.
  bool catchThrowsDuringSolve;

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

  //! Type of check to use for status tests.  See NOX::StatusTest for more details.
  NOX::StatusTest::CheckType checkType;

  //! Pointer to a user defined NOX::Observer object.
  Teuchos::RCP<NOX::Observer> observer;
};
} // namespace Solver
} // namespace NOX

#endif

