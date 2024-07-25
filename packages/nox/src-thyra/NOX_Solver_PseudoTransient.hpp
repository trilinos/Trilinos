// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_SOLVER_PSEUDO_TRANSIENT_HPP
#define NOX_SOLVER_PSEUDO_TRANSIENT_HPP

#include "NOX_Solver_Generic.H"             // base class
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp" // base class
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
  namespace Thyra {
    class Group;
  }
}

namespace Thyra {
  template <typename ScalarT> class VectorBase;
}

namespace NOX {
namespace Solver {

/*!
  \brief Pseudo-transient solver.

  Requires that the Thyra::ModelEvaluator implement transient support
  (uses alpha, beta, Jacobian evaluation and optionally x_dot in
  residual).

  Based on the 1998 Kelley Keyes paper, with minor modifications.
*/
  class PseudoTransient :
    public NOX::Solver::Generic,
    public Teuchos::ParameterListAcceptorDefaultBase
{

public:

  PseudoTransient(const Teuchos::RCP<NOX::Abstract::Group>& xGrp,
          const Teuchos::RCP<NOX::StatusTest::Generic>& tests,
          const Teuchos::RCP<Teuchos::ParameterList>& params);

  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  void reset(const NOX::Abstract::Vector& initialGuess,
         const Teuchos::RCP<NOX::StatusTest::Generic>& tests);
  void reset(const NOX::Abstract::Vector& initialGuess);
  void reset();
  NOX::StatusTest::StatusType step();
  NOX::StatusTest::StatusType solve();
  const NOX::Abstract::Group& getSolutionGroup() const;
  const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  NOX::StatusTest::StatusType getStatus() const;
  int getNumIterations() const;
  const Teuchos::ParameterList& getList() const;
  double getStepSize() const;

  Teuchos::RCP< const NOX::Abstract::Group > getSolutionGroupPtr() const;
  Teuchos::RCP< const NOX::Abstract::Group > getPreviousSolutionGroupPtr() const;
  Teuchos::RCP< const Teuchos::ParameterList > getListPtr() const;
  Teuchos::RCP<const NOX::SolverStats> getSolverStatistics() const;

protected:

  virtual void init();

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

  //! Group used to evaluate a transient residual
  Teuchos::RCP<NOX::Abstract::Group> transientResidualGroup;

  //! Current search direction pointer.
  Teuchos::RCP<NOX::Abstract::Vector> dirPtr;

  //! Stopping test.
  Teuchos::RCP<NOX::StatusTest::Generic> testPtr;

  //! Input parameters.
  //Teuchos::RCP<Teuchos::ParameterList> paramsPtr;

  //! Linesearch.
  Teuchos::RCP<NOX::LineSearch::Generic> lineSearchPtr;

  //! %Search %Direction.
  Teuchos::RCP<NOX::Direction::Generic> directionPtr;

  //! Current step.
  double stepSize;

  //! Number of nonlinear iterations.
  int nIter;

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

  //! Type of check to use for status tests.  See NOX::StatusTest for more details.
  NOX::StatusTest::CheckType checkType;

  //! Pointer to a user defined NOX::Observer object.
  Teuchos::RCP<NOX::Observer> observer;

  //! Pointer to solnPtr casted back to a thyra group
  Teuchos::RCP<NOX::Thyra::Group> thyraSolnGroup;
  //! Pointer to oldSolnPtr casted back to a thyra group
  Teuchos::RCP<NOX::Thyra::Group> thyraOldSolnGroup;
  //! Group used to evaluate a transient residual
  Teuchos::RCP<NOX::Thyra::Group> thyraTransientResidualGroup;

  //! Step size for pseudo-transient stepping
  double delta;
  //! Inverse step size for pseudo-transient stepping
  double inv_delta;
  //! Initial step size
  double deltaInit;
  //! Maximum step size
  double deltaMax;
  //! Minimum step size
  double deltaMin;
  //! Step size from previous iteration
  double deltaOld;
  //! Pseudo-transient time
  double time;
  //! solution time derivative used for scaling operator V in pseudo-transient paper
  Teuchos::RCP< ::Thyra::VectorBase<double> > x_dot;

  //! If set to true, the candidate direction will use the transient residual instead of the steady-state residual.  This is a modification of the Kelley-Keyes paper.
  bool use_transient_residual;

  //! Maximum number of iterations before pseudo-transient is disabled and the algorithm switches to a line search-based direct to steady state solve.
  int max_pseudo_transient_iterations;

  //! Parameters that are valid for this solver
  mutable Teuchos::RCP<Teuchos::ParameterList> validParameters;
};
} // namespace Solver
} // namespace NOX

#endif

