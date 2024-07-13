// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_SOLVER_SINGLESTEP_H
#define NOX_SOLVER_SINGLESTEP_H

#include "NOX_Solver_Generic.H"  // for NOX::StatusTest::StatusType
#include "NOX_Utils.H"                 // class data element
#include "NOX_Abstract_Group.H"        // argument type

// Forward declarations
namespace NOX {
  class GlobalData;
  class Observer;
}

namespace NOX {
namespace Solver {

/*!
  \brief Light wrapper over linear solver, executes exactly one iteration.

  This class exists to support applications which always instantiate a NOX
  Nonlinear solver via the NOX::Solver::Factory, but may be used to
  solve known linear problems.
  In that case, it is convenient to create a NOX::Solver object which requires
  no status tests and simply executes the linear solver.

  Valid parameters in the ParameterList are:
  \param [in] "Ignore Linear Solver Failure" (bool) - If set to true, ignore the linear solver status. This is used for iterative methods where the requested tolerance was not explicitly met, but the solve might still be good enough to use. Defaults to false.
  \param [in] "Update Jacobian" - If set to true, the Jacobian is computed at the beginning of each nonlinear iteration. In some use cases, we may want to reuse the Jacobian matrix (i.e. only the RHS is changing). Defaults to true.
  \param [in] "Print Norms" - If set to true, the step output will print the residual norm and update norm.
  \param [in] "Print Relative Norm" - If set to true, this will cache the initial residual norm and if "Print Norms" is enabled, will also print the ratio of initial to final norm for the step.
 */
class SingleStep : public Generic {

public:

  SingleStep(const Teuchos::RCP<NOX::Abstract::Group>& grp,
         const Teuchos::RCP<Teuchos::ParameterList>& params);
  virtual ~SingleStep();
  virtual void reset();
  virtual void reset(const NOX::Abstract::Vector& initial_guess);
  virtual void reset(const NOX::Abstract::Vector& initial_guess,
             const Teuchos::RCP<NOX::StatusTest::Generic>& test);
  virtual NOX::StatusTest::StatusType step();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual Teuchos::RCP< const NOX::Abstract::Group > getSolutionGroupPtr() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual Teuchos::RCP< const NOX::Abstract::Group > getPreviousSolutionGroupPtr() const;
  virtual NOX::StatusTest::StatusType getStatus() const;
  virtual int getNumIterations() const;
  virtual const Teuchos::ParameterList& getList() const;
  virtual Teuchos::RCP< const Teuchos::ParameterList > getListPtr() const;
  virtual Teuchos::RCP<const NOX::SolverStats> getSolverStatistics() const;

protected:

  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints some status information.
  virtual void printUpdate();

  //! Prints some status information.
  bool check(Abstract::Group::ReturnType ret, const std::string& task);

  //! Prints some status information.
  bool try_step();

protected:

  //! Pointer to the global data object.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Utils
  Teuchos::RCP<NOX::Utils> utilsPtr;

  //! Current solution.
  Teuchos::RCP<NOX::Abstract::Group> solnPtr;

  //! Previous solution pointer.
  Teuchos::RCP<NOX::Abstract::Group> oldSolnPtr;

  //! Input parameters.
  Teuchos::RCP<Teuchos::ParameterList> paramsPtr;

  //! Number of nonlinear iterations.
  int nIter;

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

  //! Pointer to a user defined NOX::Observer object.
  Teuchos::RCP<NOX::Observer> observer;

  //! If set to true, the return status of the linear solver is ignored, step is declared a success
  bool ignoreLinearSolverFailures;

  //! If set to true, the Jacobian is computed at the beginning of each nonlinear iteration
  bool updateJacobian;

  //! Print the norms of residual (F) and update (dx).
  bool printNorms;

  //! If print norms is enabled, this will also compute the relative norm, \f$||F||/||F_0||\f$.
  bool computeRelativeNorm;

  //! Initial residual norm used for computing the relative norm.
  double normF_0;

  //! Group that contains the "frozen Jacobian" if updateJacobian is false.
  Teuchos::RCP<NOX::Abstract::Group> frozenJacobianPtr;
};
} // namespace Solver
} // namespace NOX

#endif

