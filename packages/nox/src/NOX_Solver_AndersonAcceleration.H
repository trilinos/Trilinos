// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_SOLVER_ANDERSON_ACCELERATION_H
#define NOX_SOLVER_ANDERSON_ACCELERATION_H

#include "NOX_Solver_Generic.H"             // base class
#include "Teuchos_ParameterList.hpp"             // class data element
#include "NOX_Utils.H"                 // class data element
#include "Teuchos_RCP.hpp"       // class data element
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

// Forward declarations
namespace NOX {
  class GlobalData;
  class Observer;
  namespace LineSearch {
    class Generic;
  }
}

namespace NOX {
namespace Solver {

/*!
  \brief Nonlinear solver based on Anderson Acceleration

  Solves \f$F(x)=0\f$ using Anderson Acceleration.

  This method requires a mixing parameter \f$\beta\f$ and a preconditioner
  \f$M(x)\f$. If preconditioning is not supported, it is taken to be the identity.

  For the first iteration, the solver computes \f$x_1 = x_0 + \beta M(x_0)F(x_0)\f$.

  Each iteration afterward, the solver does the following.

  <ul>
  <li> Set \f$m_k = \min\{m,k\}\f$, where \f$m\f$ is a parameter that determines the maximum number of iterates for
       which we store iterate and function values in addition to the current iterate. \f$m_k\f$ represents the number
       of previous iterates stored currently.

  <li> Solve the least-squares problem: \f$\min_{\gamma}\|\sum_{i=0}^{m_k} \alpha_iM(x_{k-i})F(x_{k-i}) \|\f$,
       such that \f$ \sum_{i=1}^{m_k} \alpha_i = 1 \f$. In this implementation, we solve an equivalent unconstrained
       problem for which the QR factors are stored and updated.

  <li> Calculate the new iterate: $x_{k+1} = \sum_{i=0}^{m_k} \alpha_i[x_{k-i} + \beta M(x_{k-i})F(x_{k-i})]$.

  </ul>

  The iterations progress until the status tests (see NOX::StatusTest) determine either
  failure or convergence.

  <B>Input Parameters</B>

  The following parameter list entries are valid for this solver:

  - "Anderson Parameters" - Sublist of parameters specific to AA

  - "Line Search" - Sublist of the line search parameters, passed to
    the NOX::LineSearch::Manager constructor. Defaults to an empty list.

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

  \author Alex Toth and Roger Pawlowski
*/

class AndersonAcceleration : public Generic {

public:

  //! Constructor
  /*!
    See reset(NOX::Abstract::Group&, NOX::StatusTest::Generic&, Teuchos::ParameterList&) for description
   */
  AndersonAcceleration(const Teuchos::RCP<NOX::Abstract::Group>& grp,
          const Teuchos::RCP<NOX::StatusTest::Generic>& tests,
          const Teuchos::RCP<Teuchos::ParameterList>& params);

  //! Destructor
  virtual ~AndersonAcceleration();

  virtual void reset();
  virtual void reset(const NOX::Abstract::Vector& initialGuess);
  virtual void reset(const NOX::Abstract::Vector& initialGuess,
                     const Teuchos::RCP<NOX::StatusTest::Generic>& tests);
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

  //! Print out initialization information and calcuate the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

  //! Updates QR factors by addition of a column
  void qrAdd(NOX::Abstract::Vector& newCol);

  //! Updates QR factors by deletion of a column
  void qrDelete();

  //! Reorthogonalized columns of Q and updates factorization
  void reorthogonalize();

protected:

  //! Pointer to the global data object.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Utils
  Teuchos::RCP<NOX::Utils> utilsPtr;

  //! Current solution.
  Teuchos::RCP<NOX::Abstract::Group> solnPtr;

  //! Previous solution pointer.
  Teuchos::RCP<NOX::Abstract::Group> oldSolnPtr;

  //! Stopping test.
  Teuchos::RCP<NOX::StatusTest::Generic> testPtr;

  //! Input parameters.
  Teuchos::RCP<Teuchos::ParameterList> paramsPtr;

  //! Linesearch.
  Teuchos::RCP<NOX::LineSearch::Generic> lineSearchPtr;

  //! Flag for preconditioning
  bool precond;

  //! Flag to recompute Jacobian each iterate. Needed if Jacobian required to apply preconditioner.
  bool recomputeJacobian;

  //! Vector for temporary storage
  Teuchos::RCP<NOX::Abstract::Vector> workVec;

  //! Current preconditioned function value
  Teuchos::RCP<NOX::Abstract::Vector> precF;

  //! Previous preconditioned function value
  Teuchos::RCP<NOX::Abstract::Vector> oldPrecF;

  //! Maximum storage depth
  int storeParam;

  //! Current storage depth
  int nStore;

  //! Mixing parameter
  double mixParam;

  //! Current step.
  double stepSize;

  //! "Matrix" of successive iterate differences
  std::vector< Teuchos::RCP<NOX::Abstract::Vector> > xMat;

  //! Q factor the for least-squares problem
  std::vector< Teuchos::RCP<NOX::Abstract::Vector> > qMat;

  //! R factor for the least-squares problem
  Teuchos::SerialDenseMatrix<int,double> rMat;

  //! Number of nonlinear iterations.
  int nIter;

  //! Frequency for reorthogonalization of QR factors. If nonpositive, we never reorthogonalize
  int orthoFrequency;

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

  //! Type of check to use for status tests.  See NOX::StatusTest for more details.
  NOX::StatusTest::CheckType checkType;

  //! Pointer to a user defined NOX::Observer object.
  Teuchos::RCP<NOX::Observer> observer;

  //! Scratch space for condition number estimation using LAPACK.
  std::vector<double> WORK;

  //! Scratch space for condition number estimation using LAPACK.
  std::vector<int> IWORK;

  //! If true, the QR matrix will be resized if the condiiton number is greater than the dropTolerance.
  bool adjustForConditionNumber;

  //! Value of the condition number above which the QR matrix is resized.
  double dropTolerance;

  //! The nonlinear iteration where Anderson Acceleration will start.
  int accelerationStartIteration;

  //! Disables the storage depth size check - used for unit testing corner cases only.
  bool disableChecksForUnitTesting;
};
} // namespace Solver
} // namespace NOX

#endif

