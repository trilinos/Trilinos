// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_STATUSTEST_NORM_UPDATE_H
#define NOX_STATUSTEST_NORM_UPDATE_H

#include "NOX_StatusTest_Generic.H"    // base class
#include "NOX_Abstract_Vector.H"        // for NormType
#include "Teuchos_RCP.hpp"

// Forward declaration
namespace NOX {
namespace Abstract {
class Group;
}
}

namespace NOX {

namespace StatusTest {

/*!
  \brief Various convergence tests based on the norm of the change in
  the solution vector, \f$ x \f$, between outer iterations.

  If checkStatusEfficiently is called with the type set to NOX::StatusTest::None,
  then the status is set to NOX::StatusTest::Unevaluated and returned.
  (Also #normUpdate is set to -1.0.)

  If checkStatusEfficiently is called on the first iteration, then the status is
  set to NOX::StatusTest::Unconverged and returned.  (Also #normUpdate is
  set to -1.0.)

  If checkStatusEfficiently is called on a problem where the solution group does
  not have F evaluated (i.e., problem.getSolutionGroup().isF() is
  false), then the status is set to NOX::StatusTest::Unconverged and
  returned.   (Also #normUpdate is set to -1.0.)


  Finally, we return NOX::StatusTest::Converged if \f$\alpha <
  \beta\f$ and NOX::StatusTest::Uncoverged otherwise.
  Here
  \f$\alpha\f$ represents the norm of \f$ \Delta x \f$ and \f$\beta\f$
  represents the tolerance.  We define:

  \f[ \Delta x = x_k - x_{k-1} \f]

  where \f$ x_k \f$ is the solution vector of the \f$ k \f$-th nonlinear iterate.

  Let \f$\gamma\f$ denote an optional scale factor defined as

  - \f$\gamma = \frac{1}{n}\f$ if \c sType in the constructor is NOX::NormF::Scaled, and

  - \f$\gamma = 1\f$ if \c sType in the constructor is NOX::NormF::Unscaled.

  Then \f$\alpha\f$ is defined as follows:

  - If \c nType in the constructor is Abstract::Vector::TWO, then
    \f[ \alpha = \sqrt{ \gamma \sum_{i=1}^n \Delta x_i^2 } \f]

  - If \c nType in the constructor is Abstract::Vector::ONE, then
    \f[ \alpha = \gamma \sum_{i=1}^n | \Delta x_i | \f]

  - If \c nType in the constructor is Abstract::Vector::INF, then
    \f[ \alpha = \gamma \max_{i} | \Delta x_i |  \f]

  Finally, \f$\beta\f$ is set to the \c tolerance in the constructor, i.e.,

  \f[  \beta  = \mbox{tolerance} \f]
*/

/*
  NOTE: If no iterations have been performed, then there is no
  information for\f$ x_{k-1} \f$.  If the test is called while the
  number of iterations is zero (i.e. before an iteration takes place),
  the code returns the status test as Unconverged and sets the value
  to -1.0.  Once an iteration has taken place the Norm of the update
  is calculated as above.

*/
class NormUpdate : public Generic {

public:

  //! Type that determines whether to scale the norm by the problem size.
  enum ScaleType {
    //! No norm scaling
    Unscaled,
    //! Scale the norm by the length of the vector
    Scaled
  };

  //! Constructor for absolute norm.
  /*! This constructor defaults to the \c Absolute tolerance type. */
  NormUpdate(double tolerance, NOX::Abstract::Vector::NormType ntype, ScaleType stype = Scaled);

  //! Constructor for absolute norm
  /*! This constructor defaults to the \c Absolute ToleranceType and \c TWO NormType. */
  NormUpdate(double tolerance, ScaleType stype = Scaled);

  //! Destructor.
  virtual ~NormUpdate();

  virtual NOX::StatusTest::StatusType
  checkStatus(const NOX::Solver::Generic& problem,
          NOX::StatusTest::CheckType checkType);

  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual std::ostream& print(std::ostream& stream, int indent = 0) const;

  /* @name Accessor Functions
     Used to query current values of variables in the status test.
  */
  //@{

  //! Returns the value of the Update-norm computed in the last call to checkStatus.
  virtual double getNormUpdate() const;

  //! Returns the true tolerance.
  virtual double getTolerance() const;

  //@}

private:

  //! %Status
  NOX::StatusTest::StatusType status;

  //! Vector containing the update for the current outer iteration
  Teuchos::RCP<NOX::Abstract::Vector> updateVectorPtr;

  //! Type of norm to use
  NOX::Abstract::Vector::NormType normType;

  //! Scaling to use
  ScaleType scaleType;

  //! Tolerance value specified by user in constructor
  double tolerance;

  //! Norm of the update to be compared to trueTolerance
  double normUpdate;

};

} // namespace Status
} // namespace NOX

#endif
