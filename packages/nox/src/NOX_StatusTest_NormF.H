// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_STATUSTEST_NORMF_H
#define NOX_STATUSTEST_NORMF_H

#include "NOX_StatusTest_Generic.H"    // base class
#include "NOX_Abstract_Vector.H"        // for NormType
#include "NOX_Utils.H"                  // for std::cerr

// Forward declaration
namespace NOX {
namespace Abstract {
class Group;
}
}

namespace NOX {

namespace StatusTest {

//! Various convergence tests based on the norm of the residual.
/*!

  Use the constructor to define the test based on the type
  of scaling (see ScaleType) and the type of Tolerance (see
  Tolerance).

  If checkStatus is called with the type set to NOX::StatusTest::None,
  then the status is set to NOX::StatusTest::Unevaluated and returned.
  (Also #normF is set to 0.0.)

  If checkStatus is called on a problem where the solution group does
  not have F evaluated (i.e., problem.getSolutionGroup().isF() is
  false), then the status is set to NOX::StatusTest::Unconverged and
  returned.   (Also #normF is set to -1.0.)

  Finally, we return NOX::StatusTest::Converged if \f$\alpha <
  \beta\f$, and NOX::StatusTest::Unconverged otherwise.
  Here \f$\alpha\f$
  represents the norm of \f$ F(x)\f$ and \f$\beta\f$ represents the
  tolerance, as described below.

  Let \f$\gamma\f$ denote an optional scale factor defined as

  - \f$\gamma = \frac{1}{n}\f$ if \c sType in the constructor is NOX::NormF::Scaled, and


  Then \f$\alpha\f$ is defined as follows:

  - If \c nType in the constructor is Abstract::Vector::TWO, then
    \f[ \alpha = \sqrt{ \gamma \sum_{i=1}^n F_i^2 } \f]

  - If \c nType in the constructor is Abstract::Vector::ONE, then
    \f[ \alpha = \gamma \sum_{i=1}^n | F_i | \f]

  - If \c nType in the constructor is Abstract::Vector::INF, then
    \f[ \alpha = \gamma \max_{i} | F_i |  \f]

  We set \f$\beta\f$ as follows, based on the value of \c tolerance in
  the constructor.

  - If an initial guess is provided, we use a relative tolerance
  defined by \f[ \beta = \alpha_0 * \mbox{tolerance} \f] Here
  \f$\alpha_0\f$ is the \f$\alpha\f$ (as defined above) associated with
  the initial guess.

  - Otherwise, we use an absolute tolerance defined by
  \f[  \beta  = \mbox{tolerance} \f]

*/
class NormF : public Generic {

public:

  //! Type that determines whether to scale the norm by the problem size.
  enum ScaleType {
    //! No norm scaling
    Unscaled,
    //! Scale the norm by the length of the vector
    Scaled
  };

  //! Type that determines whether the norm is absolute or relative to the intial guess
  enum ToleranceType {
    //! Relative to starting guess
    Relative,
    //! Absolute
    Absolute
  };

  //! Constructor for absolute norm.
  /*! This constructor defaults to the \c Absolute tolerance type. */
  NormF(double tolerance, NOX::Abstract::Vector::NormType ntype,
    ScaleType stype = Scaled, const NOX::Utils* u = NULL);

  //! Constructor for absolute norm
  /*! This constructor defaults to the \c Absolute ToleranceType and \c TWO NormType. */
  NormF(double tolerance, ScaleType stype = Scaled,
    const NOX::Utils* u = NULL);

  //! Constructor with initial guess (for relative norms)
  /*! This constructor defaults to the \c Relative tolerance type. */
  NormF(NOX::Abstract::Group& initialGuess, double tolerance,
    NOX::Abstract::Vector::NormType ntype, ScaleType stype = Scaled,
    const NOX::Utils* u = NULL);

  //! Constructor with initial guess (for relative norms)
  /*! This constructor defaults to the \c Relative ToleranceType and \c TWO NormType. */
  NormF(NOX::Abstract::Group& initialGuess, double tolerance,
    ScaleType stype = Scaled, const NOX::Utils* u = NULL);

  //! Destructor.
  virtual ~NormF();

  // derived
  virtual NOX::StatusTest::StatusType
  checkStatus(const NOX::Solver::Generic& problem,
          NOX::StatusTest::CheckType checkType);

  // derived
  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual std::ostream& print(std::ostream& stream, int indent = 0) const;

  /* @name Reset Functions
     Used to change the tolerances in the status test after construction.
  */
  //@{

  //! Resets the user specified absolute or relative tolerance.
  virtual void reset(double tolerance);

  //! Resets the user specified relative tolerance.
  virtual void reset(NOX::Abstract::Group& initialGuess, double tolerance);

  //@}

  /* @name Accessor Functions
     Used to query current values of variables in the status test.
  */
  //@{

  //! Returns the value of the F-norm computed in the last call to checkStatus.
  virtual double getNormF() const;

  //! Returns the true tolerance.
  virtual double getTrueTolerance() const;

  //! Returns the specified tolerance set in the constructor.
  virtual double getSpecifiedTolerance() const;

  //! Returns the initial tolerance.
  virtual double getInitialTolerance() const;

  //@}


private:

  /*! \brief Calculate the norm of F for the given group according to
    the scaling type, norm type, and tolerance type.

    \note Returns -1.0 if F(x) has not been calculated for the given
    grp (i.e., grp.isF() is false).
  */
  double computeNorm(const NOX::Abstract::Group& grp);

  /*! In the case of a relative norm calculation, initializes
    \c trueTolerance based on the F-value at the initial guess.*/
  void relativeSetup(NOX::Abstract::Group& initialGuess);

private:

  //! %Status
  NOX::StatusTest::StatusType status;

  //! Type of norm to use
  NOX::Abstract::Vector::NormType normType;

  //! Scaling to use
  ScaleType scaleType;

  //! Tolerance type (i.e., relative or absolute)
  ToleranceType toleranceType;

  //! Tolerance required for convergence.
  double specifiedTolerance;

  //! Initial tolerance
  double initialTolerance;

  //! True tolerance value, i.e., specifiedTolerance / initialTolerance
  double trueTolerance;

  //! Norm of F to be compared to trueTolerance
  double normF;

  //! Ostream used to print errors
  NOX::Utils utils;

};

} // namespace Status
} // namespace NOX

#endif
