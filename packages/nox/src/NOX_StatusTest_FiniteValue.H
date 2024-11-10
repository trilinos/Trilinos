// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_STATUS_FINITEVALUE_H
#define NOX_STATUS_FINITEVALUE_H

#include "NOX_StatusTest_Generic.H"  // base class
#include "NOX_Abstract_Vector.H"     // for NOX::Abstract::Vector::NormType

namespace NOX {

namespace StatusTest {

//! Failure test based on whether the norm of a vector has a finite value.
/*!

  This test returns NOX::StatusTest::Failed if the norm of a vector is
  calssified as a NaN or Inf.  Otherwise, it returns
  NOX::StatusTest::Unconverged.  The user can specify whether to use
  the F vector or the solution vector from the current solution group
  for the check.  NOX does not have access to vector entries so the
  number used in the NaN/Inf check is based on the norm of a vector.

  \note If checkStatus is called with the type set to
  NOX::StatusTest::None, then the status is set to
  NOX::Status::Unevaluated and returned.
*/
class FiniteValue : public Generic {

public:

  //! Type of vector to run the finite value check on.
  enum VectorType {FVector, SolutionVector};

  //! Constructor. Specify which vector to check and with what norm to use.
  FiniteValue(VectorType v = FVector,
       NOX::Abstract::Vector::NormType n = NOX::Abstract::Vector::TwoNorm);

  //! Destructor.
  virtual ~FiniteValue();

  virtual NOX::StatusTest::StatusType
  checkStatus(const NOX::Solver::Generic& problem,
          NOX::StatusTest::CheckType checkType);

  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual std::ostream& print(std::ostream& stream, int indent = 0) const;

  /*! \brief The finite number test algorithm.

    Autoconf will test to see if the compiler implements the isnan()
    and isinf() functions in the cmath or math.h headers.  If so, we
    will use these.  If not, we supply a default implementation. The
    default implementation is only guaranteed to work if the code is
    IEEE 748/754 compliant.  The checks for isnan and isinf are
    separate because compilers like the old sgi platforms support one
    but not the other.  See bug 2019 for more details.

    This method is public so that other objects (solvers, line
    searches, and directions) can use this test on their own values.

    Return Values:
    - 0 = Finite Number
    - -1 = NaN
    - -2 = Inf

  */
  virtual int finiteNumberTest(double x) const;

protected:

  //! NOX's implementation of isnan().
  bool NOX_isnan(double x) const;

  //! NOX's implementation of isinf().
  bool NOX_isinf(double x) const;

private:

  //! Determines what type of vector (F or solution) to use in computing the norm.
  VectorType vectorType;

  //! Vector label used for printing.
  std::string vectorTypeLabel;

  //! Determines what norm type to use (NOX does not have access to vector entries so it must use the norm to evaluate NaN/Inf problems).
  NOX::Abstract::Vector::NormType normType;

  //! Norm type label for printing.
  std::string normTypeLabel;

  //! Current %Status
  NOX::StatusTest::StatusType status;

  //! Return code from the call:  result = fpclassify(double normValue).
  int result;

  //! Computed norm value.
  double normValue;

};

} // namespace StatusTest
} // namespace NOX

#endif
