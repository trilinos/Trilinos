// $Id$ 
// $Source$ 

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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

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
  checkStatus(const NOX::Solver::Generic& problem);

  virtual NOX::StatusTest::StatusType 
  checkStatusEfficiently(const NOX::Solver::Generic& problem, 
			 NOX::StatusTest::CheckType checkType);

  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual ostream& print(ostream& stream, int indent = 0) const;

  //! The finite number test algorithm.  
  /*! We made this method public so that other objects (solvers, line
    searches, and directions) can use this test on their own values.

    Return Values:
    - 0 = Finite Number
    - -1 = NaN
    - -2 = Inf
  */
  virtual int finiteNumberTest(double x) const;

protected:

private:

  //! Determines what type of vector (F or solution) to use in computing the norm.
  VectorType vectorType;

  //! Vector label used for printing.
  string vectorTypeLabel;
  
  //! Determines what norm type to use (NOX does not have access to vector entries so it must use the norm to evaluate NaN/Inf problems).
  NOX::Abstract::Vector::NormType normType;

  //! Norm type label for printing.
  string normTypeLabel;

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
