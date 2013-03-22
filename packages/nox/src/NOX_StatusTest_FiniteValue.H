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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
