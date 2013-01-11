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

#include "NOX_StatusTest_FiniteValue.H" // class definition
#include "NOX_Common.H"  // for std::string class
#include "NOX_Solver_Generic.H"  
#include "NOX_Abstract_Group.H"

NOX::StatusTest::FiniteValue::
FiniteValue(VectorType v, NOX::Abstract::Vector::NormType n) :
  vectorType(v),
  vectorTypeLabel("?"),
  normType(n),
  normTypeLabel("?"),
  status(Unevaluated),
  result(-1),
  normValue(-1.0)
{
  // Set the vector type label for printing
  if (vectorType == FVector) 
    vectorTypeLabel = "F";
  else
    vectorTypeLabel = "Solution";

  // Set the norm type label for printing
  if (normType == NOX::Abstract::Vector::TwoNorm)
    normTypeLabel = "Two-Norm";
  else if (normType == NOX::Abstract::Vector::OneNorm)
    normTypeLabel = "One-Norm";
  else
    normTypeLabel = "Max-Norm";

}

NOX::StatusTest::FiniteValue::~FiniteValue()
{
}

NOX::StatusTest::StatusType NOX::StatusTest::FiniteValue::
checkStatus(const Solver::Generic& problem, 
	    NOX::StatusTest::CheckType checkType)
{
  // Reset the check
  normValue = -1.0;
  const NOX::Abstract::Group& grp = problem.getSolutionGroup();

  switch (checkType)
  {
  case NOX::StatusTest::Complete:
  case NOX::StatusTest::Minimal:

  
    if (vectorType == FVector) 
    {
      if (normType == NOX::Abstract::Vector::TwoNorm)
	normValue = grp.getNormF();  // More efficient than recomputing norm
      else
	normValue = grp.getF().norm(normType);
    }
    else
      normValue = grp.getX().norm(normType);
    
    result = finiteNumberTest(normValue);

    status = (result == 0) ? Unconverged : Failed;
    break;

  case NOX::StatusTest::None:
  default:
    result = 1;
    status = Unevaluated;
    break;
  }
  
  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::FiniteValue::getStatus() const
{
  return status;
}

std::ostream& NOX::StatusTest::FiniteValue::print(std::ostream& stream, int indent) const
{


  // Set the correct label for the check result
  std::string label = "Unknown";
  if (result == 0)
    label = "Finite";
  else if (result == -1)
    label = "NaN";
  else if (result == -2) 
    label = "Infinite";

  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Finite Number Check (" << normTypeLabel;
  stream << " " << vectorTypeLabel;
  stream << ") = ";
  stream << label;
  //stream << " (" << normValue << ")";
  stream << std::endl;
  
  return stream;
}

int NOX::StatusTest::FiniteValue::finiteNumberTest(double x) const
{
  if (NOX_isnan(x))
    return -1;
 
  if (NOX_isinf(x))
    return -2;

  return 0;
}

bool NOX::StatusTest::FiniteValue::NOX_isnan(double x) const
{
#if defined(HAVE_NAN_SUPPORT)
  if (isnan(x))
#elif defined(FINITE_VALUE_HAVE_GLOBAL_ISNAN)
  if (isnan(x))
#elif defined(FINITE_VALUE_HAVE_STD_ISNAN)
  if (std::isnan(x))
#else
  if (x != x) 
#endif
    return true;

  return false;
}

bool NOX::StatusTest::FiniteValue::NOX_isinf(double x) const
{
#ifdef HAVE_INF_SUPPORT
  if (isinf(x))
#elif defined(FINITE_VALUE_HAVE_GLOBAL_ISINF)
  if (isinf(x))
#elif defined(FINITE_VALUE_HAVE_STD_ISINF)
  if (std::isinf(x))
#else
  // Use IEEE 754 definition: Inf * 0 = NaN
  double z = 0.0 * x;
  if (NOX_isnan(z)) 
#endif 
    return true;

  return false;
}
