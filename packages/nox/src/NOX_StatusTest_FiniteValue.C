// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
#include "NOX_Common.H"  // for string class
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

ostream& NOX::StatusTest::FiniteValue::print(ostream& stream, int indent) const
{


  // Set the correct label for the check result
  string label = "Unknown";
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
  stream << endl;
  
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
#ifdef HAVE_NAN_SUPPORT
  if (isnan(x))
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
#else
  // Use IEEE 754 definition: Inf * 0 = NaN
  double z = 0.0 * x;
  if (NOX_isnan(z)) 
#endif 
    return true;

  return false;
}
