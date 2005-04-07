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

#include "NOX_TestCompare.H"

NOX::TestCompare::TestCompare(std::ostream& outputStream, 
			      const NOX::Utils& utilities) :
  os(outputStream),
  utils(utilities) 
{
}

NOX::TestCompare::~TestCompare()
{
}

int NOX::TestCompare::testValue(double value, 
				double value_expected, 
				double tolerance, 
				const std::string& name, 
				NOX::TestCompare::CompareType compareType)
{
  bool passed;
  double testValue;

  // test for zero
  if (compareType == NOX::TestCompare::Relative && value_expected == 0.0)
    compareType = NOX::TestCompare::Absolute;

  if (compareType == NOX::TestCompare::Absolute)
    testValue = fabs(value-value_expected);
  else
    testValue = fabs(value-value_expected) / fabs(value_expected);

  if (testValue > tolerance)
    passed = false;
  else
    passed = true;

  if (utils.isPrintProcessAndType(NOX::Utils::TestDetails)) {
    os << std::endl
	 << "\tChecking " << name << ":  ";
    if (passed)
      os << "Passed." << std::endl;
    else
      os << "Failed." <<std:: endl;
    os << "\t\tExpected value:       " << utils.sciformat(value_expected) 
       << std::endl
       << "\t\tComputed value:       " << utils.sciformat(value) 
       << std::endl
       << "\t\tTolerance:            " << utils.sciformat(tolerance) 
       << std::endl;
    if (compareType == NOX::TestCompare::Absolute)
      os << "\t\tAbsolute ";
    else
      os << "\t\tRelative ";
    os << "Difference:  " 
       << utils.sciformat(testValue) << std::endl;
  }

  if (passed)
    return 0;
  else
    return 1;
}

int NOX::TestCompare::testVector(const NOX::Abstract::Vector& vec, 
				 const NOX::Abstract::Vector& vec_expected, 
				 double rtol, double atol, 
				 const std::string& name)
{
  bool passed;
  double testValue;

  // Compute atol + rtol*|vec_expected|
  NOX::Abstract::Vector *tmp1 = vec.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *tmp2 = vec.clone(NOX::ShapeCopy);
  tmp1->init(atol);
  tmp2->abs(vec_expected);
  tmp1->update(rtol, *tmp2, 1.0);

  // Compute 1/(atol + rtol*|vec_expected|)
  tmp2->reciprocal(*tmp1);

  // Compute |vec - vec_expected|
  tmp1->update(1.0, vec, -1.0, vec_expected, 0.0);
  tmp1->abs(*tmp1);

  // Compute |vec - vec_expected|/(atol + rtol*|vec_expected|)
  tmp1->scale(*tmp2);

  double inf_norm = tmp1->norm(NOX::Abstract::Vector::MaxNorm);

  if (inf_norm < 1)
    passed = true;
  else
    passed = false;

  delete tmp1;
  delete tmp2;

  if (utils.isPrintProcessAndType(NOX::Utils::TestDetails)) {
    os << std::endl
	 << "\tChecking " << name << ":  ";
    if (passed)
      os << "Passed." << std::endl;
    else
      os << "Failed." <<std:: endl;
    os << "\t\tComputed norm:        " << utils.sciformat(inf_norm) 
       << std::endl
       << "\t\tRelative Tolerance:   " << utils.sciformat(rtol) 
       << std::endl
       << "\t\tAbsolute Tolerance:   " << utils.sciformat(rtol) 
       << std::endl;
  }

  if (passed)
    return 0;
  else
    return 1;
}

