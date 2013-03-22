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

#include "NOX_TestCompare.H"
#include "Teuchos_RCP.hpp"

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
  double testTolerance;

  // test for zero
  if (compareType == NOX::TestCompare::Relative && ((value_expected == 0.0) ||
						    (value == 0.0)))
    compareType = NOX::TestCompare::Absolute;

  testValue = fabs(value-value_expected);
  if (compareType == NOX::TestCompare::Absolute)
    testTolerance = tolerance;
  else
    testTolerance = tolerance * fabs(value_expected);

  if (testValue > testTolerance)
    passed = false;
  else
    passed = true;

  if (utils.isPrintType(NOX::Utils::TestDetails)) {
    os << std::endl
	 << "\tChecking " << name << ":  ";
    if (passed)
      os << "Passed." << std::endl;
    else
      os << "Failed." << std::endl;
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
  double inf_norm = computeVectorNorm(vec, vec_expected, rtol, atol);

  if (inf_norm < 1)
    passed = true;
  else
    passed = false;

  if (utils.isPrintType(NOX::Utils::TestDetails)) {
    os << std::endl
	 << "\tChecking " << name << ":  ";
    if (passed)
      os << "Passed." << std::endl;
    else
      os << "Failed." << std::endl;
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

int 
NOX::TestCompare::testMatrix(
		 const NOX::Abstract::MultiVector::DenseMatrix& mat, 
		 const NOX::Abstract::MultiVector::DenseMatrix& mat_expected, 
		 double rtol, double atol, 
		 const std::string& name)
{
  bool passed;

  NOX::Abstract::MultiVector::DenseMatrix tmp(mat_expected.numRows(),
					      mat_expected.numCols());

  for (int j=0; j<mat_expected.numCols(); j++)
    for (int i=0; i<mat_expected.numRows(); i++)
      tmp(i,j) = fabs(mat(i,j)-mat_expected(i,j)) / 
	(atol + rtol * fabs(mat_expected(i,j)));
 
  double inf_norm = tmp.normInf();

  if (inf_norm < 1)
    passed = true;
  else
    passed = false;

  if (utils.isPrintType(NOX::Utils::TestDetails)) {
    os << std::endl
	 << "\tChecking " << name << ":  ";
    if (passed)
      os << "Passed." << std::endl;
    else
      os << "Failed." << std::endl;
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

int 
NOX::TestCompare::testMultiVector(
			     const NOX::Abstract::MultiVector& mvec, 
			     const NOX::Abstract::MultiVector& mvec_expected, 
			     double rtol, double atol, 
			     const std::string& name)
{
  bool passed;
  double inf_norm;
  double inf_norm_max = 0.0;

  for (int i=0; i<mvec_expected.numVectors(); i++) {
    inf_norm = computeVectorNorm(mvec[i], mvec_expected[i], rtol, atol);
    if (inf_norm > inf_norm_max)
      inf_norm_max = inf_norm;
  }

  if (inf_norm_max < 1)
    passed = true;
  else
    passed = false;

  if (utils.isPrintType(NOX::Utils::TestDetails)) {
    os << std::endl
	 << "\tChecking " << name << ":  ";
    if (passed)
      os << "Passed." << std::endl;
    else
      os << "Failed." << std::endl;
    os << "\t\tComputed norm:        " << utils.sciformat(inf_norm_max) 
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

double NOX::TestCompare::computeVectorNorm(
				   const NOX::Abstract::Vector& vec, 
				   const NOX::Abstract::Vector& vec_expected, 
				   double rtol, double atol)
{
  // Compute atol + rtol*|vec_expected|
  Teuchos::RCP<NOX::Abstract::Vector> tmp1 = vec.clone(NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::Vector> tmp2 = vec.clone(NOX::ShapeCopy);
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

  return inf_norm;
}
