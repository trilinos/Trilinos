// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "DTaylorUnitTests.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION( DTaylorOpsUnitTest );

DTaylorOpsUnitTest::DTaylorOpsUnitTest() :
  urand(0.0, 1.0), d(5), tol_a(1.0e-14), tol_r(1.0e-14) 
{
  X = new double*[2];
  X[0] = new double[d+1];
  X[1] = new double[d+1];

  Y = new double*[1];
  Y[0] = new double[d+1];
}

DTaylorOpsUnitTest::DTaylorOpsUnitTest(unsigned int degree, 
				       double absolute_tolerance, 
				       double relative_tolerance) :
  urand(0.0, 1.0), 
  d(degree), 
  tol_a(absolute_tolerance), 
  tol_r(relative_tolerance) 
{
  X = new double*[2];
  X[0] = new double[d+1];
  X[1] = new double[d+1];

  Y = new double*[1];
  Y[0] = new double[d+1];
}

DTaylorOpsUnitTest::~DTaylorOpsUnitTest()
{
  delete [] X[1];
  delete [] X[0];
  delete [] X;

  delete [] Y[0];
  delete [] Y;
}

void DTaylorOpsUnitTest::setUp() {
  double val;

  a_dtay = DTaylorType(d,0.0);
  b_dtay = DTaylorType(d,0.0);
  
  for (unsigned int i=0; i<=d; i++) {
    val = urand.number();
    a_dtay.fastAccessCoeff(i) = val;
    X[0][i] = val;

    val = urand.number();
    b_dtay.fastAccessCoeff(i) = val;
    X[1][i] = val;

    Y[0][i] = 0.0;
  }
}

void DTaylorOpsUnitTest::tearDown() {}

void DTaylorOpsUnitTest::comparePolys(const DTaylorType& x_dtay,
					double* x_adolc) {

  // Compare degrees
  CPPUNIT_ASSERT(x_dtay.degree() == d);

//   std::cout << std::endl << "DTaylor:" << x_dtay << std::endl;
//   std::cout << "ADOLC:  ";
//   print_poly(x_adolc);
//   std::cout << ":Diff:  ";
//   print_diff(x_dtay, x_adolc);

  
  // Compare coefficients
  for (unsigned int i=0; i<=d; i++) {
    compareDoubles(x_dtay.coeff(i), x_adolc[i]);
    //compareDoubles(x_dtay.fastAccessCoeff(i), x_adolc[i]);
  }
}

void DTaylorOpsUnitTest::compareDoubles(double a, double b) {
  //cout << fabs(a-b) << "   " << tol_a + tol_r*fabs(a) << endl;
  CPPUNIT_ASSERT( fabs(a-b) < tol_a + tol_r*fabs(a) );
}

void DTaylorOpsUnitTest::print_poly(double *x) {
  std::cout.setf(std::ios::fixed,std::ios::floatfield);
  std::cout.width(12);
  std::cout << "[";
      
  for (unsigned int i=0; i<=d; i++) {
    std::cout.width(12);
    std::cout << x[i];
  }

  std::cout << "]\n";
}

void DTaylorOpsUnitTest::print_diff(const DTaylorType& x_dtay,
				    double *x) {
  std::cout.setf(std::ios::scientific,std::ios::floatfield);
  //std::cout.width(12);
  std::cout << "[";
      
  for (unsigned int i=0; i<=d; i++) {
    //std::cout.width(12);
    std::cout << x_dtay.coeff(i) - x[i] << " ";
  }

  std::cout << "]\n";
}
