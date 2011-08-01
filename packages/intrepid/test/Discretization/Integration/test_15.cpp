// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER


/** \file
\brief  Unit test (CubatureDirect): correctness of
        integration of monomials for 1D reference cells.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_CubatureLineSorted.hpp"
#include "Intrepid_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace Intrepid;

/*
  Computes integrals of monomials over a given reference cell.
*/
long double evalQuad(int order, int power, EIntrepidBurkardt rule) {
 
  CubatureLineSorted<long double> lineCub(rule,order,false);
  CubatureLineSorted<long double> lineCub2(rule,order-1,false);
  lineCub.update(-1.0,lineCub2,1.0);
  int size = lineCub.getNumPoints();
  FieldContainer<long double> cubPoints(size);
  FieldContainer<long double> cubWeights(size);
  lineCub.getCubature(cubPoints,cubWeights);

  long double Q = 0.0;
  for (int i=0; i<size; i++) {
    Q += cubWeights(i)*powl(cubPoints(i),power);
  }
  return Q;

  /*
  int mid  = size/2;
  long double Q = 0.0;
  if (size%2) 
    Q = cubWeights(mid)*powl(cubPoints(mid),power);

  for (int i=0; i<mid; i++) {
    Q += cubWeights(i)*powl(cubPoints(i),power)+
      cubWeights(order-i-1)*powl(cubPoints(size-i-1),power);
  }
  return Q;
  */
}

int main(int argc, char *argv[]) {
  
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
 
  *outStream \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                         Unit Test (CubatureLineSorted)                      |\n" \
  << "|                                                                             |\n" \
  << "|     1) Computing differential integrals of monomials in 1D                  |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Drew Kouri (dpkouri@sandia.gov) or                     |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 15: differential integrals of monomials in 1D                          |\n"\
  << "===============================================================================\n";

  // internal variables:
  int         errorFlag   = 0;
  long double abstol      = 1.0e+05*INTREPID_TOL;
  int         maxDeg      = 0;    
  long double analyticInt = 0;
  long double testInt     = 0;
  int         maxOrder    = 60;
  EIntrepidBurkardt rule  = BURK_LEGENDRE;

  *outStream << "\nDifferential integrals of monomials on a reference line:\n";
  // compute and compare integrals
  try {
    *outStream << "Testing " << EIntrepidBurkardtToString(rule) << "\n";
    // compute integrals
    for (int i=2; i <= maxOrder; i++) {
      maxDeg = 2*(i-1)-1;  
      for (int j=0; j <= maxDeg; j++) {
	testInt             = evalQuad(i,j,rule);
	long double absdiff = std::fabs(testInt);
	*outStream << "Cubature order " << std::setw(2) << std::left << i << " integrating "
		   << "x^" << std::setw(2) << std::left << j <<  ":" << "   "
		   << std::scientific << std::setprecision(16) << testInt << "   " << analyticInt 
		   << "   " << std::setprecision(4) << absdiff << "   " << "<?" << "   " << abstol 
		   << "\n";
	if (absdiff > abstol) {
	  errorFlag++;
	  *outStream << std::right << std::setw(104) << "^^^^---FAILURE!\n";
	}
      } // end for j
      *outStream << "\n";
    } // end for i
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  };


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
