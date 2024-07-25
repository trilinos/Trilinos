// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
\brief  Unit test (CubatureDirect): correctness of
        integration of monomials for 1D reference cells.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid2_CubatureTensorSorted.hpp"
#include "Intrepid2_AdaptiveSparseGrid.hpp"
#include "Intrepid2_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace Intrepid2;

/*
  Computes integrals of monomials over a given reference cell.
*/
long double evalQuad(std::vector<int> power, int dimension, int order, 
		     std::vector<EIntrepidBurkardt> rule,
		     std::vector<EIntrepidGrowth> growth) {
 
  CubatureTensorSorted<long double> lineCub(0,dimension);
  AdaptiveSparseGrid<long double,std::vector<long double> >::buildSparseGrid(
						   lineCub,dimension,
						   order,rule,
						   growth,false);

  int size = lineCub.getNumPoints();
  FieldContainer<long double> cubPoints(size,dimension);
  FieldContainer<long double> cubWeights(size);
  lineCub.getCubature(cubPoints,cubWeights);

  int mid  = size/2;
  long double Q = 0.0;
  if (size%2) {
    Q  = cubWeights(mid);
    for (int i=0; i<dimension; i++) {
      Q *= powl(cubPoints(mid,i),(long double)power[i]);
    }
  }

  for (int i=0; i<mid; i++) {
    long double value1 = cubWeights(i);
    long double value2 = cubWeights(size-i-1);
    for (int j=0; j<dimension; j++) {
      value1 *= powl(cubPoints(i,j),(long double)power[j]);
      value2 *= powl(cubPoints(size-i-1,j),(long double)power[j]);
    }
    Q += value1+value2;
  }
  return Q;
}

long double evalInt(int dimension, std::vector<int> power, 
		    std::vector<EIntrepidBurkardt> rule) {
  long double I = 1.0;

  for (int i=0; i<dimension; i++) {
    if (rule[i]==BURK_CLENSHAWCURTIS||rule[i]==BURK_FEJER2||
	rule[i]==BURK_LEGENDRE||rule[i]==BURK_PATTERSON || 
	rule[i]==BURK_TRAPEZOIDAL) {
      if (power[i]%2)
	I *= 0.0;
      else 
	I *= 2.0/((long double)power[i]+1.0);
    }
    else if (rule[i]==BURK_LAGUERRE) {
      I *= tgammal((long double)(power[i]+1));
    }
    else if (rule[i]==BURK_CHEBYSHEV1) {
      long double bot, top;
      if (!(power[i]%2)) {
	top = 1; bot = 1;
	for (int j=2;j<=power[i];j+=2) {
	  top *= (long double)(j-1);
	  bot *= (long double)j;
	}
	I *= M_PI*top/bot;
      }
      else {
	I *= 0.0;
      }
    }
    else if (rule[i]==BURK_CHEBYSHEV2) {
      long double bot, top;
      if (!(power[i]%2)) {
      top = 1; bot = 1;
      for (int j=2;j<=power[i];j+=2) {
	top *= (long double)(j-1);
	bot *= (long double)j;
      }
      bot *= (long double)(power[i]+2);
      I   *= M_PI*top/bot;
      }
      else {
	I *= 0.0;
      }
    }
    else if (rule[i]==BURK_HERMITE||rule[i]==BURK_GENZKEISTER) {
      if (power[i]%2) {
	I *= 0.0;
      }
      else {  
	long double value = 1.0;
	if ((power[i]-1)>=1) {
	  int n_copy = power[i]-1;
	  while (1<n_copy) {
	    value  *= (long double)n_copy;
	    n_copy -= 2;
	  }
	}
	I *= value*sqrt(M_PI)/powl(2.0,(long double)power[i]/2.0);
      }
    }
  }
  return I;
}

int main(int argc, char *argv[]) {
  
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
Kokkos::initialize();
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
  << "|                         Unit Test (CubatureTensorSorted)                    |\n" \
  << "|                                                                             |\n" \
  << "|     1) Computing integrals of monomials in 2D                               |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Drew Kouri (dpkouri@sandia.gov) or                     |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 22: integrals of monomials in 2D - Anisotropic but no growth rules     |\n"\
  << "===============================================================================\n";


  // internal variables:
  int dimension           = 2;
  int         errorFlag   = 0;
  long double reltol      = 1.0e+02*INTREPID_TOL;
  int         maxDeg      = 0;    
  long double analyticInt = 0;
  long double testInt     = 0;
  int         maxOrder    = 6;
  std::vector<int> power(2,0);
  std::vector<EIntrepidBurkardt> rule1(2,BURK_CLENSHAWCURTIS);
  std::vector<EIntrepidGrowth> growth1(2,GROWTH_FULLEXP);

  *outStream << "\nIntegrals of monomials on a reference line (edge):\n";
  // compute and compare integrals
  try {
    for (int i=0; i<=maxOrder; i++) {
      maxDeg = i-1;
      for (int j=0; j <= maxDeg; j++) {
	power[0] = j;
	for (int k=0; k <= maxDeg; k++) {
	  power[1] = k;
	  if (j+k < maxDeg) {
	    analyticInt = evalInt(dimension, power, rule1);
	    testInt     = evalQuad(power,dimension,maxOrder,rule1,growth1);
	    
	    long double abstol  = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	    long double absdiff = std::fabs(analyticInt - testInt);
	    *outStream << "Cubature order " << std::setw(2) << std::left << i 
		       << " integrating " << "x^" << std::setw(2) 
		       << std::left << j << "y^" << std::setw(2) << std::left 
		       << k <<  ":" << "   " << std::scientific 
		       << std::setprecision(16) << testInt 
		       << "   " << analyticInt << "   " 
		       << std::setprecision(4) << absdiff << "   " 
		       << "<?" << "   " << abstol << "\n";
	    if (absdiff > abstol) {
	      errorFlag++;
	      *outStream << std::right << std::setw(104) << "^^^^---FAILURE!\n";
	    }
	  }
	} // end for k
	*outStream << "\n";
      } // end for j
      *outStream << "\n";
    } // end for i
  }
  catch (std::logic_error &err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  };


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
Kokkos::finalize();
  return errorFlag;
}
