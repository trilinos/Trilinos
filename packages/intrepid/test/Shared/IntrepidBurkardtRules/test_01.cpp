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

/** \file   test_01.cpp
    \brief  Test file for integration rules provided by John Burkardt.
	    <A HREF="http://people.sc.fsu.edu/~jburkardt/cpp_src/sandia_rules/sandia_rules.html">
	    <\A>
    \author Created by D. Kouri and D. Ridzal.
 */

#include "Intrepid_BurkardtRules.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Array.hpp"

//# include <cstdlib>
//# include <iostream>
//# include <cmath>
//# include <iomanip>

using namespace std;
using namespace Intrepid;

#define INTREPID_TEST_COMMAND( S )                                                                                  \
{                                                                                                                   \
  try {                                                                                                             \
    S ;                                                                                                             \
  }                                                                                                                 \
  catch (std::logic_error err) {                                                                                    \
      *outStream << "Expected Error ----------------------------------------------------------------\n";            \
      *outStream << err.what() << '\n';                                                                             \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";    \
  };                                                                                                                \
}

template<class Scalar>
Scalar evalQuad(int order, int power, Scalar x[], Scalar w[]) {
 
  int mid  = order/2;
  Scalar Q = 0.0;
  if (order%2) 
    Q = w[mid]*powl(x[mid],power);

  for (int i=0; i<mid; i++) {
    Q += w[i]*powl(x[i],power)+w[order-i-1]*powl(x[order-i-1],power);
  }

  return Q;
  /* 
  Scalar Q = 0.0;
  for (int i=0; i<order; i++) {
    Q += w[i]*powl(x[i],power);
  }
  return Q;
  */
}

template<class Scalar>
Scalar factorial2 (int n) {
  Scalar value = 1.0;
  if (n<1)
    return value;

  int n_copy = n;
  while (1<n_copy) {
    value  *= (Scalar)n_copy;
    n_copy -= 2;
  }
  return value;
}

template<class Scalar>
Scalar chebyshev1(int power) {
  Scalar bot, exact, top;
  if (!(power%2)) {
    top = 1; bot = 1;
    for (int i=2;i<=power;i+=2) {
      top *= (Scalar)(i-1);
      bot *= (Scalar)i;
    }
    exact = M_PI*top/bot;
  }
  else {
    exact = 0.0;
  }
  return exact;
}

template<class Scalar>
Scalar chebyshev2(int power) {
  Scalar bot, exact, top;
  if (!(power%2)) {
    top = 1; bot = 1;
    for (int i=2;i<=power;i+=2) {
      top *= (Scalar)(i-1);
      bot *= (Scalar)i;
    }
    bot *= (Scalar)(power+2);
    exact = M_PI*top/bot;
  }
  else {
    exact = 0.0;
  }
  return exact;
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
  << "|                       Unit Test (IntrepidBurkardtRules)                     |\n" \
  << "|                                                                             |\n" \
  << "|     1) the Burkardt rule tests                                              |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Drew Kouri (dpkouri@sandia.gov) or                     |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n";


  int errorFlag = 0;

  int maxOrder  = 30;
  long double reltol = 1e-8;
  long double analyticInt = 0, testInt = 0;
  // compute and compare integrals
  try {
    
    *outStream << "Gauss-Legendre Cubature \n";
    *outStream << "Integrates functions on [-1,1] weighted by w(x) = 1\n";
    for (int i = 1; i<=maxOrder; i++) {
      Teuchos::Array<long double> nodes(i), weights(i);
      IntrepidBurkardtRules::legendre_compute(i,nodes.getRawPtr(),weights.getRawPtr());
      for (int j=0; j<=2*i-1; j++) {
	if (j%2)
	  analyticInt = 0.0;
	else 
	  analyticInt = 2.0/((long double)j+1.0);
	testInt = evalQuad(i,j,nodes.getRawPtr(),weights.getRawPtr());
	long double abstol = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	long double absdiff = std::fabs(analyticInt - testInt);
	*outStream << "Cubature order " << std::setw(2) << std::left << i << " integrating "
		   << "x^" << std::setw(2) << std::left << j << ":" << "   "
		   << std::scientific << std::setprecision(16) << testInt << "   " 
		   << analyticInt << "   " << std::setprecision(4) << absdiff << "   " << "<?" 
		   << "   " << abstol << "\n";
	if (absdiff > abstol) {
	  errorFlag++;
	  *outStream << std::right << std::setw(111) << "^^^^---FAILURE!\n";
	}
      }
    }
    *outStream << "\n";    

    *outStream << "Clenshaw-Curtis Cubature \n";
    *outStream << "Integrates functions on [-1,1] weighted by w(x) = 1\n";
    for (int i = 1; i<=maxOrder; i++) {
      Teuchos::Array<long double> nodes(i), weights(i);
      IntrepidBurkardtRules::clenshaw_curtis_compute(i,nodes.getRawPtr(),weights.getRawPtr());
      for (int j=0; j<i; j++) {
	if (j%2)
	  analyticInt = 0.0;
	else 
	  analyticInt = 2.0/((long double)j+1.0);
	testInt = evalQuad(i,j,nodes.getRawPtr(),weights.getRawPtr());
	long double abstol = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	long double absdiff = std::fabs(analyticInt - testInt);
	*outStream << "Cubature order " << std::setw(2) << std::left << i << " integrating "
		   << "x^" << std::setw(2) << std::left << j << ":" << "   "
		   << std::scientific << std::setprecision(16) << testInt << "   " 
		   << analyticInt << "   " << std::setprecision(4) << absdiff << "   " << "<?" 
		   << "   " << abstol << "\n";
	if (absdiff > abstol) {
	  errorFlag++;
	  *outStream << std::right << std::setw(111) << "^^^^---FAILURE!\n";
	}
      }
    }
    *outStream << "\n";    

    *outStream << "Fejer Type 2 Cubature \n";
    *outStream << "Integrates functions on [-1,1] weighted by w(x) = 1\n";
    for (int i = 1; i<=maxOrder; i++) {
      Teuchos::Array<long double> nodes(i), weights(i);
      IntrepidBurkardtRules::fejer2_compute(i,nodes.getRawPtr(),weights.getRawPtr());
      for (int j=0; j<i; j++) {
	if (j%2)
	  analyticInt = 0.0;
	else
	  analyticInt = 2.0/((long double)j+1.0);
	testInt = evalQuad(i,j,nodes.getRawPtr(),weights.getRawPtr());
	long double abstol = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	long double absdiff = std::fabs(analyticInt - testInt);
	*outStream << "Cubature order " << std::setw(2) << std::left << i << " integrating "
		   << "x^" << std::setw(2) << std::left << j << ":" << "   "
		   << std::scientific << std::setprecision(16) << testInt << "   " 
		   << analyticInt << "   " << std::setprecision(4) << absdiff << "   " << "<?" 
		   << "   " << abstol << "\n";
	if (absdiff > abstol) {
	  errorFlag++;
	  *outStream << std::right << std::setw(111) << "^^^^---FAILURE!\n";
	}
      }
    }
    *outStream << "\n";    

    *outStream << "Gauss-Patterson Cubature \n";
    *outStream << "Integrates functions on [-1,1] weighted by w(x) = 1\n";
    for (int l = 1; l<=7; l++) {
      int i = (int)pow(2.0,(double)l+1.0)-1;
      Teuchos::Array<long double> nodes(i), weights(i);
      IntrepidBurkardtRules::patterson_lookup(i,nodes.getRawPtr(),weights.getRawPtr());
      for (int j=0; j<=(1.5*i+0.5); j++) {
	if (j%2)
	  analyticInt = 0.0;
	else
	  analyticInt = 2.0/((long double)j+1.0);
	testInt = evalQuad(i,j,nodes.getRawPtr(),weights.getRawPtr());
	long double abstol = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	long double absdiff = std::fabs(analyticInt - testInt);
	*outStream << "Cubature order " << std::setw(2) << std::left << i << " integrating "
		   << "x^" << std::setw(2) << std::left << j << ":" << "   "
		   << std::scientific << std::setprecision(16) << testInt << "   " 
		   << analyticInt << "   " << std::setprecision(4) << absdiff << "   " << "<?" 
		   << "   " << abstol << "\n";
	if (absdiff > abstol) {
	  errorFlag++;
	  *outStream << std::right << std::setw(111) << "^^^^---FAILURE!\n";
	}
      }
    }
    *outStream << "\n";

    *outStream << "Gauss-Chebyshev Type 1 Cubature \n";
    *outStream << "Integrates functions on [-1,1] weighted by w(x) = 1/sqrt(1-x^2)\n";
    for (int i = 1; i<=maxOrder; i++) {
      Teuchos::Array<long double> nodes(i), weights(i);
      IntrepidBurkardtRules::chebyshev1_compute(i,nodes.getRawPtr(),weights.getRawPtr());      
      for (int j=0; j<=2*i-1; j++) {
	analyticInt = chebyshev1<long double>(j);
	testInt = evalQuad(i,j,nodes.getRawPtr(),weights.getRawPtr());
	long double abstol = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	long double absdiff = std::fabs(analyticInt - testInt);
	*outStream << "Cubature order " << std::setw(2) << std::left << i << " integrating "
		   << "x^" << std::setw(2) << std::left << j << ":" << "   "
		   << std::scientific << std::setprecision(16) << testInt << "   " 
		   << analyticInt << "   " << std::setprecision(4) << absdiff << "   " << "<?" 
		   << "   " << abstol << "\n";
	if (absdiff > abstol) {
	  errorFlag++;
	  *outStream << std::right << std::setw(111) << "^^^^---FAILURE!\n";
	}
      }
    }
    *outStream << "\n";

    *outStream << "Gauss-Chebyshev Type 2 Cubature \n";
    *outStream << "Integrates functions on [-1,1] weighted by w(x) = sqrt(1-x^2)\n";
    for (int i = 1; i<=maxOrder; i++) {
      Teuchos::Array<long double> nodes(i), weights(i);
      IntrepidBurkardtRules::chebyshev2_compute(i,nodes.getRawPtr(),weights.getRawPtr());      
      for (int j=0; j<=2*i-1; j++) {
	analyticInt = chebyshev2<long double>(j);
	testInt = evalQuad(i,j,nodes.getRawPtr(),weights.getRawPtr());
	long double abstol = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	long double absdiff = std::fabs(analyticInt - testInt);
	*outStream << "Cubature order " << std::setw(2) << std::left << i << " integrating "
		   << "x^" << std::setw(2) << std::left << j << ":" << "   "
		   << std::scientific << std::setprecision(16) << testInt << "   " 
		   << analyticInt << "   " << std::setprecision(4) << absdiff << "   " << "<?" 
		   << "   " << abstol << "\n";
	if (absdiff > abstol) {
	  errorFlag++;
	  *outStream << std::right << std::setw(111) << "^^^^---FAILURE!\n";
	}
      }
    }
    *outStream << "\n";
    
     *outStream << "Gauss-Laguerre Cubature \n";
    *outStream << "Integrates functions on [0,oo) weighted by w(x) = exp(-x)\n";
    for (int i = 1; i<=maxOrder; i++) {
      Teuchos::Array<long double> nodes(i), weights(i);
      IntrepidBurkardtRules::laguerre_compute(i,nodes.getRawPtr(),weights.getRawPtr());
       for (int j=0; j<=2*i-1; j++) {
	analyticInt = tgammal((long double)(j+1));
	testInt = evalQuad(i,j,nodes.getRawPtr(),weights.getRawPtr());
	long double abstol = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	long double absdiff = std::fabs(analyticInt - testInt);
	*outStream << "Cubature order " << std::setw(2) << std::left << i << " integrating "
		   << "x^" << std::setw(2) << std::left << j << ":" << "   "
		   << std::scientific << std::setprecision(16) << testInt << "   " 
		   << analyticInt << "   " << std::setprecision(4) << absdiff << "   " << "<?" 
		   << "   " << abstol << "\n";
	if (absdiff > abstol) {
	  errorFlag++;
	  *outStream << std::right << std::setw(111) << "^^^^---FAILURE!\n";
	}
      }
    }
    *outStream << "\n";

    maxOrder = 10;

    *outStream << "Gauss-Hermite Cubature \n";
    *outStream << "Integrates functions on (-oo,oo) weighted by w(x) = exp(-x^2)\n";
    for (int i = 1; i<=maxOrder; i++) {
      Teuchos::Array<long double> nodes(i), weights(i);
      IntrepidBurkardtRules::hermite_compute(i,nodes.getRawPtr(),
					     weights.getRawPtr());      
      for (int j=0; j<=2*i-1; j++) { 
	if (j%2)
	  analyticInt = 0.0;
	else
	  analyticInt = factorial2<long double>(j-1)*sqrt(M_PI)/powl(2.0,(long double)j/2.0);
	testInt = evalQuad(i,j,nodes.getRawPtr(),weights.getRawPtr());
	long double abstol = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	long double absdiff = std::fabs(analyticInt - testInt);
	*outStream << "Cubature order " << std::setw(2) << std::left << i 
		   << " integrating "
		   << "x^" << std::setw(2) << std::left << j << ":" 
		   << "   "
		   << std::scientific << std::setprecision(16) << testInt 
		   << "   " 
		   << analyticInt << "   " << std::setprecision(4) 
		   << absdiff << "   " << "<?" 
		   << "   " << abstol << "\n";
	if (absdiff > abstol) {
	  errorFlag++;
	  *outStream << std::right << std::setw(111) << "^^^^---FAILURE!\n";
	}
      }
    }
    *outStream << "\n";    

    reltol = 1e-6;

    *outStream << "Hermite-Genz-Keister Cubature \n";
    *outStream << "Integrates functions on (-oo,oo) weighted by w(x) = exp(-x^2)\n";
    int order[4] = {1,3, 9,19};
    int max[4]   = {1,5,15,29};
    for (int l = 0; l<4; l++) {
      int i = order[l];
      int m = max[l];
      Teuchos::Array<long double> nodes(i), weights(i);
      IntrepidBurkardtRules::hermite_genz_keister_lookup(i,nodes.getRawPtr(),
							 weights.getRawPtr());  
      for (int j=0; j<=m; j++) { 
	if (j%2)
	  analyticInt = 0.0;
	else
	  analyticInt = factorial2<long double>(j-1)*sqrt(M_PI)/powl(2.0,(long double)j/2.0);
	if (i>=36)
	  analyticInt /= sqrt(M_PI);
	testInt = evalQuad(i,j,nodes.getRawPtr(),weights.getRawPtr());
	long double abstol = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	long double absdiff = std::fabs(analyticInt - testInt);
	*outStream << "Cubature order " << std::setw(2) << std::left << i 
		   << " integrating "
		   << "x^" << std::setw(2) << std::left << j << ":" << "   "
		   << std::scientific << std::setprecision(16) << testInt 
		   << "   " 
		   << analyticInt << "   " << std::setprecision(4) << absdiff 
		   << "   " << "<?" 
		   << "   " << abstol << "\n";
	if (absdiff > abstol) {
	  errorFlag++;
	  *outStream << std::right << std::setw(111) << "^^^^---FAILURE!\n";
	}
      }
    }
    *outStream << "\n";
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

