// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
\brief  Unit test (CubatureDirect): correctness of
        integration of monomials for 1D reference cells.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid2_CubatureTensorSorted.hpp"
//#include "Intrepid2_CubatureLineSorted.hpp"
#include "Intrepid2_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace Intrepid2;

/*
  Computes integrals of monomials over a given reference cell.
*/
long double evalQuad(std::vector<int> power, int dimension, int order, 
		     std::vector<EIntrepidBurkardt> rule, 
		     std::vector<EIntrepidGrowth> growth) {
 
  CubatureTensorSorted<long double> lineCub(dimension,order,rule,growth,false);
  int size = lineCub.getNumPoints();
  FieldContainer<long double> cubPoints(size,dimension);
  FieldContainer<long double> cubWeights(size);
  lineCub.getCubature(cubPoints,cubWeights);

  //  for (int i=0; i<size; i++) {
  //    std::cout << cubPoints(i,0) << "  " << cubPoints(i,1) << std::endl;
  //  }

  long double Q  = 0.0;
  long double Qi = 0.0;
  int l1      = growthRule1D(order,growth[0],rule[0]);
  int l2      = growthRule1D(order,growth[1],rule[1]);
  int mid2    = l2/2;
  int locMid  = 0;
  int cnt = 0;
  for (int i=0; i<l1; i++) {
    locMid = i*l1+mid2; Qi = 0.0;
    if (l2%2) {
      Qi = cubWeights(locMid)*powl(cubPoints(locMid,1),power[1]); cnt++;
      for (int j=1; j<=mid2; j++) {
	Qi += cubWeights(locMid-j)*powl(cubPoints(locMid-j,1),power[1]) 
	  +cubWeights(locMid+j)*powl(cubPoints(locMid+j,1),power[1]); cnt += 2;
      }
    }
    else {
      for (int j=0; j<mid2; j++) {
	Qi += cubWeights(locMid-j)*powl(cubPoints(locMid-j,1),power[1]) 
	  +cubWeights(locMid+j+1)*powl(cubPoints(locMid+j+1,1),power[1]); cnt += 2;
      }
    }
    Qi *= powl(cubPoints(locMid,0),power[0]);
    Q  += Qi;
  }
  return Q;
  /*
  int mid  = size/2;
  long double Q = 0.0;
  if (size%2) {
    Q  = cubWeights(mid);
    for (int i=0; i<dimension; i++) {
      Q *= powl(cubPoints(mid,i),power[i]);
    }
  }

  for (int i=0; i<mid; i++) {
    long double value1 = cubWeights(i);
    long double value2 = cubWeights(size-i-1);
    for (int j=0; j<dimension; j++) {
      value1 *= powl(cubPoints(i,j),power[j]);
      value2 *= powl(cubPoints(size-i-1,j),power[j]);
    }
    Q += value1+value2;
  }
  return Q;
  */
}

long double evalInt(int dimension, std::vector<int> power, std::vector<EIntrepidBurkardt> rule) {
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
  << "| TEST 14: integrals of monomials in 2D - Isotropic with growth rules         |\n"\
  << "===============================================================================\n";


  // internal variables:
  int dimension = 2;
  int         errorFlag   = 0;
  long double reltol      = 1.0e+06*INTREPID_TOL;
  int         maxDegx     = 0; 
  int         maxDegy     = 0; 
  long double analyticInt = 0;
  long double testInt     = 0;
  int         maxOrder    = 3;
  int         l1 = 0, l2  = 0;
  std::vector<int> power(2,0);
  std::vector<EIntrepidBurkardt> rule1(2,BURK_CLENSHAWCURTIS);
  int order = 0;;
  std::vector<EIntrepidGrowth> growth(2,GROWTH_DEFAULT);

  *outStream << "\nIntegrals of monomials on a reference line (edge):\n";
  // compute and compare integrals
  try {
    for (EIntrepidBurkardt rule=BURK_CHEBYSHEV1; rule <= BURK_LAGUERRE; rule++) {   
      // compute integrals
      rule1[0] = rule; rule1[1] = rule;
      if (rule!=BURK_PATTERSON&&rule!=BURK_GENZKEISTER&&rule!=BURK_TRAPEZOIDAL) {
	*outStream << "Testing " << EIntrepidBurkardtToString(rule) << "\n";
	for (int i=1; i <= maxOrder; i++) {
	  l1 = growthRule1D(i,growth[0],rule);
	  l2 = growthRule1D(i,growth[1],rule);
	  if (rule==BURK_CHEBYSHEV1||rule==BURK_CHEBYSHEV2||rule==BURK_LEGENDRE||
	      rule==BURK_LAGUERRE||rule==BURK_HERMITE) {
	    maxDegx = 2*l1-1;
	    maxDegy = 2*l2-1;
	  }
	  else if (rule==BURK_CLENSHAWCURTIS||rule==BURK_FEJER2) {
	    maxDegx = l1-1;
	    maxDegy = l2-1;
	  }
	  
	  order = i;
	  for (int j=0; j <= maxDegx; j++) {
	    power[0] = j;
	    for (int k=0; k <= maxDegy; k++) {
	      power[1] = k;
	      analyticInt = evalInt(dimension, power, rule1);
	      testInt     = evalQuad(power,dimension,order,rule1,growth);
	      
	      long double abstol  = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	      long double absdiff = std::fabs(analyticInt - testInt);
	      *outStream << "Cubature order (" << std::setw(2) << std::left 
			 << l1 << ", " << std::setw(2) << std::left << l2
			 << ") integrating "
			 << "x^" << std::setw(2) << std::left << j << "y^" 
			 << std::setw(2) << std::left 
			 << k <<  ":" << "   " << std::scientific 
			 << std::setprecision(16) << testInt 
			 << "   " << analyticInt << "   " 
			 << std::setprecision(4) << absdiff << "   " 
			 << "<?" << "   " << abstol << "\n";
	      if (absdiff > abstol) {
		errorFlag++;
		*outStream << std::right << std::setw(104) 
			   << "^^^^---FAILURE!\n";
	      }
	    } // end for k
	    *outStream << "\n";
	  } // end for j
	  *outStream << "\n";
	} // end for i
      }
      else if (rule==BURK_PATTERSON) {
	*outStream << "Testing " << EIntrepidBurkardtToString(rule) << "\n";
	for (int i=0; i < 3; i++) {	  
	  l1 = growthRule1D(i,growth[0],rule);
	  l2 = growthRule1D(i,growth[1],rule);
	  if (i==0) { 
	    maxDegx = 1;
	    maxDegy = 1;
	  }
	  else {
	    maxDegx = (int)(1.5*(double)l1-0.5);
	    maxDegy = (int)(1.5*(double)l2-0.5);
	  }

	  order = i;
	  for (int j=0; j <= maxDegx; j++) {	    
	    power[0] = j;
	    for (int k=0; k <= maxDegy; k++) {	      
	      power[1] = k;
	      analyticInt = evalInt(dimension, power, rule1);
	      testInt     = evalQuad(power,dimension,order,rule1,growth);
	      
	      long double abstol  = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	      long double absdiff = std::fabs(analyticInt - testInt);
	      *outStream << "Cubature order (" << std::setw(2) << std::left 
			 << l1 << ", " << std::setw(2) << std::left << l2 
			 << ") integrating "
			 << "x^" << std::setw(2) << std::left << j << "y^" 
			 << std::setw(2) << std::left 
			 << k <<  ":" << "   " << std::scientific 
			 << std::setprecision(16) << testInt 
			 << "   " << analyticInt << "   " 
			 << std::setprecision(4) << absdiff << "   " 
			 << "<?" << "   " << abstol << "\n";
	      if (absdiff > abstol) {
		errorFlag++;
		*outStream << std::right << std::setw(104) 
			   << "^^^^---FAILURE!\n";
	      }
	    } // end for k
	    *outStream << "\n";
	  } // end for j
	  *outStream << "\n";
	} // end for i
      }
      else if (rule==BURK_GENZKEISTER) {
	*outStream << "Testing " << EIntrepidBurkardtToString(rule) << "\n";
	for (int i=0; i < 3; i++) {
	  l1 = growthRule1D(i,growth[0],rule);
	  l2 = growthRule1D(i,growth[1],rule);
	  if (i==0) { 
	    maxDegx = 1;
	    maxDegy = 1;
	  }
	  else {
	    maxDegx = (int)(1.5*(double)l1-0.5);
	    maxDegy = (int)(1.5*(double)l2-0.5);
	  }

	  order = i;
	  for (int j=0; j <= maxDegx; j++) {	    
	    power[0] = j;	    
	    for (int k=0; k <= maxDegy; k++) {	      
	      power[1] = k;
	      analyticInt = evalInt(dimension, power, rule1);
	      testInt     = evalQuad(power,dimension,order,rule1,growth);
	      
	      long double abstol  = (analyticInt == 0.0 ? reltol : std::fabs(reltol*analyticInt) );
	      long double absdiff = std::fabs(analyticInt - testInt);
	      *outStream << "Cubature order (" << std::setw(2) << std::left 
			 << l1 << ", " << std::setw(2) << std::left << l2
			 << ") integrating "
			 << "x^" << std::setw(2) << std::left << j << "y^" 
			 << std::setw(2) << std::left 
			 << k <<  ":" << "   " << std::scientific 
			 << std::setprecision(16) << testInt 
			 << "   " << analyticInt << "   " 
			 << std::setprecision(4) << absdiff << "   " 
			 << "<?" << "   " << abstol << "\n";
	      if (absdiff > abstol) {
		errorFlag++;
		*outStream << std::right << std::setw(104) 
			   << "^^^^---FAILURE!\n";
	      }
	    } // end for k
	    *outStream << "\n";
	  } // end for j
	  *outStream << "\n";
	} // end for i
	}
    } // end for rule
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
Kokkos::finalize();
  return errorFlag;
}
