// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <strstream>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <complex>
using namespace std;


#include "TPetra_ScalarTraits.h"

// Local prototypes
template<class scalarType>
bool check(scalarType var, bool verbose);

int main(int argc, char *argv[]) {

  bool verbose = true;
  int ierr = 0;
  float fvar = 2.0;
  if (!check(fvar, verbose)) ierr++;
  else if (verbose) cout << "***** TPetra::ScalarTraits for type <float> check OK" << endl;

  double dvar = 2.0;
  if (!check(dvar, verbose)) ierr++;
  else if (verbose) cout << "***** TPetra::ScalarTraits for type <double> check OK" << endl;

  complex<float> cvar = complex<float>(2.0, 3.0);
  if (!check(cvar, verbose)) ierr++;
  else if (verbose) cout << "***** TPetra::ScalarTraits for type complex<float> check OK" << endl;
  
  complex<double> zvar = complex<double>(2.0, 3.0);
  if (!check(zvar, verbose)) ierr++;
  else if (verbose) cout << "***** TPetra::ScalarTraits for type complex<double> check OK" << endl;
  
  return 0; // All done
}
template<class scalarType>
bool check(scalarType var, bool verbose) {
  
  // Create shorthand for scalar and magnitude types
  typedef scalarType myScalarType;
  typedef typename TPetra::ScalarTraits<scalarType>::magnitudeType myMagnitudeType;

  myScalarType zero = TPetra::ScalarTraits<scalarType>::zero();
  myScalarType one = TPetra::ScalarTraits<scalarType>::one();
  myScalarType * randomNumbers = new myScalarType[5];
  for (int i=0; i<5; i++) randomNumbers[i] = TPetra::ScalarTraits<scalarType>::random();

  // Print out myScalarType results

  if (verbose) {
    cout << "***** Scalar Traits for type:  " << TPetra::ScalarTraits<scalarType>::name() << endl
      << "      Zero = " << zero << endl 
      << "      One  = " << one << endl
	 << "      First 5 random numbers are: " << endl;
    for (int i=0; i<5; i++) 
      cout << "        " << randomNumbers[i] << endl;
  }
  return(true);
}
