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
