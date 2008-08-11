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

#include "Tpetra_ScalarTraits.hpp"
#include "Tpetra_Version.hpp"

#define PACKETTYPE int
#define ORDINALTYPE int

// Local prototypes
template<typename ScalarType>
bool check(ScalarType var, bool verbose);

int main(int argc, char *argv[]) {
  bool verbose = false;
	bool debug = false;
  // Check if we should print results to standard out
  if (argc>1) {
		if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;
		if (argv[1][0]=='-' && argv[1][1]=='d') {
			verbose = true;
			debug = true;
		}
	}

  if (verbose)
	cout << Tpetra::Tpetra_Version() << endl << endl;

  int ierr = 0;
  float fvar = 2.0;
  if (!check(fvar, debug)) ierr++;
  else if (verbose) cout << "ScalarTraits for type <float> successful." << endl;

  double dvar = 2.0;
  if (!check(dvar, debug)) ierr++;
  else if (verbose) cout << "ScalarTraits for type <double> successful." << endl;
  
  std::complex<float> cvar = std::complex<float>(2.0, 3.0);
  if (!check(cvar, debug)) ierr++;
  else if (verbose) cout << "ScalarTraits for type complex<float> successful." << endl;
  
  std::complex<double> zvar = std::complex<double>(2.0, 3.0);
  if (!check(zvar, debug)) ierr++;
  else if (verbose) cout << "ScalarTraits for type complex<double> successful." << endl;



	if(verbose) cout << "ScalarTraits test successful." << endl;
  
  return 0; // All done
}
template<typename ScalarType>
bool check(ScalarType var, bool verbose) {
  
  // Create shorthand for scalar and magnitude types
  typedef ScalarType myScalarType;
  typedef typename Tpetra::ScalarTraits<ScalarType>::magnitudeType myMagnitudeType;

  myScalarType zero = Tpetra::ScalarTraits<ScalarType>::zero();
  myScalarType one = Tpetra::ScalarTraits<ScalarType>::one();
  myScalarType* randomNumbers = new myScalarType[5];
  for(int i=0; i<5; i++) 
		randomNumbers[i] = Tpetra::ScalarTraits<ScalarType>::random();

  // Print out myScalarType results
  if (verbose) {
    cout << "\n***** Scalar Traits for type:  " << Tpetra::ScalarTraits<ScalarType>::name() 
	 << " *****" << endl << endl
	 << "      Zero = " << zero << endl 
	 << "      One  = " << one << endl
	 << "      First 5 random numbers are: " << endl;
    for (int i=0; i<5; i++) 
      cout << "        " << randomNumbers[i] << endl;
  }
    
  // Confirm magnitude traits
  myMagnitudeType eps = Tpetra::ScalarTraits<myMagnitudeType>::eps();
  myMagnitudeType sfmin = Tpetra::ScalarTraits<myMagnitudeType>::sfmin();
  myMagnitudeType base = Tpetra::ScalarTraits<myMagnitudeType>::base();
  myMagnitudeType prec = Tpetra::ScalarTraits<myMagnitudeType>::prec();
  myMagnitudeType t = Tpetra::ScalarTraits<myMagnitudeType>::t();
  myMagnitudeType rnd = Tpetra::ScalarTraits<myMagnitudeType>::rnd();
  myMagnitudeType emin = Tpetra::ScalarTraits<myMagnitudeType>::emin();
  myMagnitudeType rmin = Tpetra::ScalarTraits<myMagnitudeType>::rmin();
  myMagnitudeType emax = Tpetra::ScalarTraits<myMagnitudeType>::emax();
  myMagnitudeType rmax = Tpetra::ScalarTraits<myMagnitudeType>::rmax();
  if (verbose) {
    cout << "\n***** Magnitude Traits for type:  " << Tpetra::ScalarTraits<ScalarType>::name() 
	 << " *****" << endl << endl
      << "      Relative Machine Epsilon                            (eps)   =  " << eps << endl 
      << "      Safe minimum, such that 1/sfmin does not overflow   (sfmin) =  " << sfmin << endl 
      << "      Base of the machine                                 (base)  =  " << base << endl 
      << "      eps*base                                            (prec)  =  " << prec << endl 
      << "      Number of (base) digits in the mantissa             (t)     =  " << t << endl 
      << "      1.0 when rounding occurs in addition, 0.0 otherwise (rnd)   =  " << rnd << endl 
      << "      Minimum exponent before (gradual) underflow         (emin)  =  " << emin << endl 
      << "      Underflow threshold - base**(emin-1)                (rmin)  =  " << rmin << endl 
      << "      Largest exponent before overflow                    (emax)  =  " << emax << endl 
      << "      Overflow threshold  - (base**emax)*(1-eps)          (rmax)  =  " << rmax << endl;
  }

  return(true);
}
