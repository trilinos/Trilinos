#include <cstdlib>
#include <cassert>
#include <iostream>
#include <strstream>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <complex>
using namespace std;


#include "Petra_Comm.h"
#include "Petra_Time.h"
#include "Petra_Map.h" 

#include "TPetra_ScalarTraits.h"

// Local prototypes
template<class scalarType>
bool check(scalarType var, bool verbose);

int main(int argc, char *argv[]) {

#ifdef PETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Petra_Comm Comm(MPI_COMM_WORLD);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  Petra_Comm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
		    << " is alive."<<endl;

  bool verbose1 = verbose;
  verbose = (MyPID==0);  // Only print most results on PE 0

  int ierr = 0;
  float fvar = 2.0;
  if (!check(fvar, verbose)) ierr++;
  else if (verbose) cout << "***** TPetra::ScalarTraits for type <float> check OK" << endl;

  double dvar = 2.0;
  if (!check(dvar, verbose)) ierr++;
  else if (verbose) cout << "***** TPetra::ScalarTraits for type <double> check OK" << endl;
  
  std::complex<float> cvar = std::complex<float>(2.0, 3.0);
  if (!check(cvar, verbose)) ierr++;
  else if (verbose) cout << "***** TPetra::ScalarTraits for type complex<float> check OK" << endl;
  
  std::complex<double> zvar = std::complex<double>(2.0, 3.0);
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
    cout << "\n***** Scalar Traits for type:  " << TPetra::ScalarTraits<scalarType>::name() 
	 << " *****" << endl << endl
	 << "      Zero = " << zero << endl 
	 << "      One  = " << one << endl
	 << "      First 5 random numbers are: " << endl;
    for (int i=0; i<5; i++) 
      cout << "        " << randomNumbers[i] << endl;
  }
    

  // Confirm magnitude traits
  myMagnitudeType eps = TPetra::ScalarTraits<myMagnitudeType>::eps();
  myMagnitudeType sfmin = TPetra::ScalarTraits<myMagnitudeType>::sfmin();
  myMagnitudeType base = TPetra::ScalarTraits<myMagnitudeType>::base();
  myMagnitudeType prec = TPetra::ScalarTraits<myMagnitudeType>::prec();
  myMagnitudeType t = TPetra::ScalarTraits<myMagnitudeType>::t();
  myMagnitudeType rnd = TPetra::ScalarTraits<myMagnitudeType>::rnd();
  myMagnitudeType emin = TPetra::ScalarTraits<myMagnitudeType>::emin();
  myMagnitudeType rmin = TPetra::ScalarTraits<myMagnitudeType>::rmin();
  myMagnitudeType emax = TPetra::ScalarTraits<myMagnitudeType>::emax();
  myMagnitudeType rmax = TPetra::ScalarTraits<myMagnitudeType>::rmax();
  if (verbose) {
    cout << "\n***** Magnitude Traits for type:  " << TPetra::ScalarTraits<scalarType>::name() 
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
