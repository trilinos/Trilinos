// Kris
// 07.24.03 -- Initial checkin
// 08.08.03 -- All test suites except for TRSM are finished.

/*

This test program is intended to check an experimental default type (e.g. mp_real) against an "officialy supported" control type (e.g. double). For each test, the program will generate the appropriate scalars and randomly-sized vectors and matrices of random numbers for the routine being tested. All of the input data for the experimental type is casted into the control type, so in theory both BLAS routines should get the same input data. Upon return, the experimental output data is casted back into the control type, and the results are compared; if they are equal (within a user-definable tolerance) the test is considered successful.

The test routine for TRSM is still being developed; all of the others are more or less finalized.

*/

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include "Teuchos_BLAS.hpp"
#include "mp/mpreal.h"

using namespace std;
using namespace Teuchos;

// SType1 and SType2 define the datatypes for which BLAS output will be compared.
// SType2 should generally be a control datatype "officially" supported by the BLAS; SType1 should be the experimental type being checked.
#define SType1     float
#define SType2     double
// TOL defines the tolerance allowed for differences in BLAS output. Set to 0 for strict comparisons.
#define TOL        1
// MVMIN/MAX define the minimum and maximum dimensions of generated matrices and vectors, respectively.
#define MVMIN      2
#define MVMAX      20
// SCALARMAX defines the maximum positive value (with a little leeway) generated for matrix and vector elements and scalars:
// random numbers in [-SCALARMAX, SCALARMAX] will be generated.
// Set SCALARMAX to a floating-point value (e.g. 10.0) to enable floating-point random number generation, such that
// random numbers in (-SCALARMAX - 1, SCALARMAX + 1) will be generated.
// Large values of SCALARMAX may cause problems with SType2 = int, as large integer values will overflow floating-point types.
#define SCALARMAX  10
// These define the number of tests to be run for each individual BLAS routine.
#define ASUMTESTS  0
#define AXPYTESTS  0
#define COPYTESTS  0
#define DOTTESTS   0
#define IAMAXTESTS 0
#define NRM2TESTS  0
#define SCALTESTS  0
#define GEMVTESTS  0
#define GERTESTS   0
#define TRMVTESTS  0
#define GEMMTESTS  0
#define SYMMTESTS  0
#define TRMMTESTS  0
#define TRSMTESTS  1

// Returns ScalarTraits<TYPE>::random() (the input parameters are ignored)
template<typename TYPE>
TYPE GetRandom(TYPE, TYPE);

// Returns a random integer between the two input parameters, inclusive
template<>
int GetRandom(int, int);

// Returns a random double between the two input parameters, plus or minus a random number between 0 and 1
template<>
double GetRandom(double, double);

template<typename TYPE>
void PrintVector(TYPE*, int, string, bool = 0);

template<typename TYPE>
void PrintMatrix(TYPE*, int, int, string, bool = 0);

template<typename TYPE1, typename TYPE2>
bool CompareScalars(TYPE1, TYPE2, double = 0);

template<typename TYPE1, typename TYPE2>
bool CompareVectors(TYPE1*, TYPE2*, int, double = 0);

template<typename TYPE1, typename TYPE2>
bool CompareMatrices(TYPE1*, TYPE2*, int, int, double = 0);

// For most types, this function is just a wrapper for static_cast(), but for mp_real/double, it calls mp::dble()
// The second input parameter is not used; it is only needed to determine what type to convert *to*
template<typename TYPE1, typename TYPE2>
TYPE2 ConvertType(TYPE1, TYPE2);

template<>
double ConvertType(mp_real, double);

// These functions return a random character appropriate for the BLAS arguments that share their names (uses GetRandom())
char RandomSIDE();
char RandomUPLO();
char RandomTRANS();
char RandomDIAG();

int main(int argc, char *argv[])
{
  bool verbose = 0;
  bool debug = 0;
  bool matlab = 0;
  bool InvalidCmdLineArgs = 0;
  int i, j;
  for(i = 1; i < argc; i++)
    {
      if(argv[i][0] == '-')
	{
	  switch(argv[i][1])
	    {
	    case 'v':
	      if(!verbose)
		{
		  verbose = 1;
		}
	      else
		{
		  InvalidCmdLineArgs = 1;
		}
	      break;
	    case 'd':
	      if(!debug)
		{
		  debug = 1;
		}
	      else
		{
		  InvalidCmdLineArgs = 1;
		}
	      break;
	    case 'm':
	      if(!matlab)
		{
		  matlab = 1;
		}
	      else
		{
		  InvalidCmdLineArgs = 1;
		}
	      break;
	    default:
	      InvalidCmdLineArgs = 1;
	      break;
	    }
	}
    }
  if(InvalidCmdLineArgs || (argc > 4))
    {
      cout << "Invalid command line arguments detected. Use the following flags:" << endl
	   << "\t -v enables verbose mode (reports number of failed/successful tests)" << endl
	   << "\t -d enables debug mode (same as verbose with output of each test, not recommended for large numbers of tests)" << endl
	   << "\t -m enables matlab-style output; only has an effect if debug mode is enabled" << endl;
      return 1;
    }
  BLAS<int, SType1> SType1BLAS;
  BLAS<int, SType2> SType2BLAS;
  SType1 SType1zero = ScalarTraits<SType1>::zero();
  SType1 SType1one = ScalarTraits<SType1>::one();
  SType2 SType2zero = ScalarTraits<SType2>::zero();
  SType2 SType2one = ScalarTraits<SType2>::one();
  SType1* SType1A;
  SType1* SType1B;
  SType1* SType1C;
  SType1* SType1x;
  SType1* SType1y;
  SType1 SType1alpha, SType1beta;
  SType2* SType2A;
  SType2* SType2B;
  SType2* SType2C;
  SType2* SType2x;
  SType2* SType2y; 
  SType2 SType2alpha, SType2beta;
  SType1 SType1ASUMresult, SType1DOTresult, SType1NRM2result;
  SType2 SType2ASUMresult, SType2DOTresult, SType2NRM2result;
  int SType1IAMAXresult;
  int SType2IAMAXresult;
  int TotalTestCount = 1, GoodTestSubcount, GoodTestCount = 0, M, N, P, LDA, LDB;
  char UPLO, SIDE, TRANSA, TRANSB, DIAG;
  SType2 convertTo;

  srand(time(NULL));
  mp::mp_init(200);

  // Begin ASUM Tests
  GoodTestSubcount = 0;
  for(i = 0; i < ASUMTESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      SType1x = new SType1[M];
      SType2x = new SType2[M];
      for(j = 0; j < M; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  PrintVector(SType1x, M, "SType1x", matlab);
	  PrintVector(SType2x, M, "SType2x", matlab);
	}
      TotalTestCount++;
      SType1ASUMresult = SType1BLAS.ASUM(M, SType1x, 1);
      SType2ASUMresult = SType2BLAS.ASUM(M, SType2x, 1);
      if(debug)
	{
	  cout << "SType1 ASUM result: " << SType1ASUMresult << endl;
	  cout << "SType2 ASUM result: " << SType2ASUMresult << endl;
	}
      GoodTestSubcount += CompareScalars(SType1ASUMresult, SType2ASUMresult, TOL);
      delete [] SType1x;
      delete [] SType2x;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) cout << "ASUM: " << GoodTestSubcount << " of " << ASUMTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End ASUM Tests

  // Begin AXPY Tests
  GoodTestSubcount = 0;
  for(i = 0; i < AXPYTESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      SType1x = new SType1[M];
      SType1y = new SType1[M];
      SType2x = new SType2[M];
      SType2y = new SType2[M]; 
      for(j = 0; j < M; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	  SType1y[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2y[j] = ConvertType(SType1y[j], convertTo);
	}
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  cout << "SType1alpha = "  << SType1alpha << endl;
	  cout << "SType2alpha = " << SType2alpha << endl;
	  PrintVector(SType1x, M, "SType1x", matlab);
	  PrintVector(SType1y, M, "SType1y_before_operation", matlab);
	  PrintVector(SType2x, M, "SType2x", matlab);
	  PrintVector(SType2y, M, "SType2y_before_operation",  matlab);
	}
      TotalTestCount++;
      SType1BLAS.AXPY(M, SType1alpha, SType1x, 1, SType1y, 1);
      SType2BLAS.AXPY(M, SType2alpha, SType2x, 1, SType2y, 1);
      if(debug)
	{
	  PrintVector(SType1y, M, "SType1y_after_operation", matlab);
	  PrintVector(SType2y, M, "SType2y_after_operation", matlab);
	}
      GoodTestSubcount += CompareVectors(SType1y, SType2y, M, TOL);
      delete [] SType1x;
      delete [] SType1y;
      delete [] SType2x;
      delete [] SType2y;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) cout << "AXPY: " << GoodTestSubcount << " of " << AXPYTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End AXPY Tests

  // Begin COPY Tests
  GoodTestSubcount = 0;
  for(i = 0; i < COPYTESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      SType1x = new SType1[M];
      SType1y = new SType1[M];
      SType2x = new SType2[M];
      SType2y = new SType2[M]; 
      for(j = 0; j < M; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	  SType1y[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2y[j] = ConvertType(SType1y[j], convertTo);
	}
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  PrintVector(SType1x, M, "SType1x", matlab);
	  PrintVector(SType1y, M, "SType1y_before_operation", matlab);
	  PrintVector(SType2x, M, "SType2x", matlab);
	  PrintVector(SType2y, M, "SType2y_before_operation", matlab);
	}
      TotalTestCount++;
      SType1BLAS.COPY(M, SType1x, 1, SType1y, 1);
      SType2BLAS.COPY(M, SType2x, 1, SType2y, 1);
      if(debug)
	{
	  PrintVector(SType1y, M, "SType1y_after_operation", matlab);
	  PrintVector(SType2y, M, "SType2y_after_operation", matlab);
	}
      GoodTestSubcount += CompareVectors(SType1y, SType2y, M, TOL);
      delete [] SType1x;
      delete [] SType1y;
      delete [] SType2x;
      delete [] SType2y;
    }
   GoodTestCount += GoodTestSubcount; if(verbose || debug) cout << "COPY: " << GoodTestSubcount << " of " << COPYTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End COPY Tests

  // Begin DOT Tests
  GoodTestSubcount = 0;
  for(i = 0; i < DOTTESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      SType1x = new SType1[M];
      SType1y = new SType1[M];
      SType2x = new SType2[M];
      SType2y = new SType2[M]; 
      for(j = 0; j < M; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	  SType1y[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2y[j] = ConvertType(SType1y[j], convertTo);
	}
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  PrintVector(SType1x, M, "SType1x", matlab);
	  PrintVector(SType1y, M, "SType1y", matlab);
	  PrintVector(SType2x, M, "SType2x", matlab);
	  PrintVector(SType2y, M, "SType2y", matlab);
	}
      TotalTestCount++;
      SType1DOTresult = SType1BLAS.DOT(M, SType1x, 1, SType1y, 1);
      SType2DOTresult = SType2BLAS.DOT(M, SType2x, 1, SType2y, 1);
      if(debug)
	{
	  cout << "SType1 DOT result: " << SType1DOTresult << endl;
	  cout << "SType2 DOT result: " << SType2DOTresult << endl;
	}
      GoodTestSubcount += CompareScalars(SType1DOTresult, SType2DOTresult, TOL);
      delete [] SType1x;
      delete [] SType1y;
      delete [] SType2x;
      delete [] SType2y;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) cout << "DOT: " << GoodTestSubcount << " of " << DOTTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End DOT Tests

  // Begin NRM2 Tests
  GoodTestSubcount = 0;
  for(i = 0; i < NRM2TESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      SType1x = new SType1[M];
      SType2x = new SType2[M];
      for(j = 0; j < M; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  PrintVector(SType1x, M, "SType1x", matlab);
	  PrintVector(SType2x, M, "SType2x", matlab);
	}
      TotalTestCount++;
      SType1NRM2result = SType1BLAS.NRM2(M, SType1x, 1);
      SType2NRM2result = SType2BLAS.NRM2(M, SType2x, 1);
      if(debug)
	{
	  cout << "SType1 NRM2 result: " << SType1NRM2result << endl;
	  cout << "SType2 NRM2 result: " << SType2NRM2result << endl;
	}
      GoodTestSubcount += CompareScalars(SType1NRM2result, SType2NRM2result, TOL);
      delete [] SType1x;
      delete [] SType2x;
    }
   GoodTestCount += GoodTestSubcount; if(verbose || debug) cout << "NRM2: " << GoodTestSubcount << " of " << NRM2TESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End NRM2 Tests

  // Begin SCAL Tests
  GoodTestSubcount = 0;
  for(i = 0; i < SCALTESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      SType1x = new SType1[M];
      SType2x = new SType2[M];
      for(j = 0; j < M; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  cout << "SType1alpha = " << SType1alpha << endl;
	  cout << "SType2alpha = " << SType2alpha << endl;
	  PrintVector(SType1x, M, "SType1x_before_operation", matlab);
	  PrintVector(SType2x, M, "SType2x_before_operation", matlab);
	}
      TotalTestCount++;
      SType1BLAS.SCAL(M, SType1alpha, SType1x, 1);
      SType2BLAS.SCAL(M, SType2alpha, SType2x, 1);
      if(debug)
	{
	  PrintVector(SType1x, M, "SType1x_after_operation", matlab);
	  PrintVector(SType2x, M, "SType2x_after_operation", matlab);
	}
      GoodTestSubcount += CompareVectors(SType1x, SType2x, M, TOL);
      delete [] SType1x;
      delete [] SType2x;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) cout << "SCAL: " << GoodTestSubcount << " of " << SCALTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End SCAL Tests

  // Begin IAMAX Tests
  GoodTestSubcount = 0;
  for(i = 0; i < IAMAXTESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      SType1x = new SType1[M];
      SType2x = new SType2[M];
      for(j = 0; j < M; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  PrintVector(SType1x, M, "SType1x", matlab);
	  PrintVector(SType2x, M, "SType2x", matlab);
	}
      TotalTestCount++;
      SType1IAMAXresult = SType1BLAS.IAMAX(M, SType1x, 1);
      SType2IAMAXresult = SType2BLAS.IAMAX(M, SType2x, 1);
      if(debug)
	{
	  cout << "SType1 IAMAX result: " << SType1IAMAXresult << endl;
	  cout << "SType2 IAMAX result: " << SType2IAMAXresult << endl;
	}
      GoodTestSubcount += (SType1IAMAXresult == SType2IAMAXresult);
      delete [] SType1x;
      delete [] SType2x;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) cout << "IAMAX: " << GoodTestSubcount << " of " << IAMAXTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End IAMAX Tests

  // Begin GEMV Tests
  GoodTestSubcount = 0;
  for(i = 0; i < GEMVTESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      N = GetRandom(MVMIN, MVMAX);
      SType1A = new SType1[M * N];
      SType1x = new SType1[N];
      SType1y = new SType1[M];
      SType2A = new SType2[M * N];
      SType2x = new SType2[N];
      SType2y = new SType2[M]; 
      for(j = 0; j < M * N; j++)
	{
	  SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2A[j] = ConvertType(SType1A[j], convertTo);
	}
      for(j = 0; j < N; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      for(j = 0; j < M; j++)
	{
	  SType1y[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2y[j] = ConvertType(SType1y[j], convertTo);
	}
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType1beta = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      SType2beta = ConvertType(SType1beta, convertTo);
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  cout << "SType1alpha = " << SType1alpha << endl;
	  cout << "SType2alpha = " << SType2alpha << endl;
	  cout << "SType1beta = " << SType1beta << endl;
	  cout << "SType2beta = " << SType2beta << endl;
	  PrintMatrix(SType1A, M, N, "SType1A", matlab);
	  PrintVector(SType1x, N, "SType1x", matlab);
	  PrintVector(SType1y, M, "SType1y_before_operation", matlab);
	  PrintMatrix(SType2A, M, N, "SType2A", matlab);
	  PrintVector(SType2x, N, "SType2x", matlab);
	  PrintVector(SType2y, M, "SType2y_before_operation", matlab);
	}
      TotalTestCount++;
      SType1BLAS.GEMV('N', M, N, SType1alpha, SType1A, M, SType1x, 1, SType1beta, SType1y, 1);
      SType2BLAS.GEMV('N', M, N, SType2alpha, SType2A, M, SType2x, 1, SType2beta, SType2y, 1);
      if(debug)
	{
	  PrintVector(SType1y, M, "SType1y_after_operation", matlab);
	  PrintVector(SType2y, M, "SType2y_after_operation", matlab);
	}
      GoodTestSubcount += CompareVectors(SType1y, SType2y, M, TOL);
      delete [] SType1A;
      delete [] SType1x;
      delete [] SType1y;
      delete [] SType2A;
      delete [] SType2x;
      delete [] SType2y;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) cout << "GEMV: " << GoodTestSubcount << " of " << GEMVTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End GEMV Tests

  // Begin TRMV Tests
  GoodTestSubcount = 0;
  for(i = 0; i < TRMVTESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      SType1A = new SType1[M * M];
      SType1x = new SType1[M];
      SType2A = new SType2[M * M];
      SType2x = new SType2[M];
      for(j = 0; j < M * M; j++)
	{
	  SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2A[j] = ConvertType(SType1A[j], convertTo);
	}
      for(j = 0; j < M; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  PrintMatrix(SType1A, M, M, "SType1A", matlab);
	  PrintVector(SType1x, M, "SType1x_before_operation", matlab);
	  PrintMatrix(SType2A, M, M, "SType2A", matlab);
	  PrintVector(SType2x, M, "SType2x_before_operation", matlab);
	}
      TotalTestCount++;
      SType1BLAS.TRMV('U', 'N', 'U', M, SType1A, M, SType1x, 1);
      SType2BLAS.TRMV('U', 'N', 'U', M, SType2A, M, SType2x, 1);
      if(debug)
	{
	  PrintVector(SType1x, M, "SType1x_after_operation", matlab);
	  PrintVector(SType2x, M, "SType2x_after_operation", matlab);
	}
      GoodTestSubcount += CompareVectors(SType1x, SType2x, M, TOL);
      delete [] SType1A;
      delete [] SType1x;
      delete [] SType2A;
      delete [] SType2x;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) cout << "TRMV: " << GoodTestSubcount << " of " << TRMVTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End TRMV Tests

  // Begin GER Tests
  GoodTestSubcount = 0;
  for(i = 0; i < GERTESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      N = GetRandom(MVMIN, MVMAX);
      SType1A = new SType1[M * N];
      SType1x = new SType1[M];
      SType1y = new SType1[N];
      SType2A = new SType2[M * N];
      SType2x = new SType2[M];
      SType2y = new SType2[N];
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      for(j = 0; j < M * N; j++)
	{
	  SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2A[j] = ConvertType(SType1A[j], convertTo);
	}
      for(j = 0; j < M; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      for(j = 0; j < N; j++)
	{
	  SType1y[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2y[j] = ConvertType(SType1y[j], convertTo);
	}
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  cout << "SType1alpha = " << SType1alpha << endl;
	  cout << "SType2alpha = " << SType2alpha << endl;
	  PrintMatrix(SType1A, M, N, "SType1A_before_operation", matlab);
	  PrintVector(SType1x, M, "SType1x", matlab);
	  PrintVector(SType1y, N, "SType1y", matlab);
	  PrintMatrix(SType2A, M, N, "SType2A_before_operation", matlab);
	  PrintVector(SType2x, M, "SType2x", matlab);
	  PrintVector(SType2y, N, "SType2y", matlab);
	}
      TotalTestCount++;
      SType1BLAS.GER(M, N, SType1alpha, SType1x, 1, SType1y, 1, SType1A, M);
      SType2BLAS.GER(M, N, SType2alpha, SType2x, 1, SType2y, 1, SType2A, M);
      if(debug)
	{
	  PrintMatrix(SType1A, M, N, "SType1A_after_operation", matlab);
	  PrintMatrix(SType2A, M, N, "SType2A_after_operation", matlab);
	}
      GoodTestSubcount += CompareMatrices(SType1A, SType2A, M, N, TOL);
      delete [] SType1A;
      delete [] SType1x;
      delete [] SType1y;
      delete [] SType2A;
      delete [] SType2x;
      delete [] SType2y;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) cout << "GER: " << GoodTestSubcount << " of " << GERTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End GER Tests

  // Begin GEMM Tests
  GoodTestSubcount = 0;
  for(i = 0; i < GEMMTESTS; i++)
    { 
      M = GetRandom(MVMIN, MVMAX);
      N = GetRandom(MVMIN, MVMAX);
      P = GetRandom(MVMIN, MVMAX);
      SType1A = new SType1[M * P];
      SType1B = new SType1[P * N];
      SType1C = new SType1[M * N];
      SType2A = new SType2[M * P];
      SType2B = new SType2[P * N];
      SType2C = new SType2[M * N];
      for(j = 0; j < M * P; j++)
	{
	  SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2A[j] = ConvertType(SType1A[j], convertTo);
	}
      for(j = 0; j < P * N; j++)
	{
	  SType1B[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2B[j] = ConvertType(SType1B[j], convertTo);
	}
      for(j = 0; j < M * N; j++)
	{
	  SType1C[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2C[j] = ConvertType(SType1C[j], convertTo);
	}
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType1beta = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      SType2beta = ConvertType(SType1beta, convertTo);
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  PrintMatrix(SType1A, M, P, "SType1A", matlab);
	  PrintMatrix(SType1B, P, N, "SType1B", matlab);
	  PrintMatrix(SType1C, M, N, "SType1C_before_operation", matlab);
	  PrintMatrix(SType2A, M, P, "SType2A", matlab);
	  PrintMatrix(SType2B, P, N, "SType2B", matlab);
	  PrintMatrix(SType2C, M, N, "SType2C_before_operation", matlab);
	}
      TotalTestCount++;
      TRANSA = RandomTRANS();
      TRANSB = RandomTRANS();
      if(TRANSA == 'N')
	{
	  LDA = M;
	}
      else
	{
	  LDA = P;
	}
      if(TRANSB == 'N')
	{
	  LDB = P;
	}
      else
	{
	  LDB = N;
	}
      SType1BLAS.GEMM(TRANSA, TRANSB, M, N, P, SType1alpha, SType1A, LDA, SType1B, LDB, SType1beta, SType1C, M);
      SType2BLAS.GEMM(TRANSA, TRANSB, M, N, P, SType2alpha, SType2A, LDA, SType2B, LDB, SType2beta, SType2C, M);
      if(debug)
	{
	  PrintMatrix(SType1C, M, N, "SType1C_after_operation", matlab);
	  PrintMatrix(SType2C, M, N, "SType2C_after_operation", matlab);
	}
      GoodTestSubcount += CompareMatrices(SType1C, SType2C, M, N, TOL);
      delete [] SType1A;
      delete [] SType1B;
      delete [] SType1C;
      delete [] SType2A;
      delete [] SType2B;
      delete [] SType2C;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) cout << "GEMM: " << GoodTestSubcount << " of " << GEMMTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End GEMM Tests

  // Begin SYMM Tests
  GoodTestSubcount = 0;
  for(i = 0; i < SYMMTESTS; i++)
    { 
      M = GetRandom(MVMIN, MVMAX);
      N = GetRandom(MVMIN, MVMAX);
      SType1B = new SType1[M * N];
      SType1C = new SType1[M * N];
      SType2B = new SType2[M * N];
      SType2C = new SType2[M * N];

      SIDE = RandomSIDE();
      if(SIDE == 'L')
	{
	  SType1A = new SType1[M * M];
	  SType2A = new SType2[M * M];
	  for(j = 0; j < M * M; j++)
	    {
	      SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	      SType2A[j] = ConvertType(SType1A[j], convertTo);
	    }
	  LDA = M;
	}
      else
	{
	  SType1A = new SType1[N * N];
	  SType2A = new SType2[N * N];
	  for(j = 0; j < N * N; j++)
	    {
	      SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	      SType2A[j] = ConvertType(SType1A[j], convertTo);
	    }
	  LDA = N;
	}
      for(j = 0; j < M * N; j++)
	{
	  SType1B[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2B[j] = ConvertType(SType1B[j], convertTo);
	  SType1C[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2C[j] = ConvertType(SType1C[j], convertTo);
	}
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType1beta = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      SType2beta = ConvertType(SType1beta, convertTo);
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  cout << "SType1alpha = " << SType1alpha << endl;
	  cout << "SType2alpha = " << SType2alpha << endl;
	  cout << "SType1beta = " << SType1beta << endl;
	  cout << "SType2beta = " << SType2beta << endl;
	  if(SIDE == 'L')
	    {
	      PrintMatrix(SType1A, M, M, "SType1A", matlab);
	      PrintMatrix(SType2A, M, M, "SType2A", matlab);
	    }
	  else
	    {
	      PrintMatrix(SType1A, N, N, "SType1A", matlab);
	      PrintMatrix(SType2A, N, N, "SType2A", matlab);
	    }
	  PrintMatrix(SType1B, M, N, "SType1B", matlab);
	  PrintMatrix(SType1C, M, N, "SType1C_before_operation", matlab);
	  PrintMatrix(SType2B, M, N, "SType2B", matlab);
	  PrintMatrix(SType2C, M, N, "SType2C_before_operation", matlab);
	}
      TotalTestCount++;
      UPLO = RandomUPLO();
      SType1BLAS.SYMM(SIDE, UPLO, M, N, SType1alpha, SType1A, LDA, SType1B, M, SType1beta, SType1C, M);
      SType2BLAS.SYMM(SIDE, UPLO, M, N, SType2alpha, SType2A, LDA, SType2B, M, SType2beta, SType2C, M);
      if(debug)
	{
	  PrintMatrix(SType1C, M, N, "SType1C_after_operation", matlab);
	  PrintMatrix(SType2C, M, N, "SType2C_after_operation", matlab);
	}
      GoodTestSubcount += CompareMatrices(SType1C, SType2C, M, N, TOL);
      delete [] SType1A;
      delete [] SType1B;
      delete [] SType1C;
      delete [] SType2A;
      delete [] SType2B;
      delete [] SType2C;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) cout << "SYMM: " << GoodTestSubcount << " of " << SYMMTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End SYMM Tests

  // Begin TRMM Tests
  GoodTestSubcount = 0;
  for(i = 0; i < TRMMTESTS; i++)
    { 
      M = GetRandom(MVMIN, MVMAX);
      N = GetRandom(MVMIN, MVMAX);
      SType1B = new SType1[M * N];
      SType2B = new SType2[M * N];

      SIDE = RandomSIDE();
      if(SIDE == 'L')
	{
	  SType1A = new SType1[M * M];
	  SType2A = new SType2[M * M];
	  for(j = 0; j < M * M; j++)
	    {
	      SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	      SType2A[j] = ConvertType(SType1A[j], convertTo);
	    }
	  LDA = M;
	}
      else
	{
	  SType1A = new SType1[N * N];
	  SType2A = new SType2[N * N];
	  for(j = 0; j < N * N; j++)
	    {
	      SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	      SType2A[j] = ConvertType(SType1A[j], convertTo);
	    }
	  LDA = N;
	}
      for(j = 0; j < M * N; j++)
	{
	  SType1B[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2B[j] = ConvertType(SType1B[j], convertTo);
	}
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      UPLO = RandomUPLO();
      TRANSA = RandomTRANS();
      DIAG = RandomDIAG();
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  cout << "SType1alpha = " << SType1alpha << endl;
	  cout << "SType2alpha = " << SType2alpha << endl;
	  PrintMatrix(SType1A, LDA, LDA, "SType1A", matlab);
	  PrintMatrix(SType2A, LDA, LDA, "SType2A", matlab);
	  PrintMatrix(SType1B, M, N, "SType1B_before_operation", matlab);
	  PrintMatrix(SType2B, M, N, "SType2B_before_operation", matlab);
	}
      TotalTestCount++;
      SType1BLAS.TRMM(SIDE, UPLO, TRANSA, DIAG, M, N, SType1alpha, SType1A, LDA, SType1B, M);
      SType2BLAS.TRMM(SIDE, UPLO, TRANSA, DIAG, M, N, SType2alpha, SType2A, LDA, SType2B, M);
      if(debug)
	{
	  PrintMatrix(SType1B, M, N, "SType1B_after_operation", matlab);
	  PrintMatrix(SType2B, M, N, "SType2B_after_operation", matlab);
	}
      GoodTestSubcount += CompareMatrices(SType1B, SType2B, M, N, TOL);
      delete [] SType1A;
      delete [] SType1B;
      delete [] SType2A;
      delete [] SType2B;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) cout << "TRMM: " << GoodTestSubcount << " of " << TRMMTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End TRMM Tests

  // Begin TRSM Tests

  GoodTestSubcount = 0;
  for(i = 0; i < TRSMTESTS; i++)
    { 
      M = GetRandom(MVMIN, MVMAX);
      N = GetRandom(MVMIN, MVMAX);

      M = 3;
      N = 3;

      SType1B = new SType1[M * N];
      SType2B = new SType2[M * N];

      SIDE = RandomSIDE();

      SIDE = 'R';

      if(SIDE == 'L')
	{
	  SType1A = new SType1[M * M];
	  SType2A = new SType2[M * M];
	  for(j = 0; j < M * M; j++)
	    {
	     
	      SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	      SType2A[j] = ConvertType(SType1A[j], convertTo);
	    }
	  LDA = M;
	}
      else
	{
	  SType1A = new SType1[N * N];
	  SType2A = new SType2[N * N];
	  for(j = 0; j < N * N; j++)
	    {
	      SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	      SType2A[j] = ConvertType(SType1A[j], convertTo);
	    }
	  LDA = N;
	}

      for(j = 0; j < M * N; j++)
	{
	  SType1B[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2B[j] = ConvertType(SType1B[j], convertTo);
	}
      
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType1alpha = 1;
      SType2alpha = ConvertType(SType1alpha, convertTo);
 
      UPLO = RandomUPLO();
      TRANSA = RandomTRANS();
      DIAG = RandomDIAG();

      UPLO = 'U';
      TRANSA = 'N';
      DIAG = 'N';
 
      if(debug)
	{
	  cout << "Test #" << TotalTestCount << " --" << endl;
	  cout << "SType1alpha = " << SType1alpha << endl;
	  cout << "SType2alpha = " << SType2alpha << endl;
	  PrintMatrix(SType1A, LDA, LDA, "SType1A", matlab);
	  PrintMatrix(SType2A, LDA, LDA, "SType2A", matlab);
	  PrintMatrix(SType1B, M, N, "SType1B_before_operation", matlab);
	  PrintMatrix(SType2B, M, N, "SType2B_before_operation", matlab);
	}
      TotalTestCount++;

      SType1BLAS.TRSM(SIDE, UPLO, TRANSA, DIAG, M, N, SType1alpha, SType1A, LDA, SType1B, M);
      SType2BLAS.TRSM(SIDE, UPLO, TRANSA, DIAG, M, N, SType2alpha, SType2A, LDA, SType2B, M);
 
      if(debug)
	{
	  PrintMatrix(SType1B, M, N, "SType1B_after_operation", matlab);
	  PrintMatrix(SType2B, M, N, "SType2B_after_operation", matlab);
	}

      GoodTestSubcount += CompareMatrices(SType1B, SType2B, M, N, TOL);

      delete [] SType1A;
      delete [] SType1B;
      delete [] SType2A;
      delete [] SType2B;
    }
  GoodTestCount += GoodTestSubcount; 
  if(verbose || debug) cout << "TRSM: " << GoodTestSubcount << " of " << TRSMTESTS << " tests were successful." << endl;
  if(debug) cout << endl;
  // End TRSM Tests
  
  mp::mp_finalize();

  if((((TotalTestCount - 1) - GoodTestCount) != 0) || (verbose) || (debug))
    {
      cout << GoodTestCount << " of " << (TotalTestCount - 1) << " total tests were successful." << endl;
    }

  return 0;
}

template<typename TYPE>
TYPE GetRandom(TYPE Low, TYPE High)
{
  return ScalarTraits<TYPE>::random();
}

template<>
int GetRandom(int Low, int High)
{
  return ((int)((double)((1.0 * ScalarTraits<int>::random()) / RAND_MAX) * (High - Low + 1)) + Low);
}

template<>
double GetRandom(double Low, double High)
{
  return (((double)((1.0 * ScalarTraits<int>::random()) / RAND_MAX) * (High - Low + 1)) + Low + ScalarTraits<double>::random());
}

template<typename TYPE>
void PrintVector(TYPE* Vector, int Size, string Name, bool Matlab = 0)
{
  cout << Name << " =" << endl;
  int i;
  if(Matlab) cout << "[";
  for(i = 0; i < Size; i++)
    {
      cout << Vector[i] << " ";
    }
  if(Matlab) cout << "]";
  if(!Matlab)
    {
      cout << endl << endl;
    }
  else
    {
      cout << ";" << endl;
    }
}

template<typename TYPE>
void PrintMatrix(TYPE* Matrix, int Rows, int Columns, string Name, bool Matlab = 0)
{
  if(!Matlab)
    {
      cout << Name << " =" << endl;
      int i, j;
      for(j = 0; j < Rows; j++)
	{
	  for(i = 0; i < Columns; i++)
	    {
	      cout << Matrix[j + (i * Rows)] << " ";
	    }
	  cout << endl;
	}
      cout << endl;
    }
  else
    {
      cout << Name << " = ";
      int i, j;
      cout << "[";
      for(i = 0; i < Rows; i++)
	{
	  cout << "[";
	  for(j = 0; j < Columns; j++)
	    {
	      cout << Matrix[i + (j * Rows)] << " ";
	    }
	  cout << "];";
	}
      cout << "];" << endl;
    }
}

template<typename TYPE1, typename TYPE2>
bool CompareScalars(TYPE1 Scalar1, TYPE2 Scalar2, double Tolerance = 0)
{
  TYPE2 convertTo;
  return(ScalarTraits<TYPE2>::magnitude(ScalarTraits<TYPE2>::magnitude(ConvertType(Scalar1, convertTo)) - ScalarTraits<TYPE2>::magnitude(Scalar2)) <= Tolerance);
}

template<typename TYPE1, typename TYPE2>
bool CompareVectors(TYPE1* Vector1, TYPE2* Vector2, int Size, double Tolerance = 0)
{
  TYPE2 convertTo;
  int i;
  for(i = 0; i < Size; i++)
    {
      // if(ScalarTraits<TYPE1>::magnitude(ScalarTraits<TYPE1>::magnitude(Vector1[i]) - ScalarTraits<TYPE1>::magnitude((TYPE1)Vector2[i])) > Tolerance)
      if(ScalarTraits<TYPE2>::magnitude(ScalarTraits<TYPE2>::magnitude(ConvertType(Vector1[i], convertTo)) - ScalarTraits<TYPE2>::magnitude(Vector2[i])) > Tolerance)
	{
	  return 0;
	}
    }
  return 1;
}

template<typename TYPE1, typename TYPE2>
bool CompareMatrices(TYPE1* Matrix1, TYPE2* Matrix2, int Rows, int Columns, double Tolerance = 0)
{
  TYPE2 convertTo;
  int i;
  for(i = 0; i < (Rows * Columns); i++)
    {
      if(ScalarTraits<TYPE2>::magnitude(ScalarTraits<TYPE2>::magnitude(ConvertType(Matrix1[i], convertTo)) - ScalarTraits<TYPE2>::magnitude(Matrix2[i])) > Tolerance)
	{
	  return 0;
	}
    }
  return 1;
}

template<typename TYPE1, typename TYPE2>
TYPE2 ConvertType(TYPE1 T1, TYPE2 T2)
{
  return static_cast<TYPE2>(T1);
}

template<>
double ConvertType(mp_real T1, double T2)
{
  return dble(T1);
}

char RandomSIDE()
{
  char result = 'X';
  int r = GetRandom(1, 2);
  if(r == 1)
    {
      result = 'L';
    }
  else
    {
      result = 'R';
    }
  return result;
}

char RandomUPLO()
{
  char result = 'X';
  int r = GetRandom(1, 2);
  if(r == 1)
    {
      result = 'U';
    }
  else
    {
      result = 'L';
    }
  return result;
}

char RandomTRANS()
{
  char result = 'X';
  int r = GetRandom(1, 4);
  if(r == 1 || r == 2)
    {
      result = 'N';
    }
  else if(r == 3)
    {
      result = 'T';
    }
  else
    {
      result = 'C';
    }
  return result;
}

char RandomDIAG()
{
  char result = 'X';
  int r = GetRandom(1, 2);
  if(r == 1)
    {
      result = 'N';
    }
  else
    {
      result = 'U';
    }
  return result;
}

