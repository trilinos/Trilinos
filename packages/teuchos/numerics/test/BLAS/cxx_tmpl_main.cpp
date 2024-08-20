// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Kris
// 07.24.03 -- Initial checkin
// 08.08.03 -- All test suites except for TRSM are finished.
// 08.14.03 -- The test suite for TRSM is finished (Heidi).
/*

This test program is intended to check an experimental default type (e.g. mp_real) against an "officialy supported" control type (e.g. double). For each test, the program will generate the appropriate scalars and randomly-sized vectors and matrices of random numbers for the routine being tested. All of the input data for the experimental type is casted into the control type, so in theory both BLAS routines should get the same input data. Upon return, the experimental output data is casted back into the control type, and the results are compared; if they are equal (within a user-definable tolerance) the test is considered successful.

The test routine for TRSM is still being developed; all of the others are more or less finalized.

*/

#include "Teuchos_BLAS.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using Teuchos::BLAS;
using Teuchos::ScalarTraits;

// SType1 and SType2 define the datatypes for which BLAS output will be compared.
// SType2 should generally be a control datatype "officially" supported by the BLAS; SType1 should be the experimental type being checked.

// Define the first scalar type
#ifdef HAVE_TEUCHOS_COMPLEX
#define SType1     std::complex<float>
#else
#define SType1     float
#endif

// Define the second scalar type
#ifdef HAVE_TEUCHOS_COMPLEX
#define SType2     std::complex<double>
#else
#define SType2     double
#endif

// Define the ordinal type
#define OType	   long int

// MVMIN/MAX define the minimum and maximum dimensions of generated matrices and vectors, respectively.
#define MVMIN      2
#define MVMAX      20
// SCALARMAX defines the maximum positive value (with a little leeway) generated for matrix and std::vector elements and scalars:
// random numbers in [-SCALARMAX, SCALARMAX] will be generated.
// Set SCALARMAX to a floating-point value (e.g. 10.0) to enable floating-point random number generation, such that
// random numbers in (-SCALARMAX - 1, SCALARMAX + 1) will be generated.
// Large values of SCALARMAX may cause problems with SType2 = int, as large integer values will overflow floating-point types.
#define SCALARMAX  10
// These define the number of tests to be run for each individual BLAS routine.
#define ROTGTESTS  5
#define ROTTESTS   5
#define ASUMTESTS  5
#define AXPYTESTS  5
#define COPYTESTS  5
#define DOTTESTS   5
#define IAMAXTESTS 5
#define NRM2TESTS  5
#define SCALTESTS  5
#define GEMVTESTS  5
#define GERTESTS   5
#define TRMVTESTS  5
#define GEMMTESTS  5
#define SYMMTESTS  5
#define SYRKTESTS  5
#define TRMMTESTS  5
#define TRSMTESTS  5

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
void PrintVector(TYPE* Vector, OType Size, std::string Name, bool Matlab = 0);

template<typename TYPE>
void PrintMatrix(TYPE* Matrix, OType Rows, OType Columns, OType LDM, std::string Name, bool Matlab = 0);

template<typename TYPE1, typename TYPE2>
bool CompareScalars(TYPE1 Scalar1, TYPE2 Scalar2, typename ScalarTraits<TYPE2>::magnitudeType Tolerance );

template<typename TYPE1, typename TYPE2>
bool CompareVectors(TYPE1* Vector1, TYPE2* Vector2, OType Size, typename ScalarTraits<TYPE2>::magnitudeType Tolerance );

template<typename TYPE1, typename TYPE2>
bool CompareMatrices(TYPE1* Matrix1, TYPE2* Matrix2, OType Rows, OType Columns, OType LDM, typename ScalarTraits<TYPE2>::magnitudeType Tolerance );

// For most types, this function is just a wrapper for static_cast(), but for mp_real/double, it calls mp::dble()
// The second input parameter is not used; it is only needed to determine what type to convert *to*
template<typename TYPE1, typename TYPE2>
TYPE2 ConvertType(TYPE1, TYPE2);

// These functions return a random character appropriate for the BLAS arguments that share their names (uses GetRandom())
Teuchos::ESide RandomSIDE();
Teuchos::EUplo RandomUPLO();
Teuchos::ETransp RandomTRANS();
Teuchos::EDiag RandomDIAG();

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  bool verbose = 0;
  bool debug = 0;
  bool matlab = 0;
  bool InvalidCmdLineArgs = 0;
  OType i, j, k;
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

  if (verbose)
    std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  if(InvalidCmdLineArgs || (argc > 4))
    {
      std::cout << "Invalid command line arguments detected. Use the following flags:" << std::endl
	   << "\t -v enables verbose mode (reports number of failed/successful tests)" << std::endl
	   << "\t -d enables debug mode (same as verbose with output of each test, not recommended for large numbers of tests)" << std::endl
	   << "\t -m enables matlab-style output; only has an effect if debug mode is enabled" << std::endl;
      return 1;
    }
  typedef ScalarTraits<SType1>::magnitudeType MType1;
  typedef ScalarTraits<SType2>::magnitudeType MType2;
  BLAS<OType, SType1> SType1BLAS;
  BLAS<OType, SType2> SType2BLAS;
  SType1 SType1zero = ScalarTraits<SType1>::zero();
  SType1 SType1one = ScalarTraits<SType1>::one();
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
  SType1 SType1ASUMresult, SType1DOTresult, SType1NRM2result, SType1SINresult;
  SType2 SType2ASUMresult, SType2DOTresult, SType2NRM2result, SType2SINresult;
  MType1 SType1COSresult;
  MType2 SType2COSresult;
  OType incx, incy;
  OType SType1IAMAXresult;
  OType SType2IAMAXresult;
  OType TotalTestCount = 1, GoodTestSubcount, GoodTestCount = 0, M, M2, N, N2, P, K, LDA, LDB, LDC, Mx, My;
  Teuchos::EUplo UPLO;
  Teuchos::ESide SIDE;
  Teuchos::ETransp TRANS, TRANSA, TRANSB;
  Teuchos::EDiag DIAG;
  SType2 convertTo = ScalarTraits<SType2>::zero();
  MType2 mConvertTo = ScalarTraits<MType2>::zero();
  MType2 TOL = 1e-5*ScalarTraits<MType2>::one();

  std::srand(time(NULL));

  //--------------------------------------------------------------------------------
  // BEGIN LEVEL 1 BLAS TESTS
  //--------------------------------------------------------------------------------
  // Begin ROTG Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < ROTGTESTS; i++)
    {
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      SType1beta = GetRandom(-SCALARMAX, SCALARMAX);
      SType2beta = ConvertType(SType1beta, convertTo);
      SType1COSresult = ScalarTraits<MType1>::zero();
      SType2COSresult = ConvertType(SType1COSresult, ScalarTraits<MType2>::zero());
      SType1SINresult = ScalarTraits<SType1>::zero();
      SType2SINresult = ConvertType(SType1SINresult, convertTo);

      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  std::cout << "SType1alpha = "  << SType1alpha << std::endl;
	  std::cout << "SType2alpha = " << SType2alpha << std::endl;
	  std::cout << "SType1beta = "  << SType1beta << std::endl;
	  std::cout << "SType2beta = " << SType2beta << std::endl;
	}
      TotalTestCount++;
      SType1BLAS.ROTG(&SType1alpha, &SType1beta, &SType1COSresult, &SType1SINresult);
      SType2BLAS.ROTG(&SType2alpha, &SType2beta, &SType2COSresult, &SType2SINresult);
      if(debug)
	{
	  std::cout << "SType1 ROTG COS result: " << SType1COSresult << std::endl;
	  std::cout << "SType2 ROTG COS result: " << SType2COSresult << std::endl;
	  std::cout << "SType1 ROTG SIN result: " << SType1SINresult << std::endl;
	  std::cout << "SType2 ROTG SIN result: " << SType2SINresult << std::endl;
	}
      GoodTestSubcount += ( CompareScalars(SType1COSresult, SType2COSresult, TOL) &&
			    CompareScalars(SType1SINresult, SType2SINresult, TOL) );
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "ROTG: " << GoodTestSubcount << " of " << ROTGTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End ROTG Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin ROT Tests
  //--------------------------------------------------------------------------------
  typedef Teuchos::ScalarTraits<SType1>::magnitudeType MType1;
  typedef Teuchos::ScalarTraits<SType2>::magnitudeType MType2;
  GoodTestSubcount = 0;
  for(i = 0; i < ROTTESTS; i++)
  {
    incx = GetRandom(-5,5);
    incy = GetRandom(-5,5);
    if (incx == 0) incx = 1;
    if (incy == 0) incy = 1;
    M = GetRandom(MVMIN, MVMIN+8);
    Mx = M*std::abs(incx);
    My = M*std::abs(incy);
    if (Mx == 0) { Mx = 1; }
    if (My == 0) { My = 1; }
    SType1x = new SType1[Mx];
    SType1y = new SType1[My];
    SType2x = new SType2[Mx];
    SType2y = new SType2[My];
    for(j = 0; j < Mx; j++)
    {
      SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
      SType2x[j] = ConvertType(SType1x[j], convertTo);
    }
    for(j = 0; j < My; j++)
    {
      SType1y[j] = GetRandom(-SCALARMAX, SCALARMAX);
      SType2y[j] = ConvertType(SType1y[j], convertTo);
    }
    MType1 c1 = Teuchos::ScalarTraits<SType1>::magnitude(cos(static_cast<double>(GetRandom(-SCALARMAX,SCALARMAX))));
    MType2 c2 = ConvertType(c1, mConvertTo);
    SType1 s1 = sin(static_cast<double>(GetRandom(-SCALARMAX,SCALARMAX)));
    SType2 s2 = ConvertType(s1, convertTo);
    if(debug)
    {
      std::cout << "Test #" << TotalTestCount << " --" << std::endl;
      std::cout << "c1 = "  << c1 << ", s1 = " << s1 << std::endl;
      std::cout << "c2 = " << c2 << ", s2 = " << s2 << std::endl;
      std::cout << "incx = " << incx << ", incy = " << incy << std::endl;
      PrintVector(SType1x, Mx, "SType1x", matlab);
      PrintVector(SType1y, My, "SType1y_before_operation", matlab);
      PrintVector(SType2x, Mx, "SType2x", matlab);
      PrintVector(SType2y, My, "SType2y_before_operation",  matlab);
    }
    TotalTestCount++;
    SType1BLAS.ROT(M, SType1x, incx, SType1y, incy, &c1, &s1);
    SType2BLAS.ROT(M, SType2x, incx, SType2y, incy, &c2, &s2);
    if(debug)
    {
      PrintVector(SType1y, My, "SType1y_after_operation", matlab);
      PrintVector(SType2y, My, "SType2y_after_operation", matlab);
    }
    GoodTestSubcount += ( CompareVectors(SType1x, SType2x, Mx, TOL) &&
                          CompareVectors(SType1y, SType2y, My, TOL) );
    delete [] SType1x;
    delete [] SType1y;
    delete [] SType2x;
    delete [] SType2y;
  }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "ROT: " << GoodTestSubcount << " of " << ROTTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End ROT Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin ASUM Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  ScalarTraits<int>::seedrandom(0);
  for(i = 0; i < ASUMTESTS; i++)
    {
      incx = GetRandom(1, SCALARMAX);
      M = GetRandom(MVMIN, MVMAX);
      M2 = M*incx;
      SType1x = new SType1[M2];
      SType2x = new SType2[M2];
      for(j = 0; j < M2; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  PrintVector(SType1x, M2, "SType1x", matlab);
	  PrintVector(SType2x, M2, "SType2x", matlab);
	}
      TotalTestCount++;
      SType1ASUMresult = SType1BLAS.ASUM(M, SType1x, incx);
      SType2ASUMresult = SType2BLAS.ASUM(M, SType2x, incx);
      if(debug)
	{
	  std::cout << "SType1 ASUM result: " << SType1ASUMresult << std::endl;
	  std::cout << "SType2 ASUM result: " << SType2ASUMresult << std::endl;
	}
      GoodTestSubcount += CompareScalars(SType1ASUMresult, SType2ASUMresult, TOL);
      delete [] SType1x;
      delete [] SType2x;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "ASUM: " << GoodTestSubcount << " of " << ASUMTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;

  //--------------------------------------------------------------------------------
  // End ASUM Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin AXPY Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < AXPYTESTS; i++)
    {
      incx = GetRandom(1, MVMAX);
      incy = GetRandom(1, MVMAX);
      M = GetRandom(MVMIN, MVMAX);
      Mx = M*std::abs(incx);
      My = M*std::abs(incy);
      if (Mx == 0) { Mx = 1; }
      if (My == 0) { My = 1; }
      SType1x = new SType1[Mx];
      SType1y = new SType1[My];
      SType2x = new SType2[Mx];
      SType2y = new SType2[My];
      for(j = 0; j < Mx; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      for(j = 0; j < My; j++)
	{
	  SType1y[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2y[j] = ConvertType(SType1y[j], convertTo);
	}
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  std::cout << "SType1alpha = "  << SType1alpha << std::endl;
	  std::cout << "SType2alpha = " << SType2alpha << std::endl;
	  PrintVector(SType1x, Mx, "SType1x", matlab);
	  PrintVector(SType1y, My, "SType1y_before_operation", matlab);
	  PrintVector(SType2x, Mx, "SType2x", matlab);
	  PrintVector(SType2y, My, "SType2y_before_operation",  matlab);
	}
      TotalTestCount++;
      SType1BLAS.AXPY(M, SType1alpha, SType1x, incx, SType1y, incy);
      SType2BLAS.AXPY(M, SType2alpha, SType2x, incx, SType2y, incy);
      if(debug)
	{
	  PrintVector(SType1y, My, "SType1y_after_operation", matlab);
	  PrintVector(SType2y, My, "SType2y_after_operation", matlab);
	}
      GoodTestSubcount += CompareVectors(SType1y, SType2y, My, TOL);
      delete [] SType1x;
      delete [] SType1y;
      delete [] SType2x;
      delete [] SType2y;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "AXPY: " << GoodTestSubcount << " of " << AXPYTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End AXPY Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin COPY Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < COPYTESTS; i++)
    {
      incx = GetRandom(1, MVMAX);
      incy = GetRandom(1, MVMAX);
      M = GetRandom(MVMIN, MVMAX);
      Mx = M*std::abs(incx);
      My = M*std::abs(incy);
      if (Mx == 0) { Mx = 1; }
      if (My == 0) { My = 1; }
      SType1x = new SType1[Mx];
      SType1y = new SType1[My];
      SType2x = new SType2[Mx];
      SType2y = new SType2[My];
      for(j = 0; j < Mx; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      for(j = 0; j < My; j++)
	{
	  SType1y[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2y[j] = ConvertType(SType1y[j], convertTo);
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  PrintVector(SType1x, Mx, "SType1x", matlab);
	  PrintVector(SType1y, My, "SType1y_before_operation", matlab);
	  PrintVector(SType2x, Mx, "SType2x", matlab);
	  PrintVector(SType2y, My, "SType2y_before_operation", matlab);
	}
      TotalTestCount++;
      SType1BLAS.COPY(M, SType1x, incx, SType1y, incy);
      SType2BLAS.COPY(M, SType2x, incx, SType2y, incy);
      if(debug)
	{
	  PrintVector(SType1y, My, "SType1y_after_operation", matlab);
	  PrintVector(SType2y, My, "SType2y_after_operation", matlab);
	}
      GoodTestSubcount += CompareVectors(SType1y, SType2y, My, TOL);
      delete [] SType1x;
      delete [] SType1y;
      delete [] SType2x;
      delete [] SType2y;
    }
   GoodTestCount += GoodTestSubcount; if(verbose || debug) std::cout << "COPY: " << GoodTestSubcount << " of " << COPYTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End COPY Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin DOT Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < DOTTESTS; i++)
    {
      incx = GetRandom(1, MVMAX);
      incy = GetRandom(1, MVMAX);
      M = GetRandom(MVMIN, MVMAX);
      Mx = M*std::abs(incx);
      My = M*std::abs(incy);
      if (Mx == 0) { Mx = 1; }
      if (My == 0) { My = 1; }
      SType1x = new SType1[Mx];
      SType1y = new SType1[My];
      SType2x = new SType2[Mx];
      SType2y = new SType2[My];
      for(j = 0; j < Mx; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      for(j = 0; j < My; j++)
	{
	  SType1y[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2y[j] = ConvertType(SType1y[j], convertTo);
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  PrintVector(SType1x, Mx, "SType1x", matlab);
	  PrintVector(SType1y, My, "SType1y", matlab);
	  PrintVector(SType2x, Mx, "SType2x", matlab);
	  PrintVector(SType2y, My, "SType2y", matlab);
	}
      TotalTestCount++;
      SType1DOTresult = SType1BLAS.DOT(M, SType1x, incx, SType1y, incy);
      SType2DOTresult = SType2BLAS.DOT(M, SType2x, incx, SType2y, incy);
      if(debug)
	{
	  std::cout << "SType1 DOT result: " << SType1DOTresult << std::endl;
	  std::cout << "SType2 DOT result: " << SType2DOTresult << std::endl;
	}
      GoodTestSubcount += CompareScalars(SType1DOTresult, SType2DOTresult, TOL);
      delete [] SType1x;
      delete [] SType1y;
      delete [] SType2x;
      delete [] SType2y;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "DOT: " << GoodTestSubcount << " of " << DOTTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End DOT Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin NRM2 Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < NRM2TESTS; i++)
    {
      incx = GetRandom(1, SCALARMAX);
      M = GetRandom(MVMIN, MVMAX);
      M2 = M*incx;
      SType1x = new SType1[M2];
      SType2x = new SType2[M2];
      for(j = 0; j < M2; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  PrintVector(SType1x, M2, "SType1x", matlab);
	  PrintVector(SType2x, M2, "SType2x", matlab);
	}
      TotalTestCount++;
      SType1NRM2result = SType1BLAS.NRM2(M, SType1x, incx);
      SType2NRM2result = SType2BLAS.NRM2(M, SType2x, incx);
      if(debug)
	{
	  std::cout << "SType1 NRM2 result: " << SType1NRM2result << std::endl;
	  std::cout << "SType2 NRM2 result: " << SType2NRM2result << std::endl;
	}
      GoodTestSubcount += CompareScalars(SType1NRM2result, SType2NRM2result, TOL);
      delete [] SType1x;
      delete [] SType2x;
    }
   GoodTestCount += GoodTestSubcount; if(verbose || debug) std::cout << "NRM2: " << GoodTestSubcount << " of " << NRM2TESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End NRM2 Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin SCAL Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < SCALTESTS; i++)
    {
      // These will only test for the case that the increment is > 0, the
      // templated case can handle when incx < 0, but the blas library doesn't
      // seem to be able to on some machines.
      incx = GetRandom(1, SCALARMAX);
      M = GetRandom(MVMIN, MVMAX);
      M2 = M*incx;
      SType1x = new SType1[M2];
      SType2x = new SType2[M2];
      for(j = 0; j < M2; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  std::cout << "SType1alpha = " << SType1alpha << std::endl;
	  std::cout << "SType2alpha = " << SType2alpha << std::endl;
	  PrintVector(SType1x, M2, "SType1x_before_operation", matlab);
	  PrintVector(SType2x, M2, "SType2x_before_operation", matlab);
	}
      TotalTestCount++;
      SType1BLAS.SCAL(M, SType1alpha, SType1x, incx);
      SType2BLAS.SCAL(M, SType2alpha, SType2x, incx);
      if(debug)
	{
	  PrintVector(SType1x, M2, "SType1x_after_operation", matlab);
	  PrintVector(SType2x, M2, "SType2x_after_operation", matlab);
	}
      GoodTestSubcount += CompareVectors(SType1x, SType2x, M2, TOL);
      delete [] SType1x;
      delete [] SType2x;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "SCAL: " << GoodTestSubcount << " of " << SCALTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End SCAL Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin IAMAX Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < IAMAXTESTS; i++)
    {
      incx = GetRandom(1, SCALARMAX);
      M = GetRandom(MVMIN, MVMAX);
      M2 = M*incx;
      SType1x = new SType1[M2];
      SType2x = new SType2[M2];
      for(j = 0; j < M2; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  PrintVector(SType1x, M2, "SType1x", matlab);
	  PrintVector(SType2x, M2, "SType2x", matlab);
	}
      TotalTestCount++;
      SType1IAMAXresult = SType1BLAS.IAMAX(M, SType1x, incx);
      SType2IAMAXresult = SType2BLAS.IAMAX(M, SType2x, incx);
      if(debug)
	{
	  std::cout << "SType1 IAMAX result: " << SType1IAMAXresult << std::endl;
	  std::cout << "SType2 IAMAX result: " << SType2IAMAXresult << std::endl;
	}
      GoodTestSubcount += (SType1IAMAXresult == SType2IAMAXresult);
      delete [] SType1x;
      delete [] SType2x;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "IAMAX: " << GoodTestSubcount << " of " << IAMAXTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End IAMAX Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // BEGIN LEVEL 2 BLAS TESTS
  //--------------------------------------------------------------------------------
  // Begin GEMV Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < GEMVTESTS; i++)
    {
      // The parameters used to construct the test problem are chosen to be
      // valid parameters, so the GEMV routine won't bomb out.
      incx = GetRandom(1, MVMAX);
      while (incx == 0) {
      	  incx = GetRandom(1, MVMAX);
      }
      incy = GetRandom(1, MVMAX);
      while (incy == 0) {
      	  incy = GetRandom(1, MVMAX);
      }
      M = GetRandom(MVMIN, MVMAX);
      N = GetRandom(MVMIN, MVMAX);

      TRANS = RandomTRANS();
      if (Teuchos::ETranspChar[TRANS] == 'N') {	
      	M2 = M*std::abs(incy);
      	N2 = N*std::abs(incx);
      } else {
	M2 = N*std::abs(incy);
	N2 = M*std::abs(incx);
      }

      LDA = GetRandom(MVMIN, MVMAX);
      while (LDA < M) {
          LDA = GetRandom(MVMIN, MVMAX);
      }

      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType1beta = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      SType2beta = ConvertType(SType1beta, convertTo);

      SType1A = new SType1[LDA * N];
      SType1x = new SType1[N2];
      SType1y = new SType1[M2];
      SType2A = new SType2[LDA * N];
      SType2x = new SType2[N2];
      SType2y = new SType2[M2];

      for(j = 0; j < LDA * N; j++)
	{
	  SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2A[j] = ConvertType(SType1A[j], convertTo);
	}
      for(j = 0; j < N2; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      for(j = 0; j < M2; j++)
	{
	  SType1y[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2y[j] = ConvertType(SType1y[j], convertTo);
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  std::cout << "SType1alpha = " << SType1alpha << std::endl;
	  std::cout << "SType2alpha = " << SType2alpha << std::endl;
	  std::cout << "SType1beta = " << SType1beta << std::endl;
	  std::cout << "SType2beta = " << SType2beta << std::endl;
	  PrintMatrix(SType1A, M, N, LDA, "SType1A", matlab);
	  PrintVector(SType1x, N2, "SType1x", matlab);
	  PrintVector(SType1y, M2, "SType1y_before_operation", matlab);
	  PrintMatrix(SType2A, M, N, LDA, "SType2A", matlab);
	  PrintVector(SType2x, N2, "SType2x", matlab);
	  PrintVector(SType2y, M2, "SType2y_before_operation", matlab);
	}
      TotalTestCount++;
      SType1BLAS.GEMV(TRANS, M, N, SType1alpha, SType1A, LDA, SType1x, incx, SType1beta, SType1y, incy);
      SType2BLAS.GEMV(TRANS, M, N, SType2alpha, SType2A, LDA, SType2x, incx, SType2beta, SType2y, incy);
      if(debug)
	{
	  PrintVector(SType1y, M2, "SType1y_after_operation", matlab);
	  PrintVector(SType2y, M2, "SType2y_after_operation", matlab);
	}
      GoodTestSubcount += CompareVectors(SType1y, SType2y, M2, TOL);
      delete [] SType1A;
      delete [] SType1x;
      delete [] SType1y;
      delete [] SType2A;
      delete [] SType2x;
      delete [] SType2y;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "GEMV: " << GoodTestSubcount << " of " << GEMVTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End GEMV Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin TRMV Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < TRMVTESTS; i++)
    {
      UPLO = RandomUPLO();
      TRANSA = RandomTRANS();

      // Since the entries are integers, we don't want to use the unit diagonal feature,
      // this creates ill-conditioned, nearly-singular matrices.
      //DIAG = RandomDIAG();
      DIAG = Teuchos::NON_UNIT_DIAG;

      N = GetRandom(MVMIN, MVMAX);
      incx = GetRandom(1, MVMAX);
      while (incx == 0) {
      	  incx = GetRandom(1, MVMAX);
      }
      N2 = N*std::abs(incx);
      SType1x = new SType1[N2];
      SType2x = new SType2[N2];

      for(j = 0; j < N2; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}

      LDA = GetRandom(MVMIN, MVMAX);
      while (LDA < N) {
	LDA = GetRandom(MVMIN, MVMAX);
      }
      SType1A = new SType1[LDA * N];
      SType2A = new SType2[LDA * N];

      for(j = 0; j < N; j++)
	{	
	  if(Teuchos::EUploChar[UPLO] == 'U') {
	    // The operator is upper triangular, make sure that the entries are
	    // only in the upper triangular part of A and the diagonal is non-zero.
	    for(k = 0; k < N; k++)
	    {
	      if(k < j) {
		SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
	      } else {
		SType1A[j*LDA+k] = SType1zero;
	      }
	      SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
	      if(k == j) {
		if (Teuchos::EDiagChar[DIAG] == 'N') {
		  SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		  while (SType1A[j*LDA+k] == SType1zero) {
		    SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		  }
		  SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		} else {
		  SType1A[j*LDA+k] = SType1one;
		  SType2A[j*LDA+k] = SType2one;
		}
	      }			
	    }
	  } else {
	    // The operator is lower triangular, make sure that the entries are
	    // only in the lower triangular part of A and the diagonal is non-zero.
	    for(k = 0; k < N; k++)
	      {
		if(k > j) {
		  SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		} else {
		  SType1A[j*LDA+k] = SType1zero;
		}
		SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		if(k == j) {
		  if (Teuchos::EDiagChar[DIAG] == 'N') {
		    SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		    while (SType1A[j*LDA+k] == SType1zero) {
		      SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		    }
		    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		  } else {
		    SType1A[j*LDA+k] = SType1one;
		    SType2A[j*LDA+k] = SType2one;
		  }
		}			
	      } // end for(k=0 ...		
	  } // end if(UPLO == 'U') ...
	} // end for(j=0 ...      for(j = 0; j < N*N; j++)

      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  PrintMatrix(SType1A, N, N, LDA,"SType1A", matlab);
	  PrintVector(SType1x, N2, "SType1x_before_operation", matlab);
	  PrintMatrix(SType2A, N, N, LDA, "SType2A", matlab);
	  PrintVector(SType2x, N2, "SType2x_before_operation", matlab);
	}
      TotalTestCount++;
      SType1BLAS.TRMV(UPLO, TRANSA, DIAG, N, SType1A, LDA, SType1x, incx);
      SType2BLAS.TRMV(UPLO, TRANSA, DIAG, N, SType2A, LDA, SType2x, incx);
      if(debug)
	{
	  PrintVector(SType1x, N2, "SType1x_after_operation", matlab);
	  PrintVector(SType2x, N2, "SType2x_after_operation", matlab);
	}
      GoodTestSubcount += CompareVectors(SType1x, SType2x, N2, TOL);
      delete [] SType1A;
      delete [] SType1x;
      delete [] SType2A;
      delete [] SType2x;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "TRMV: " << GoodTestSubcount << " of " << TRMVTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End TRMV Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin GER Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < GERTESTS; i++)
    {
      incx = GetRandom(1, MVMAX);
      while (incx == 0) {
      	  incx = GetRandom(1, MVMAX);
      }
      incy = GetRandom(1, MVMAX);
      while (incy == 0) {
      	  incy = GetRandom(1, MVMAX);
      }
      M = GetRandom(MVMIN, MVMAX);
      N = GetRandom(MVMIN, MVMAX);

      M2 = M*std::abs(incx);
      N2 = N*std::abs(incy);

      LDA = GetRandom(MVMIN, MVMAX);
      while (LDA < M) {
          LDA = GetRandom(MVMIN, MVMAX);
      }

      SType1A = new SType1[LDA * N];
      SType1x = new SType1[M2];
      SType1y = new SType1[N2];
      SType2A = new SType2[LDA * N];
      SType2x = new SType2[M2];
      SType2y = new SType2[N2];
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      for(j = 0; j < LDA * N; j++)
	{
	  SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2A[j] = ConvertType(SType1A[j], convertTo);
	}
      for(j = 0; j < M2; j++)
	{
	  SType1x[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2x[j] = ConvertType(SType1x[j], convertTo);
	}
      for(j = 0; j < N2; j++)
	{
	  SType1y[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2y[j] = ConvertType(SType1y[j], convertTo);
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  std::cout << "SType1alpha = " << SType1alpha << std::endl;
	  std::cout << "SType2alpha = " << SType2alpha << std::endl;
	  PrintMatrix(SType1A, M, N, LDA,"SType1A_before_operation", matlab);
	  PrintVector(SType1x, M2, "SType1x", matlab);
	  PrintVector(SType1y, N2, "SType1y", matlab);
	  PrintMatrix(SType2A, M, N, LDA,"SType2A_before_operation", matlab);
	  PrintVector(SType2x, M2, "SType2x", matlab);
	  PrintVector(SType2y, N2, "SType2y", matlab);
	}
      TotalTestCount++;
      SType1BLAS.GER(M, N, SType1alpha, SType1x, incx, SType1y, incy, SType1A, LDA);
      SType2BLAS.GER(M, N, SType2alpha, SType2x, incx, SType2y, incy, SType2A, LDA);
      if(debug)
	{
	  PrintMatrix(SType1A, M, N, LDA, "SType1A_after_operation", matlab);
	  PrintMatrix(SType2A, M, N, LDA, "SType2A_after_operation", matlab);
	}
      GoodTestSubcount += CompareMatrices(SType1A, SType2A, M, N, LDA, TOL);
      delete [] SType1A;
      delete [] SType1x;
      delete [] SType1y;
      delete [] SType2A;
      delete [] SType2x;
      delete [] SType2y;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "GER: " << GoodTestSubcount << " of " << GERTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End GER Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // BEGIN LEVEL 3 BLAS TESTS
  //--------------------------------------------------------------------------------
  // Begin GEMM Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < GEMMTESTS; i++)
    {
      TRANSA = RandomTRANS();
      TRANSB = RandomTRANS();
      M = GetRandom(MVMIN, MVMAX);
      N = GetRandom(MVMIN, MVMAX);
      P = GetRandom(MVMIN, MVMAX);

      if(debug)	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
      }
      LDA = GetRandom(MVMIN, MVMAX);
      if (Teuchos::ETranspChar[TRANSA] == 'N') {
	while (LDA < M) {  LDA = GetRandom(MVMIN, MVMAX); }
	SType1A = new SType1[LDA * P];
	SType2A = new SType2[LDA * P];
	for(j = 0; j < LDA * P; j++)
	{
	  SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2A[j] = ConvertType(SType1A[j], convertTo);
	}
	if (debug) {
	PrintMatrix(SType1A, M, P, LDA, "SType1A", matlab);
	PrintMatrix(SType2A, M, P, LDA, "SType2A", matlab);
	}
      } else {
	while (LDA < P) {  LDA = GetRandom(MVMIN, MVMAX); }
	SType1A = new SType1[LDA * M];
	SType2A = new SType2[LDA * M];
	for(j = 0; j < LDA * M; j++)
	{
	  SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2A[j] = ConvertType(SType1A[j], convertTo);
	}
	if (debug) {
	PrintMatrix(SType1A, P, M, LDA, "SType1A", matlab);
	PrintMatrix(SType2A, P, M, LDA, "SType2A", matlab);
	}
      }

      LDB = GetRandom(MVMIN, MVMAX);
      if (Teuchos::ETranspChar[TRANSB] == 'N') {
	while (LDB < P) {  LDB = GetRandom(MVMIN, MVMAX); }
	SType1B = new SType1[LDB * N];
	SType2B = new SType2[LDB * N];
	for(j = 0; j < LDB * N; j++)
	{
	  SType1B[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2B[j] = ConvertType(SType1B[j], convertTo);
	}
	if (debug) {
	  PrintMatrix(SType1B, P, N, LDB,"SType1B", matlab);
	  PrintMatrix(SType2B, P, N, LDB,"SType2B", matlab);
	}
      } else {
	while (LDB < N) {  LDB = GetRandom(MVMIN, MVMAX); }
	SType1B = new SType1[LDB * P];
	SType2B = new SType2[LDB * P];
	for(j = 0; j < LDB * P; j++)
	{
	  SType1B[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2B[j] = ConvertType(SType1B[j], convertTo);
	}	
	if (debug) {
	  PrintMatrix(SType1B, N, P, LDB,"SType1B", matlab);
	  PrintMatrix(SType2B, N, P, LDB,"SType2B", matlab);
	}
      }

      LDC = GetRandom(MVMIN, MVMAX);
      while (LDC < M) {  LDC = GetRandom(MVMIN, MVMAX); }
      SType1C = new SType1[LDC * N];
      SType2C = new SType2[LDC * N];
      for(j = 0; j < LDC * N; j++) {
	  SType1C[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2C[j] = ConvertType(SType1C[j], convertTo);
      }
      if(debug)
	{
	  PrintMatrix(SType1C, M, N, LDC, "SType1C_before_operation", matlab);
	  PrintMatrix(SType2C, M, N, LDC, "SType2C_before_operation", matlab);
	}
	
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType1beta = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      SType2beta = ConvertType(SType1beta, convertTo);

      TotalTestCount++;
      SType1BLAS.GEMM(TRANSA, TRANSB, M, N, P, SType1alpha, SType1A, LDA, SType1B, LDB, SType1beta, SType1C, LDC);
      SType2BLAS.GEMM(TRANSA, TRANSB, M, N, P, SType2alpha, SType2A, LDA, SType2B, LDB, SType2beta, SType2C, LDC);
      if(debug)
	{
	  PrintMatrix(SType1C, M, N, LDC, "SType1C_after_operation", matlab);
	  PrintMatrix(SType2C, M, N, LDC, "SType2C_after_operation", matlab);
	}
      GoodTestSubcount += CompareMatrices(SType1C, SType2C, M, N, LDC, TOL);
      delete [] SType1A;
      delete [] SType1B;
      delete [] SType1C;
      delete [] SType2A;
      delete [] SType2B;
      delete [] SType2C;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "GEMM: " << GoodTestSubcount << " of " << GEMMTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End GEMM Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin SYMM Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < SYMMTESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      N = GetRandom(MVMIN, MVMAX);
      SIDE = RandomSIDE();
      UPLO = RandomUPLO();

      LDA = GetRandom(MVMIN, MVMAX);
      if(Teuchos::ESideChar[SIDE] == 'L') {
	while (LDA < M) { LDA = GetRandom(MVMIN, MVMAX); }
	SType1A = new SType1[LDA * M];
	SType2A = new SType2[LDA * M];
	for(j = 0; j < LDA * M; j++) {
	  SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2A[j] = ConvertType(SType1A[j], convertTo);
	}
      } else {
	while (LDA < N) { LDA = GetRandom(MVMIN, MVMAX); }
	SType1A = new SType1[LDA * N];
	SType2A = new SType2[LDA * N];
	for(j = 0; j < LDA * N; j++) {
	  SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2A[j] = ConvertType(SType1A[j], convertTo);
	}
      }

      LDB = GetRandom(MVMIN, MVMAX);
      while (LDB < M) {  LDB = GetRandom(MVMIN, MVMAX); }
      SType1B = new SType1[LDB * N];
      SType2B = new SType2[LDB * N];
      for(j = 0; j < LDB * N; j++) {
	  SType1B[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2B[j] = ConvertType(SType1B[j], convertTo);
      }

      LDC = GetRandom(MVMIN, MVMAX);
      while (LDC < M) {  LDC = GetRandom(MVMIN, MVMAX); }
      SType1C = new SType1[LDC * N];
      SType2C = new SType2[LDC * N];
      for(j = 0; j < LDC * N; j++) {
	  SType1C[j] = GetRandom(-SCALARMAX, SCALARMAX);
	  SType2C[j] = ConvertType(SType1C[j], convertTo);
      }

      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType1beta = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      SType2beta = ConvertType(SType1beta, convertTo);
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  std::cout << "SType1alpha = " << SType1alpha << std::endl;
	  std::cout << "SType2alpha = " << SType2alpha << std::endl;
	  std::cout << "SType1beta = " << SType1beta << std::endl;
	  std::cout << "SType2beta = " << SType2beta << std::endl;
	  if (Teuchos::ESideChar[SIDE] == 'L') {
	      PrintMatrix(SType1A, M, M, LDA,"SType1A", matlab);
	      PrintMatrix(SType2A, M, M, LDA,"SType2A", matlab);
	  } else {
	    PrintMatrix(SType1A, N, N, LDA, "SType1A", matlab);
	    PrintMatrix(SType2A, N, N, LDA, "SType2A", matlab);
	  }
	  PrintMatrix(SType1B, M, N, LDB,"SType1B", matlab);
	  PrintMatrix(SType1C, M, N, LDC,"SType1C_before_operation", matlab);
	  PrintMatrix(SType2B, M, N, LDB,"SType2B", matlab);
	  PrintMatrix(SType2C, M, N, LDC,"SType2C_before_operation", matlab);
	}
      TotalTestCount++;

      SType1BLAS.SYMM(SIDE, UPLO, M, N, SType1alpha, SType1A, LDA, SType1B, LDB, SType1beta, SType1C, LDC);
      SType2BLAS.SYMM(SIDE, UPLO, M, N, SType2alpha, SType2A, LDA, SType2B, LDB, SType2beta, SType2C, LDC);
      if(debug)
	{
	  PrintMatrix(SType1C, M, N, LDC,"SType1C_after_operation", matlab);
	  PrintMatrix(SType2C, M, N, LDC,"SType2C_after_operation", matlab);
	}
      GoodTestSubcount += CompareMatrices(SType1C, SType2C, M, N, LDC, TOL);

      delete [] SType1A;
      delete [] SType1B;
      delete [] SType1C;
      delete [] SType2A;
      delete [] SType2B;
      delete [] SType2C;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "SYMM: " << GoodTestSubcount << " of " << SYMMTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End SYMM Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin SYRK Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < SYRKTESTS; i++)
  {
    N = GetRandom(MVMIN, MVMAX);
    K = GetRandom(MVMIN, MVMAX);
    while (K > N) { K = GetRandom(MVMIN, MVMAX); }

    UPLO = RandomUPLO();
    TRANS = RandomTRANS();
#ifdef HAVE_TEUCHOS_COMPLEX
    while (TRANS == Teuchos::CONJ_TRANS) { TRANS = RandomTRANS(); }
#endif

    LDA = GetRandom(MVMIN, MVMAX);
    if(Teuchos::ETranspChar[TRANS] == 'N') {
      while (LDA < N) { LDA = GetRandom(MVMIN, MVMAX); }
      SType1A = new SType1[LDA * K];
      SType2A = new SType2[LDA * K];
      for(j = 0; j < LDA * K; j++) {
        SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
        SType2A[j] = ConvertType(SType1A[j], convertTo);
      }
    } else {
      while (LDA < K) { LDA = GetRandom(MVMIN, MVMAX); }
      SType1A = new SType1[LDA * N];
      SType2A = new SType2[LDA * N];
      for(j = 0; j < LDA * N; j++) {
        SType1A[j] = GetRandom(-SCALARMAX, SCALARMAX);
        SType2A[j] = ConvertType(SType1A[j], convertTo);
      }
    }

    LDC = GetRandom(MVMIN, MVMAX);
    while (LDC < N) {  LDC = GetRandom(MVMIN, MVMAX); }
    SType1C = new SType1[LDC * N];
    SType2C = new SType2[LDC * N];
    for(j = 0; j < N; j++) {

      if(Teuchos::EUploChar[UPLO] == 'U') {
        // The operator is upper triangular, make sure that the entries are
        // only in the upper triangular part of C.
        for(k = 0; k < N; k++)
        {
          if(k <= j) {
            SType1C[j*LDC+k] = GetRandom(-SCALARMAX, SCALARMAX);
          } else {
            SType1C[j*LDC+k] = SType1zero;
          }
          SType2C[j*LDC+k] = ConvertType(SType1C[j*LDC+k], convertTo);
        }
      }
      else {
        for(k = 0; k < N; k++)
        {
          if(k >= j) {
            SType1C[j*LDC+k] = GetRandom(-SCALARMAX, SCALARMAX);
          } else {
            SType1C[j*LDC+k] = SType1zero;
          }
          SType2C[j*LDC+k] = ConvertType(SType1C[j*LDC+k], convertTo);
        }
      }
    }

    SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
    SType1beta = GetRandom(-SCALARMAX, SCALARMAX);
    SType2alpha = ConvertType(SType1alpha, convertTo);
    SType2beta = ConvertType(SType1beta, convertTo);
    if(debug)
    {
      std::cout << "Test #" << TotalTestCount << " --" << std::endl;
      std::cout << "SType1alpha = " << SType1alpha << std::endl;
      std::cout << "SType2alpha = " << SType2alpha << std::endl;
      std::cout << "SType1beta = " << SType1beta << std::endl;
      std::cout << "SType2beta = " << SType2beta << std::endl;
      if (Teuchos::ETranspChar[TRANS] == 'N') {
        PrintMatrix(SType1A, N, K, LDA,"SType1A", matlab);
        PrintMatrix(SType2A, N, K, LDA,"SType2A", matlab);
      } else {
        PrintMatrix(SType1A, K, N, LDA, "SType1A", matlab);
        PrintMatrix(SType2A, K, N, LDA, "SType2A", matlab);
      }
      PrintMatrix(SType1C, N, N, LDC,"SType1C_before_operation", matlab);
      PrintMatrix(SType2C, N, N, LDC,"SType2C_before_operation", matlab);
    }
    TotalTestCount++;

    SType1BLAS.SYRK(UPLO, TRANS, N, K, SType1alpha, SType1A, LDA, SType1beta, SType1C, LDC);
    SType2BLAS.SYRK(UPLO, TRANS, N, K, SType2alpha, SType2A, LDA, SType2beta, SType2C, LDC);
    if(debug)
    {
      PrintMatrix(SType1C, N, N, LDC,"SType1C_after_operation", matlab);
      PrintMatrix(SType2C, N, N, LDC,"SType2C_after_operation", matlab);
    }
    GoodTestSubcount += CompareMatrices(SType1C, SType2C, N, N, LDC, TOL);

    delete [] SType1A;
    delete [] SType1C;
    delete [] SType2A;
    delete [] SType2C;
  }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "SYRK: " << GoodTestSubcount << " of " << SYRKTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End SYRK Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin TRMM Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < TRMMTESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      N = GetRandom(MVMIN, MVMAX);

      LDB = GetRandom(MVMIN, MVMAX);
      while (LDB < M) {
	  LDB = GetRandom(MVMIN, MVMAX);
      }

      SType1B = new SType1[LDB * N];
      SType2B = new SType2[LDB * N];

      SIDE = RandomSIDE();
      UPLO = RandomUPLO();
      TRANSA = RandomTRANS();
      DIAG = RandomDIAG();

      if(Teuchos::ESideChar[SIDE] == 'L')  // The operator is on the left side
	{
          LDA = GetRandom(MVMIN, MVMAX);
      	  while (LDA < M) {
	      LDA = GetRandom(MVMIN, MVMAX);
       	  }

	  SType1A = new SType1[LDA * M];
	  SType2A = new SType2[LDA * M];

	  for(j = 0; j < M; j++)
	    {	
	      if(Teuchos::EUploChar[UPLO] == 'U') {
		// The operator is upper triangular, make sure that the entries are
		// only in the upper triangular part of A and the diagonal is non-zero.
		for(k = 0; k < M; k++)
		{
		    if(k < j) {
	      		SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			SType1A[j*LDA+k] = SType1zero;
		    }
	      	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    if(k == j) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
	      		    SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (SType1A[j*LDA+k] == SType1zero) {
				SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);	
		    	} else {
	      		    SType1A[j*LDA+k] = SType1one;
	      	    	    SType2A[j*LDA+k] = SType2one;
		    	}			
		    }
		}
	      } else {
		// The operator is lower triangular, make sure that the entries are
		// only in the lower triangular part of A and the diagonal is non-zero.
		for(k = 0; k < M; k++)
		{
		    if(k > j) {
	      		SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			SType1A[j*LDA+k] = SType1zero;
		    }
	      	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    if(k == j) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
      			    SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (SType1A[j*LDA+k] == SType1zero) {
				SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    	} else {
	      		    SType1A[j*LDA+k] = SType1one;
	      	    	    SType2A[j*LDA+k] = SType2one;
			}
		    }			
		} // end for(k=0 ...		
	      } // end if(UPLO == 'U') ...
	    } // end for(j=0 ...
	} // if(SIDE == 'L') ...
      else // The operator is on the right side
	{
          LDA = GetRandom(MVMIN, MVMAX);
      	  while (LDA < N) {
	      LDA = GetRandom(MVMIN, MVMAX);
       	  }

	  SType1A = new SType1[LDA * N];
	  SType2A = new SType2[LDA * N];

	  for(j = 0; j < N; j++)
	    {	
	      if(Teuchos::EUploChar[UPLO] == 'U') {
		// The operator is upper triangular, make sure that the entries are
		// only in the upper triangular part of A and the diagonal is non-zero.
		for(k = 0; k < N; k++)
		{
		    if(k < j) {
	      		SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			SType1A[j*LDA+k] = SType1zero;
		    }
	      	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    if(k == j) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
	      		    SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (SType1A[j*LDA+k] == SType1zero) {
				SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    	} else {
	      		    SType1A[j*LDA+k] = SType1one;
	      	    	    SType2A[j*LDA+k] = SType2one;
			}
		    }			
		}
	      } else {
		// The operator is lower triangular, make sure that the entries are
		// only in the lower triangular part of A and the diagonal is non-zero.
		for(k = 0; k < N; k++)
		{
		    if(k > j) {
	      		SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			SType1A[j*LDA+k] = SType1zero;
		    }
	      	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    if(k == j) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
	      		    SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (SType1A[j*LDA+k] == SType1zero) {
				SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    	} else {
	      		    SType1A[j*LDA+k] = SType1one;
	      	    	    SType2A[j*LDA+k] = SType2one;
			}
		    }			
		} // end for(k=0 ...		
	      } // end if(UPLO == 'U') ...
	    } // end for(j=0 ...
	} // end if(SIDE == 'L') ...

      // Fill in the right hand side block B.
      for(j = 0; j < N; j++) {
	  for(k = 0; k < M; k++) {
	    SType1B[j*LDB+k] = GetRandom(-SCALARMAX, SCALARMAX);
	    SType2B[j*LDB+k] = ConvertType(SType1B[j*LDB+k], convertTo);
	  }
      }
      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  std::cout << "SType1alpha = " << SType1alpha << std::endl;
	  std::cout << "SType2alpha = " << SType2alpha << std::endl;
          if(Teuchos::ESideChar[SIDE] == 'L') {
	    PrintMatrix(SType1A, M, M, LDA, "SType1A", matlab);
	    PrintMatrix(SType2A, M, M, LDA, "SType2A", matlab);
	  } else {
	    PrintMatrix(SType1A, N, N, LDA, "SType1A", matlab);
	    PrintMatrix(SType2A, N, N, LDA, "SType2A", matlab);
	  }
	  PrintMatrix(SType1B, M, N, LDB,"SType1B_before_operation", matlab);
	  PrintMatrix(SType2B, M, N, LDB,"SType2B_before_operation", matlab);
	}
      TotalTestCount++;
      SType1BLAS.TRMM(SIDE, UPLO, TRANSA, DIAG, M, N, SType1alpha, SType1A, LDA, SType1B, LDB);
      SType2BLAS.TRMM(SIDE, UPLO, TRANSA, DIAG, M, N, SType2alpha, SType2A, LDA, SType2B, LDB);
      if(debug)
	{
	  PrintMatrix(SType1B, M, N, LDB, "SType1B_after_operation", matlab);
	  PrintMatrix(SType2B, M, N, LDB, "SType2B_after_operation", matlab);
	}
      GoodTestSubcount += CompareMatrices(SType1B, SType2B, M, N, LDB, TOL);
      delete [] SType1A;
      delete [] SType1B;
      delete [] SType2A;
      delete [] SType2B;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "TRMM: " << GoodTestSubcount << " of " << TRMMTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End TRMM Tests
  //--------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------
  // Begin TRSM Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < TRSMTESTS; i++)
    {
      M = GetRandom(MVMIN, MVMAX);
      N = GetRandom(MVMIN, MVMAX);

      LDB = GetRandom(MVMIN, MVMAX);
      while (LDB < M) {
	  LDB = GetRandom(MVMIN, MVMAX);
      }

      SType1B = new SType1[LDB * N];
      SType2B = new SType2[LDB * N];

      SIDE = RandomSIDE();
      UPLO = RandomUPLO();
      TRANSA = RandomTRANS();
      // Since the entries are integers, we don't want to use the unit diagonal feature,
      // this creates ill-conditioned, nearly-singular matrices.
      //DIAG = RandomDIAG();
      DIAG = Teuchos::NON_UNIT_DIAG;

      if(Teuchos::ESideChar[SIDE] == 'L')  // The operator is on the left side
	{
          LDA = GetRandom(MVMIN, MVMAX);
      	  while (LDA < M) {
	      LDA = GetRandom(MVMIN, MVMAX);
       	  }

	  SType1A = new SType1[LDA * M];
	  SType2A = new SType2[LDA * M];

	  for(j = 0; j < M; j++)
	    {	
	      if(Teuchos::EUploChar[UPLO] == 'U') {
		// The operator is upper triangular, make sure that the entries are
		// only in the upper triangular part of A and the diagonal is non-zero.
		for(k = 0; k < M; k++)
		{
		    if(k < j) {
	      		SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			SType1A[j*LDA+k] = SType1zero;
		    }
	      	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    if(k == j) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
	      		    SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (SType1A[j*LDA+k] == SType1zero) {
				SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);	
		    	} else {
	      		    SType1A[j*LDA+k] = SType1one;
	      	    	    SType2A[j*LDA+k] = SType2one;
		    	}			
		    }
		}
	      } else {
		// The operator is lower triangular, make sure that the entries are
		// only in the lower triangular part of A and the diagonal is non-zero.
		for(k = 0; k < M; k++)
		{
		    if(k > j) {
	      		SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			SType1A[j*LDA+k] = SType1zero;
		    }
	      	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    if(k == j) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
      			    SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (SType1A[j*LDA+k] == SType1zero) {
				SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    	} else {
	      		    SType1A[j*LDA+k] = SType1one;
	      	    	    SType2A[j*LDA+k] = SType2one;
			}
		    }			
		} // end for(k=0 ...		
	      } // end if(UPLO == 'U') ...
	    } // end for(j=0 ...
	} // if(SIDE == 'L') ...
      else // The operator is on the right side
	{
          LDA = GetRandom(MVMIN, MVMAX);
      	  while (LDA < N) {
	      LDA = GetRandom(MVMIN, MVMAX);
       	  }

	  SType1A = new SType1[LDA * N];
	  SType2A = new SType2[LDA * N];

	  for(j = 0; j < N; j++)
	    {	
	      if(Teuchos::EUploChar[UPLO] == 'U') {
		// The operator is upper triangular, make sure that the entries are
		// only in the upper triangular part of A and the diagonal is non-zero.
		for(k = 0; k < N; k++)
		{
		    if(k < j) {
	      		SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			SType1A[j*LDA+k] = SType1zero;
		    }
	      	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    if(k == j) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
	      		    SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (SType1A[j*LDA+k] == SType1zero) {
				SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    	} else {
	      		    SType1A[j*LDA+k] = SType1one;
	      	    	    SType2A[j*LDA+k] = SType2one;
			}
		    }			
		}
	      } else {
		// The operator is lower triangular, make sure that the entries are
		// only in the lower triangular part of A and the diagonal is non-zero.
		for(k = 0; k < N; k++)
		{
		    if(k > j) {
	      		SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			SType1A[j*LDA+k] = SType1zero;
		    }
	      	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    if(k == j) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
	      		    SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (SType1A[j*LDA+k] == SType1zero) {
				SType1A[j*LDA+k] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    SType2A[j*LDA+k] = ConvertType(SType1A[j*LDA+k], convertTo);
		    	} else {
	      		    SType1A[j*LDA+k] = SType1one;
	      	    	    SType2A[j*LDA+k] = SType2one;
			}
		    }			
		} // end for(k=0 ...		
	      } // end if(UPLO == 'U') ...
	    } // end for(j=0 ...
	} // end if(SIDE == 'L') ...

      // Fill in the right hand side block B.
      for(j = 0; j < N; j++)
	{
	  for(k = 0; k < M; k++)
	    {
	  	SType1B[j*LDB+k] = GetRandom(-SCALARMAX, SCALARMAX);
	  	SType2B[j*LDB+k] = ConvertType(SType1B[j*LDB+k], convertTo);
	    }
	}

      SType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      SType2alpha = ConvertType(SType1alpha, convertTo);

      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " --" << std::endl;
	  std::cout << Teuchos::ESideChar[SIDE] << "\t"
	       << Teuchos::EUploChar[UPLO] << "\t"
	       << Teuchos::ETranspChar[TRANSA] << "\t"
	       << Teuchos::EDiagChar[DIAG] << std::endl;
	  std::cout << "M="<<M << "\t" << "N="<<N << "\t" << "LDA="<<LDA << "\t" << "LDB="<<LDB << std::endl;
	  std::cout << "SType1alpha = " << SType1alpha << std::endl;
	  std::cout << "SType2alpha = " << SType2alpha << std::endl;
	  if (Teuchos::ESideChar[SIDE] == 'L') {
	      PrintMatrix(SType1A, M, M, LDA, "SType1A", matlab);
	      PrintMatrix(SType2A, M, M, LDA, "SType2A", matlab);
	  } else {
	      PrintMatrix(SType1A, N, N, LDA, "SType1A", matlab);
	      PrintMatrix(SType2A, N, N, LDA, "SType2A", matlab);
	  }
	  PrintMatrix(SType1B, M, N, LDB, "SType1B_before_operation", matlab);
	  PrintMatrix(SType2B, M, N, LDB, "SType2B_before_operation", matlab);
	}
      TotalTestCount++;

      SType1BLAS.TRSM(SIDE, UPLO, TRANSA, DIAG, M, N, SType1alpha, SType1A, LDA, SType1B, LDB);
      SType2BLAS.TRSM(SIDE, UPLO, TRANSA, DIAG, M, N, SType2alpha, SType2A, LDA, SType2B, LDB);

      if(debug)
	{
	  PrintMatrix(SType1B, M, N, LDB, "SType1B_after_operation", matlab);
	  PrintMatrix(SType2B, M, N, LDB, "SType2B_after_operation", matlab);
	}

      if (CompareMatrices(SType1B, SType2B, M, N, LDB, TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareMatrices(SType1B, SType2B, M, N, LDB, TOL);

      delete [] SType1A;
      delete [] SType1B;
      delete [] SType2A;
      delete [] SType2B;
    }
  GoodTestCount += GoodTestSubcount;
  if(verbose || debug) std::cout << "TRSM: " << GoodTestSubcount << " of " << TRSMTESTS << " tests were successful." << std::endl;
  if(debug) std::cout << std::endl;
  //--------------------------------------------------------------------------------
  // End TRSM Tests
  //--------------------------------------------------------------------------------

  if((((TotalTestCount - 1) - GoodTestCount) != 0) || (verbose) || (debug))
    {
      std::cout << GoodTestCount << " of " << (TotalTestCount - 1) << " total tests were successful." << std::endl;
    }

  if ((TotalTestCount-1) == GoodTestCount) {
    std::cout << "End Result: TEST PASSED" << std::endl;
    return 0;
  }

  std::cout << "End Result: TEST FAILED" << std::endl;
  return (TotalTestCount-GoodTestCount-1);
}

template<typename TYPE>
TYPE GetRandom(TYPE Low, TYPE High)
{
  return ((TYPE)((double)((1.0 * ScalarTraits<int>::random()) / RAND_MAX) * (High - Low + 1)) + Low);
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
void PrintVector(TYPE* Vector, OType Size, std::string Name, bool Matlab)
{
  std::cout << Name << " =" << std::endl;
  OType i;
  if(Matlab) std::cout << "[";
  for(i = 0; i < Size; i++)
    {
      std::cout << Vector[i] << " ";
    }
  if(Matlab) std::cout << "]";
  if(!Matlab)
    {
      std::cout << std::endl << std::endl;
    }
  else
    {
      std::cout << ";" << std::endl;
    }
}

template<typename TYPE>
void PrintMatrix(TYPE* Matrix, OType Rows, OType Columns, OType LDM, std::string Name, bool Matlab)
{
  if(!Matlab)
    {
      std::cout << Name << " =" << std::endl;
      OType i, j;
      for(i = 0; i < Rows; i++)
	{
      	  for(j = 0; j < Columns; j++)
	    {
	      std::cout << Matrix[i + (j * LDM)] << " ";
	    }
	  std::cout << std::endl;
	}
      std::cout << std::endl;
    }
  else
    {
      std::cout << Name << " = ";
      OType i, j;
      std::cout << "[";
      for(i = 0; i < Rows; i++)
        {
	  std::cout << "[";
      	  for(j = 0; j < Columns; j++)
	    {
	      std::cout << Matrix[i + (j * LDM)] << " ";
	    }
	  std::cout << "];";
	}
      std::cout << "];" << std::endl;
    }
}

template<typename TYPE1, typename TYPE2>
bool CompareScalars(TYPE1 Scalar1, TYPE2 Scalar2, typename ScalarTraits<TYPE2>::magnitudeType Tolerance)
{
  TYPE2 convertTo = ScalarTraits<TYPE2>::zero();
  typename ScalarTraits<TYPE2>::magnitudeType temp = ScalarTraits<TYPE2>::magnitude(Scalar2);
  typename ScalarTraits<TYPE2>::magnitudeType temp2 = ScalarTraits<TYPE2>::magnitude(ConvertType(Scalar1, convertTo) - Scalar2);
  if (temp != ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::zero()) {
    temp2 /= temp;
  }
  return( temp2 < Tolerance );
}


/*  Function:  CompareVectors
    Purpose:   Compares the difference between two vectors using relative euclidean-norms, i.e. ||v_1-v_2||_2/||v_2||_2
*/
template<typename TYPE1, typename TYPE2>
bool CompareVectors(TYPE1* Vector1, TYPE2* Vector2, OType Size, typename ScalarTraits<TYPE2>::magnitudeType Tolerance)
{
  TYPE2 convertTo = ScalarTraits<TYPE2>::zero();
  TYPE2 temp = ScalarTraits<TYPE2>::zero();
  typename ScalarTraits<TYPE2>::magnitudeType temp2 = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::zero();
  typename ScalarTraits<TYPE2>::magnitudeType temp3 = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::zero();
  typename ScalarTraits<TYPE2>::magnitudeType sum = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::zero();
  typename ScalarTraits<TYPE2>::magnitudeType sum2 = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::zero();
  OType i;
  for(i = 0; i < Size; i++)
    {
      sum2 += ScalarTraits<TYPE2>::magnitude(ScalarTraits<TYPE2>::conjugate(Vector2[i])*Vector2[i]);
      temp = ConvertType(Vector1[i], convertTo) - Vector2[i];
      sum += ScalarTraits<TYPE2>::magnitude(ScalarTraits<TYPE2>::conjugate(temp)*temp);
    }
  temp2 = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::squareroot(sum2);
  if (temp2 != ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::zero())
    temp3 = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::squareroot(sum)/temp2;
  else
    temp3 = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::squareroot(sum);
  if (temp3 > Tolerance )
    return false;
  else
    return true;
}

/*  Function:  CompareMatrices
    Purpose:   Compares the difference between two matrices using relative frobenius-norms, i.e. ||M_1-M_2||_F/||M_2||_F
*/
template<typename TYPE1, typename TYPE2>
bool CompareMatrices(TYPE1* Matrix1, TYPE2* Matrix2, OType Rows, OType Columns, OType LDM, typename ScalarTraits<TYPE2>::magnitudeType Tolerance)
{
  TYPE2 convertTo = ScalarTraits<TYPE2>::zero();
  TYPE2 temp = ScalarTraits<TYPE2>::zero();
  typename ScalarTraits<TYPE2>::magnitudeType temp2 = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::zero();
  typename ScalarTraits<TYPE2>::magnitudeType temp3 = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::zero();
  typename ScalarTraits<TYPE2>::magnitudeType sum = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::zero();
  typename ScalarTraits<TYPE2>::magnitudeType sum2 = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::zero();
  OType i,j;
  for(j = 0; j < Columns; j++)
    {
      for(i = 0; i < Rows; i++)
	{
	  sum2 = ScalarTraits<TYPE2>::magnitude(ScalarTraits<TYPE2>::conjugate(Matrix2[j*LDM + i])*Matrix2[j*LDM + i]);
	  temp = ConvertType(Matrix1[j*LDM + i],convertTo) - Matrix2[j*LDM + i];
	  sum = ScalarTraits<TYPE2>::magnitude(ScalarTraits<TYPE2>::conjugate(temp)*temp);
	}
    }
  temp2 = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::squareroot(sum2);
  if (temp2 != ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::zero())
    temp3 = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::squareroot(sum)/temp2;
  else
    temp3 = ScalarTraits<typename ScalarTraits<TYPE2>::magnitudeType>::squareroot(sum);
  if (temp3 > Tolerance)
    return false;
  else
    return true;
}


template<typename TYPE1, typename TYPE2>
TYPE2 ConvertType(TYPE1 T1, TYPE2 T2)
{
  return static_cast<TYPE2>(T1);
}

Teuchos::ESide RandomSIDE()
{
  Teuchos::ESide result;
  int r = GetRandom(1, 2);
  if(r == 1)
    {
      result = Teuchos::LEFT_SIDE;
    }
  else
    {
      result = Teuchos::RIGHT_SIDE;
    }
  return result;
}

Teuchos::EUplo RandomUPLO()
{
  Teuchos::EUplo result;
  int r = GetRandom(1, 2);
  if(r == 1)
    {
      result = Teuchos::UPPER_TRI;
    }
  else
    {
      result = Teuchos::LOWER_TRI;
    }
  return result;
}

Teuchos::ETransp RandomTRANS()
{
  Teuchos::ETransp result;
  int r = GetRandom(1, 4);
  if(r == 1 || r == 2)
    {
      result = Teuchos::NO_TRANS;
    }
  else if(r == 3)
    {
      result = Teuchos::TRANS;
    }
  else
    {
      result = Teuchos::CONJ_TRANS;
    }
  return result;
}

Teuchos::EDiag RandomDIAG()
{
  Teuchos::EDiag result;
  int r = GetRandom(1, 2);
  if(r == 1)
    {
      result = Teuchos::NON_UNIT_DIAG;
    }
  else
    {
      result = Teuchos::UNIT_DIAG;
    }
  return result;
}

