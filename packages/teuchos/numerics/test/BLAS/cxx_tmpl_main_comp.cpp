// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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

// OType1 and OType2 define the ordinal datatypes for which BLAS output will be compared.
// The difference in OType should enable the comparison of the templated routines with the "officially" supported BLAS.

// Define the scalar type
#ifdef HAVE_TEUCHOS_COMPLEX
#define SType     std::complex<double>
#else
#define SType     double
#endif

// Define the ordinal type
#define OType1	   long int 
#define OType2	   int 

// MVMIN/MAX define the minimum and maximum dimensions of generated matrices and vectors, respectively.
// These are well within the range of OType1 and OType2
#define MVMIN      2
#define MVMAX      20
// SCALARMAX defines the maximum positive value (with a little leeway) generated for matrix and std::vector elements and scalars:
// random numbers in [-SCALARMAX, SCALARMAX] will be generated.
// Set SCALARMAX to a floating-point value (e.g. 10.0) to enable floating-point random number generation, such that
// random numbers in (-SCALARMAX - 1, SCALARMAX + 1) will be generated.
#ifdef HAVE_TEUCHOS_COMPLEX
#define SCALARMAX  SType(10,0)
#else
#define SCALARMAX  SType(10)
#endif
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

template<typename T>
std::complex<T> GetRandom( std::complex<T>, std::complex<T> );

template<typename TYPE, typename OTYPE>
void PrintVector(TYPE* Vector, OTYPE Size, std::string Name, bool Matlab = 0);

template<typename TYPE, typename OTYPE>
void PrintMatrix(TYPE* Matrix, OTYPE Rows, OTYPE Columns, OTYPE LDM, std::string Name, bool Matlab = 0);

template<typename TYPE>
bool CompareScalars(TYPE Scalar1, TYPE Scalar2, typename ScalarTraits<TYPE>::magnitudeType Tolerance ); 

template<typename TYPE, typename OTYPE1, typename OTYPE2>
bool CompareVectors(TYPE* Vector1, OTYPE1 Size1, TYPE* Vector2, OTYPE2 Size2, typename ScalarTraits<TYPE>::magnitudeType Tolerance ); 

template<typename TYPE, typename OTYPE1, typename OTYPE2>
bool CompareMatrices(TYPE* Matrix1, OTYPE1 Rows1, OTYPE1 Columns1, OTYPE1 LDM1, 
                     TYPE* Matrix2, OTYPE2 Rows2, OTYPE2 Columns2, OTYPE2 LDM2,
                     typename ScalarTraits<TYPE>::magnitudeType Tolerance ); 

// Use this to convert the larger ordinal type to the smaller one (nothing inherently makes sure of this).
template<typename OTYPE1, typename OTYPE2>
OTYPE2 ConvertType(OTYPE1 T1, OTYPE2 T2)
{
  return static_cast<OTYPE2>(T1);
}

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
  int i;
  OType1 j1, k1;
  OType2 j2, k2;
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
  typedef ScalarTraits<SType>::magnitudeType MType;
  BLAS<OType1, SType> OType1BLAS;
  BLAS<OType2, SType> OType2BLAS;
  SType STypezero = ScalarTraits<SType>::zero();
  SType STypeone = ScalarTraits<SType>::one();
  SType OType1alpha, OType1beta;
  SType OType2alpha, OType2beta;
  SType *OType1A, *OType1B, *OType1C, *OType1x, *OType1y;
  SType *OType2A, *OType2B, *OType2C, *OType2x, *OType2y; 
  SType OType1ASUMresult, OType1DOTresult, OType1NRM2result, OType1SINresult;
  SType OType2ASUMresult, OType2DOTresult, OType2NRM2result, OType2SINresult;
  MType OType1COSresult, OType2COSresult;
  OType1 incx1, incy1;
  OType2 incx2, incy2;
  OType1 OType1IAMAXresult;
  OType2 OType2IAMAXresult;
  OType1 TotalTestCount = 1, GoodTestSubcount, GoodTestCount = 0, M1, N1, P1, K1, LDA1, LDB1, LDC1, Mx1, My1;
  OType2 M2, N2, P2, K2, LDA2, LDB2, LDC2, Mx2, My2;
  Teuchos::EUplo UPLO;
  Teuchos::ESide SIDE;
  Teuchos::ETransp TRANS, TRANSA, TRANSB;
  Teuchos::EDiag DIAG;
  MType TOL = 1e-8*ScalarTraits<MType>::one();
  
  std::srand(time(NULL));

  //--------------------------------------------------------------------------------
  // BEGIN LEVEL 1 BLAS TESTS
  //--------------------------------------------------------------------------------
  // Begin ROTG Tests
  //--------------------------------------------------------------------------------
  GoodTestSubcount = 0;
  for(i = 0; i < ROTGTESTS; i++)
    {
      OType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      OType2alpha = OType1alpha;
      OType1beta = GetRandom(-SCALARMAX, SCALARMAX);
      OType2beta = OType1beta;
      OType1COSresult = ScalarTraits<MType>::zero();
      OType2COSresult = OType1COSresult;
      OType1SINresult = ScalarTraits<SType>::zero();
      OType2SINresult = OType1SINresult;
      
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- ROTG -- " << std::endl;
	  std::cout << "OType1alpha = "  << OType1alpha << std::endl;
	  std::cout << "OType2alpha = " << OType2alpha << std::endl;
	  std::cout << "OType1beta = "  << OType1beta << std::endl;
	  std::cout << "OType2beta = " << OType2beta << std::endl;
	}
      TotalTestCount++;
      OType1BLAS.ROTG(&OType1alpha, &OType1beta, &OType1COSresult, &OType1SINresult);
      OType2BLAS.ROTG(&OType2alpha, &OType2beta, &OType2COSresult, &OType2SINresult);
      if(debug)
	{
	  std::cout << "OType1 ROTG COS result: " << OType1COSresult << std::endl;
	  std::cout << "OType2 ROTG COS result: " << OType2COSresult << std::endl;
	  std::cout << "OType1 ROTG SIN result: " << OType1SINresult << std::endl;
	  std::cout << "OType2 ROTG SIN result: " << OType2SINresult << std::endl;
	}
      if ( !CompareScalars(OType1COSresult, OType2COSresult, TOL) || !CompareScalars(OType1SINresult, OType2SINresult, TOL) )
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += ( CompareScalars(OType1COSresult, OType2COSresult, TOL) && 
			    CompareScalars(OType1SINresult, OType2SINresult, TOL) );
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
  GoodTestSubcount = 0;
  for(i = 0; i < ROTTESTS; i++)
  {
    incx1 = GetRandom(-5,5);
    incy1 = GetRandom(-5,5);
    if (incx1 == 0) incx1 = 1;
    if (incy1 == 0) incy1 = 1;
    incx2 = ConvertType( incx1, incx2 ); 
    incy2 = ConvertType( incy1, incy2 ); 
    M1 = GetRandom(MVMIN, MVMIN+8);
    M2 = ConvertType( M1, M2 );
    Mx1 = M1*std::abs(incx1);
    My1 = M1*std::abs(incy1);
    if (Mx1 == 0) { Mx1 = 1; }
    if (My1 == 0) { My1 = 1; }
    Mx2 = ConvertType( Mx1, Mx2 ); 
    My2 = ConvertType( My1, My2 );
    OType1x = new SType[Mx1];
    OType1y = new SType[My1];
    OType2x = new SType[Mx2];
    OType2y = new SType[My2];
    for(j1 = 0, j2 = 0; j1 < Mx1; j1++, j2++)
    {
      OType1x[j1] = GetRandom(-SCALARMAX, SCALARMAX);
      OType2x[j2] = OType1x[j1];
    }
    for(j1 = 0, j2 = 0; j1 < My1; j1++, j2++)
    {
      OType1y[j1] = GetRandom(-SCALARMAX, SCALARMAX);
      OType2y[j2] = OType1y[j1];
    }
    MType c1 = cos(ScalarTraits<SType>::magnitude(GetRandom(-SCALARMAX,SCALARMAX)));
    MType c2 = c1;
    SType s1 = sin(ScalarTraits<SType>::magnitude(GetRandom(-SCALARMAX,SCALARMAX)));
    SType s2 = s1;
    if(debug)
    {
      std::cout << "Test #" << TotalTestCount << " -- ROT -- " << std::endl;
      std::cout << "c1 = "  << c1 << ", s1 = " << s1 << std::endl;
      std::cout << "c2 = " << c2 << ", s2 = " << s2 << std::endl;
      std::cout << "incx1 = " << incx1 << ", incy1 = " << incy1 << std::endl;
      std::cout << "incx2 = " << incx2 << ", incy2 = " << incy2 << std::endl;
      PrintVector(OType1x, Mx1, "OType1x", matlab);
      PrintVector(OType1y, My1, "OType1y_before_operation", matlab);
      PrintVector(OType2x, Mx2, "OType2x", matlab);
      PrintVector(OType2y, My2, "OType2y_before_operation",  matlab);
    }
    TotalTestCount++;
    OType1BLAS.ROT(M1, OType1x, incx1, OType1y, incy1, &c1, &s1);
    OType2BLAS.ROT(M2, OType2x, incx2, OType2y, incy2, &c2, &s2);
    if(debug)
    {
      PrintVector(OType1y, My1, "OType1y_after_operation", matlab);
      PrintVector(OType2y, My2, "OType2y_after_operation", matlab);
    }
    if ( !CompareVectors(OType1x, Mx1, OType2x, Mx2, TOL) || !CompareVectors(OType1y, My1, OType2y, My2, TOL) )
	std::cout << "FAILED TEST!!!!!!" << std::endl;
    GoodTestSubcount += ( CompareVectors(OType1x, Mx1, OType2x, Mx2, TOL) &&
                          CompareVectors(OType1y, My1, OType2y, My2, TOL) );
    delete [] OType1x;
    delete [] OType1y;
    delete [] OType2x;
    delete [] OType2y;
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
      incx1 = GetRandom(1, MVMAX);
      incx2 = ConvertType( incx1, incx2 );
      M1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      OType1x = new SType[M1*incx1];
      OType2x = new SType[M2*incx2];
      for(j1 = 0, j2 = 0; j2 < M2*incx2; j1++, j2++)
	{
	  OType1x[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2x[j2] = OType1x[j1];
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- ASUM -- " << std::endl;
          std::cout << "incx1 = " << incx1 << "\t" << "incx2 = " << incx2
                    << "\t" << "M1 = " << M1 << "\t" << "M2 = " << M2 << std::endl;
	  PrintVector(OType1x, M1*incx1, "OType1x", matlab);
	  PrintVector(OType2x, M2*incx2, "OType2x", matlab);
	}
      TotalTestCount++;
      OType1ASUMresult = OType1BLAS.ASUM(M1, OType1x, incx1);
      OType2ASUMresult = OType2BLAS.ASUM(M2, OType2x, incx2);
      if(debug)
	{
	  std::cout << "OType1 ASUM result: " << OType1ASUMresult << std::endl;
	  std::cout << "OType2 ASUM result: " << OType2ASUMresult << std::endl;
	}
      if (CompareScalars(OType1ASUMresult, OType2ASUMresult, TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareScalars(OType1ASUMresult, OType2ASUMresult, TOL);

      delete [] OType1x;
      delete [] OType2x;
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
      incx1 = GetRandom(1, MVMAX);
      incy1 = GetRandom(1, MVMAX);
      incx2 = ConvertType( incx1, incx2 );
      incy2 = ConvertType( incy1, incy2 );
      M1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      Mx1 = M1*std::abs(incx1);
      My1 = M1*std::abs(incy1);
      if (Mx1 == 0) { Mx1 = 1; }
      if (My1 == 0) { My1 = 1; }
      Mx2 = ConvertType( Mx1, Mx2 );
      My2 = ConvertType( My1, My2 );
      OType1x = new SType[Mx1];
      OType1y = new SType[My1];
      OType2x = new SType[Mx2];
      OType2y = new SType[My2]; 
      for(j1 = 0, j2 = 0; j1 < Mx1; j1++, j2++)
	{
	  OType1x[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2x[j2] = OType1x[j1];
	}
      for(j1 = 0, j2 = 0; j1 < My1; j1++, j2++)
	{
	  OType1y[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2y[j2] = OType1y[j1];
	}
      OType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      OType2alpha = OType1alpha;
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- AXPY -- " << std::endl;
	  std::cout << "OType1alpha = "  << OType1alpha << std::endl;
	  std::cout << "OType2alpha = " << OType2alpha << std::endl;
	  PrintVector(OType1x, Mx1, "OType1x", matlab);
	  PrintVector(OType1y, My1, "OType1y_before_operation", matlab);
	  PrintVector(OType2x, Mx2, "OType2x", matlab);
	  PrintVector(OType2y, My2, "OType2y_before_operation",  matlab);
	}
      TotalTestCount++;
      OType1BLAS.AXPY(M1, OType1alpha, OType1x, incx1, OType1y, incy1);
      OType2BLAS.AXPY(M2, OType2alpha, OType2x, incx2, OType2y, incy2);
      if(debug)
	{
	  PrintVector(OType1y, My1, "OType1y_after_operation", matlab);
	  PrintVector(OType2y, My2, "OType2y_after_operation", matlab);
	}
      if (CompareVectors(OType1y, My1, OType2y, My2, TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareVectors(OType1y, My1, OType2y, My2, TOL);

      delete [] OType1x;
      delete [] OType1y;
      delete [] OType2x;
      delete [] OType2y;
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
      incx1 = GetRandom(1, MVMAX);
      incy1 = GetRandom(1, MVMAX);
      incx2 = ConvertType( incx1, incx2 );
      incy2 = ConvertType( incy1, incy2 );
      M1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      Mx1 = M1*std::abs(incx1);
      My1 = M1*std::abs(incy1);
      if (Mx1 == 0) { Mx1 = 1; }
      if (My1 == 0) { My1 = 1; }
      Mx2 = ConvertType( Mx1, Mx2 );
      My2 = ConvertType( My1, My2 );
      OType1x = new SType[Mx1];
      OType1y = new SType[My1];
      OType2x = new SType[Mx2];
      OType2y = new SType[My2]; 
      for(j1 = 0, j2 = 0; j1 < Mx1; j1++, j2++)
	{
	  OType1x[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2x[j2] = OType1x[j1];
	}
      for(j1 = 0, j2 = 0; j1 < My1; j1++, j2++)
	{
	  OType1y[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2y[j2] = OType1y[j1];
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- COPY -- " << std::endl;
	  PrintVector(OType1x, Mx1, "OType1x", matlab);
	  PrintVector(OType1y, My1, "OType1y_before_operation", matlab);
	  PrintVector(OType2x, Mx2, "OType2x", matlab);
	  PrintVector(OType2y, My2, "OType2y_before_operation", matlab);
	}
      TotalTestCount++;
      OType1BLAS.COPY(M1, OType1x, incx1, OType1y, incy1);
      OType2BLAS.COPY(M2, OType2x, incx2, OType2y, incy2);
      if(debug)
	{
	  PrintVector(OType1y, My1, "OType1y_after_operation", matlab);
	  PrintVector(OType2y, My2, "OType2y_after_operation", matlab);
	}
      if (CompareVectors(OType1y, My1, OType2y, My2, TOL) == 0 )
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareVectors(OType1y, My1, OType2y, My2, TOL);

      delete [] OType1x;
      delete [] OType1y;
      delete [] OType2x;
      delete [] OType2y;
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
      incx1 = GetRandom(1, MVMAX);
      incy1 = GetRandom(1, MVMAX);
      incx2 = ConvertType( incx1, incx2 );
      incy2 = ConvertType( incy1, incy2 );
      M1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      Mx1 = M1*std::abs(incx1);
      My1 = M1*std::abs(incy1);
      if (Mx1 == 0) { Mx1 = 1; }
      if (My1 == 0) { My1 = 1; }
      Mx2 = ConvertType( Mx1, Mx2 );
      My2 = ConvertType( My1, My2 );
      OType1x = new SType[Mx1];
      OType1y = new SType[My1];
      OType2x = new SType[Mx2];
      OType2y = new SType[My2]; 
      for(j1 = 0, j2 = 0; j1 < Mx1; j1++, j2++)
	{
	  OType1x[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2x[j2] = OType1x[j1];
	}
      for(j1 = 0, j2 = 0; j1 < My1; j1++, j2++)
	{
	  OType1y[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2y[j2] = OType1y[j1];
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- DOT -- " << std::endl;
	  PrintVector(OType1x, Mx1, "OType1x", matlab);
	  PrintVector(OType1y, My1, "OType1y", matlab);
	  PrintVector(OType2x, Mx2, "OType2x", matlab);
	  PrintVector(OType2y, My2, "OType2y", matlab);
	}
      TotalTestCount++;
      OType1DOTresult = OType1BLAS.DOT(M1, OType1x, incx1, OType1y, incy1);
      OType2DOTresult = OType2BLAS.DOT(M2, OType2x, incx2, OType2y, incy2);
      if(debug)
	{
	  std::cout << "OType1 DOT result: " << OType1DOTresult << std::endl;
	  std::cout << "OType2 DOT result: " << OType2DOTresult << std::endl;
	}
      if (CompareScalars(OType1DOTresult, OType2DOTresult, TOL) == 0) {
	std::cout << "DOT test " << i+1 << " of " << DOTTESTS << " FAILED!  "
		  << "SType = " << Teuchos::TypeNameTraits<SType>::name () << ".  "
		  << "The two results are " << OType1DOTresult << " and " 
		  << OType2DOTresult << ".  incx1 = " << incx1 << ", incy1 = " 
		  << incy1 << ", incx2 = " << incx2 << ", and incy2 = " 
		  << incy2 << std::endl;
      }

      GoodTestSubcount += CompareScalars(OType1DOTresult, OType2DOTresult, TOL);

      delete [] OType1x;
      delete [] OType1y;
      delete [] OType2x;
      delete [] OType2y;
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
      incx1 = GetRandom(1, MVMAX);
      incx2 = ConvertType( incx1, incx2 );
      M1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      OType1x = new SType[M1*incx1];
      OType2x = new SType[M2*incx2];
      for(j1 = 0, j2 = 0; j1 < M1*incx1; j1++, j2++)
	{
	  OType1x[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2x[j2] = OType1x[j1];
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- NRM2 -- " << std::endl;
	  PrintVector(OType1x, M1*incx1, "OType1x", matlab);
	  PrintVector(OType2x, M2*incx2, "OType2x", matlab);
	}
      TotalTestCount++;
      OType1NRM2result = OType1BLAS.NRM2(M1, OType1x, incx1);
      OType2NRM2result = OType2BLAS.NRM2(M2, OType2x, incx2);
      if(debug)
	{
	  std::cout << "OType1 NRM2 result: " << OType1NRM2result << std::endl;
	  std::cout << "OType2 NRM2 result: " << OType2NRM2result << std::endl;
	}
      if (CompareScalars(OType1NRM2result, OType2NRM2result, TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareScalars(OType1NRM2result, OType2NRM2result, TOL);

      delete [] OType1x;
      delete [] OType2x;
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
      incx1 = GetRandom(1, MVMAX);
      incx2 = ConvertType( incx1, incx2 );
      M1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      OType1x = new SType[M1*incx1];
      OType2x = new SType[M2*incx2];
      for(j1 = 0, j2 = 0; j1 < M1*incx1; j1++, j2++)
	{
	  OType1x[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2x[j2] = OType1x[j1];
	}
      OType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      OType2alpha = OType1alpha;
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- SCAL -- " << std::endl;
	  std::cout << "OType1alpha = " << OType1alpha << std::endl;
	  std::cout << "OType2alpha = " << OType2alpha << std::endl;
	  PrintVector(OType1x, M1*incx1, "OType1x_before_operation", matlab);
	  PrintVector(OType2x, M2*incx2, "OType2x_before_operation", matlab);
	}
      TotalTestCount++;
      OType1BLAS.SCAL(M1, OType1alpha, OType1x, incx1);
      OType2BLAS.SCAL(M2, OType2alpha, OType2x, incx2);
      if(debug)
	{
	  PrintVector(OType1x, M1*incx1, "OType1x_after_operation", matlab);
	  PrintVector(OType2x, M2*incx2, "OType2x_after_operation", matlab);
	}
      if (CompareVectors(OType1x, M1*incx1, OType2x, M2*incx2, TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareVectors(OType1x, M1*incx1, OType2x, M2*incx2, TOL);

      delete [] OType1x;
      delete [] OType2x;
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
      incx1 = GetRandom(1, MVMAX);
      incx2 = ConvertType( incx1, incx2 );
      M1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      OType1x = new SType[M1*incx1];
      OType2x = new SType[M2*incx2];
      for(j1 = 0, j2 = 0; j1 < M1*incx1; j1++, j2++)
	{
	  OType1x[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2x[j2] = OType1x[j1];
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- IAMAX -- " << std::endl;
	  PrintVector(OType1x, M1*incx1, "OType1x", matlab);
	  PrintVector(OType2x, M2*incx2, "OType2x", matlab);
	}
      TotalTestCount++;
      OType1IAMAXresult = OType1BLAS.IAMAX(M1, OType1x, incx1);
      OType2IAMAXresult = OType2BLAS.IAMAX(M2, OType2x, incx2);
      if(debug)
	{
	  std::cout << "OType1 IAMAX result: " << OType1IAMAXresult << std::endl;
	  std::cout << "OType2 IAMAX result: " << OType2IAMAXresult << std::endl;
	}
      if (OType1IAMAXresult != OType2IAMAXresult)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += (OType1IAMAXresult == OType2IAMAXresult);

      delete [] OType1x;
      delete [] OType2x;
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
      incx1 = GetRandom(1, MVMAX);
      while (incx1 == 0) {
      	  incx1 = GetRandom(1, MVMAX);
      }
      incy1 = GetRandom(1, MVMAX);
      while (incy1 == 0) {
      	  incy1 = GetRandom(1, MVMAX);
      }
      incx2 = ConvertType( incx1, incx2 );   
      incy2 = ConvertType( incy1, incy2 );   
      M1 = GetRandom(MVMIN, MVMAX);
      N1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      N2 = ConvertType( N1, N2 );

      TRANS = RandomTRANS();
      OType1 M2_1 = 0, N2_1 = 0;
      if (Teuchos::ETranspChar[TRANS] == 'N') {	
      	M2_1 = M1*std::abs(incy1);
      	N2_1 = N1*std::abs(incx1);   
      } else {
	M2_1 = N1*std::abs(incy1);
	N2_1 = M1*std::abs(incx1);
      }
      OType2 M2_2 = ConvertType( M2_1, M2_2 );
      OType2 N2_2 = ConvertType( N2_1, N2_2 );

      LDA1 = GetRandom(MVMIN, MVMAX);
      while (LDA1 < M1) {
          LDA1 = GetRandom(MVMIN, MVMAX);
      }   
      LDA2 = ConvertType( LDA1, LDA2 );

      OType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      OType1beta = GetRandom(-SCALARMAX, SCALARMAX);
      OType2alpha = OType1alpha;
      OType2beta = OType1beta;

      OType1A = new SType[LDA1 * N1];
      OType1x = new SType[N2_1];
      OType1y = new SType[M2_1];
      OType2A = new SType[LDA2 * N2];
      OType2x = new SType[N2_2];
      OType2y = new SType[M2_2]; 

      for(j1 = 0, j2 = 0; j1 < LDA1 * N1; j1++, j2++)
	{
	  OType1A[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2A[j2] = OType1A[j1];
	}
      for(j1 = 0, j2 = 0; j1 < N2_1; j1++, j2++)
	{
	  OType1x[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2x[j2] = OType1x[j1];
	}
      for(j1 = 0, j2 = 0; j1 < M2_1; j1++, j2++)
	{
	  OType1y[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2y[j2] = OType1y[j1];
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- GEMV -- " << std::endl;
	  std::cout << "TRANS = " << Teuchos::ETranspChar[TRANS] << std::endl;
	  std::cout << "OType1alpha = " << OType1alpha << std::endl;
	  std::cout << "OType2alpha = " << OType2alpha << std::endl;
	  std::cout << "OType1beta = " << OType1beta << std::endl;
	  std::cout << "OType2beta = " << OType2beta << std::endl;
	  PrintMatrix(OType1A, M1, N1, LDA1, "OType1A", matlab);
	  PrintVector(OType1x, N2_1, "OType1x", matlab);
	  PrintVector(OType1y, M2_1, "OType1y_before_operation", matlab);
	  PrintMatrix(OType2A, M2, N2, LDA2, "OType2A", matlab);
	  PrintVector(OType2x, N2_2, "OType2x", matlab);
	  PrintVector(OType2y, M2_2, "OType2y_before_operation", matlab);
	}
      TotalTestCount++;
      OType1BLAS.GEMV(TRANS, M1, N1, OType1alpha, OType1A, LDA1, OType1x, incx1, OType1beta, OType1y, incy1);
      OType2BLAS.GEMV(TRANS, M2, N2, OType2alpha, OType2A, LDA2, OType2x, incx2, OType2beta, OType2y, incy2);
      if(debug)
	{
	  PrintVector(OType1y, M2_1, "OType1y_after_operation", matlab);
	  PrintVector(OType2y, M2_2, "OType2y_after_operation", matlab);
	}
      if (CompareVectors(OType1y, M2_1, OType2y, M2_2, TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareVectors(OType1y, M2_1, OType2y, M2_2, TOL);

      delete [] OType1A;
      delete [] OType1x;
      delete [] OType1y;
      delete [] OType2A;
      delete [] OType2x;
      delete [] OType2y;
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

      N1 = GetRandom(MVMIN, MVMAX);
      N2 = ConvertType( N1, N2 );
      incx1 = GetRandom(1, MVMAX);
      while (incx1 == 0) {
      	  incx1 = GetRandom(1, MVMAX);
      }
      incx2 = ConvertType( incx1, incx2 );
      OType1x = new SType[N1*std::abs(incx1)];
      OType2x = new SType[N2*std::abs(incx2)];

      for(j1 = 0, j2 = 0; j1 < N1*std::abs(incx1); j1++, j2++)
	{
	  OType1x[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2x[j2] = OType1x[j1];
	}

      LDA1 = GetRandom(MVMIN, MVMAX);
      while (LDA1 < N1) {
	LDA1 = GetRandom(MVMIN, MVMAX);
      }
      LDA2 = ConvertType( LDA1, LDA2 );
      OType1A = new SType[LDA1 * N1];
      OType2A = new SType[LDA2 * N2];

      for(j1 = 0, j2 = 0; j1 < N1; j1++, j2++)
	{	     
	  if(Teuchos::EUploChar[UPLO] == 'U') {
	    // The operator is upper triangular, make sure that the entries are
	    // only in the upper triangular part of A and the diagonal is non-zero.
	    for(k1 = 0, k2 = 0; k1 < N1; k1++, k2++) 
	    {
	      if(k1 < j1) {
		OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
	      } else {
		OType1A[j1*LDA1+k1] = STypezero;
	      }
	      OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
	      if(k1 == j1) {
		if (Teuchos::EDiagChar[DIAG] == 'N') {
		  OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		  while (OType1A[j1*LDA1+k1] == STypezero) {
		    OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		  }
		  OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		} else {
		  OType1A[j1*LDA1+k1] = STypeone;
		  OType2A[j2*LDA2+k2] = STypeone;
		}
	      }			
	    }
	  } else {
	    // The operator is lower triangular, make sure that the entries are
	    // only in the lower triangular part of A and the diagonal is non-zero.
	    for(k1 = 0, k2 = 0; k1 < N1; k1++, k2++) 
	      {
		if(k1 > j1) {
		  OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		} else {
		  OType1A[j1*LDA1+k1] = STypezero;
		}
		OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		if(k1 == j1) {
		  if (Teuchos::EDiagChar[DIAG] == 'N') {
		    OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		    while (OType1A[j1*LDA1+k1] == STypezero) {
		      OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		    }
		    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		  } else {
		    OType1A[j1*LDA1+k1] = STypeone;
		    OType2A[j2*LDA2+k2] = STypeone;
		  }
		}			
	      } // end for(k=0 ...		
	  } // end if(UPLO == 'U') ...
	} // end for(j=0 ...      for(j = 0; j < N*N; j++)
      
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- TRMV -- " << std::endl;
	  std::cout << "UPLO = " << Teuchos::EUploChar[UPLO] << "\t" 
	       << "TRANSA = " << Teuchos::ETranspChar[TRANSA] << "\t" 
	       << "DIAG = " << Teuchos::EDiagChar[DIAG] << std::endl;
	  PrintMatrix(OType1A, N1, N1, LDA1,"OType1A", matlab);
	  PrintVector(OType1x, N1*incx1, "OType1x_before_operation", matlab);
	  PrintMatrix(OType2A, N2, N2, LDA2, "OType2A", matlab);
	  PrintVector(OType2x, N2*incx2, "OType2x_before_operation", matlab);
	}
      TotalTestCount++;
      OType1BLAS.TRMV(UPLO, TRANSA, DIAG, N1, OType1A, LDA1, OType1x, incx1);
      OType2BLAS.TRMV(UPLO, TRANSA, DIAG, N2, OType2A, LDA2, OType2x, incx2);
      if(debug)
	{
	  PrintVector(OType1x, N1*incx1, "OType1x_after_operation", matlab);
	  PrintVector(OType2x, N2*incx2, "OType2x_after_operation", matlab);
	}
      if (CompareVectors(OType1x, std::abs(N1*incx1), OType2x, std::abs(N2*incx2), TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareVectors(OType1x, std::abs(N1*incx1), OType2x, std::abs(N2*incx2), TOL);

      delete [] OType1A;
      delete [] OType1x;
      delete [] OType2A;
      delete [] OType2x;
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
      incx1 = GetRandom(1, MVMAX);
      while (incx1 == 0) {
      	  incx1 = GetRandom(1, MVMAX);
      }   
      incy1 = GetRandom(1, MVMAX);
      while (incy1 == 0) {
      	  incy1 = GetRandom(1, MVMAX);
      }   
      incx2 = ConvertType( incx1, incx2 );
      incy2 = ConvertType( incy1, incy2 );
      M1 = GetRandom(MVMIN, MVMAX);
      N1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      N2 = ConvertType( N1, N2 );

      LDA1 = GetRandom(MVMIN, MVMAX);
      while (LDA1 < M1) {
          LDA1 = GetRandom(MVMIN, MVMAX);
      }   
      LDA2 = ConvertType( LDA1, LDA2 );

      OType1A = new SType[LDA1 * N1];
      OType1x = new SType[M1*std::abs(incx1)];
      OType1y = new SType[N1*std::abs(incy1)];
      OType2A = new SType[LDA2 * N2];
      OType2x = new SType[M2*std::abs(incx2)];
      OType2y = new SType[N2*std::abs(incy2)];
      OType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      OType2alpha = OType1alpha;
      for(j1 = 0, j2 = 0; j1 < LDA1 * N1; j1++, j2++)
	{
	  OType1A[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2A[j2] = OType1A[j1];
	}
      for(j1 = 0, j2 = 0; j1 < std::abs(M1*incx1); j1++, j2++)
	{
	  OType1x[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2x[j2] = OType1x[j1];
	}
      for(j1 = 0, j2 = 0; j1 < std::abs(N1*incy1); j1++, j2++)
	{
	  OType1y[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2y[j2] = OType1y[j1];
	}
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- GER -- " << std::endl;
	  std::cout << "OType1alpha = " << OType1alpha << std::endl;
	  std::cout << "OType2alpha = " << OType2alpha << std::endl;
	  PrintMatrix(OType1A, M1, N1, LDA1,"OType1A_before_operation", matlab);
	  PrintVector(OType1x, std::abs(M1*incx1), "OType1x", matlab);
	  PrintVector(OType1y, std::abs(N1*incy1), "OType1y", matlab);
	  PrintMatrix(OType2A, M2, N2, LDA2,"OType2A_before_operation", matlab);
	  PrintVector(OType2x, std::abs(M2*incx2), "OType2x", matlab);
	  PrintVector(OType2y, std::abs(N2*incy2), "OType2y", matlab);
	}
      TotalTestCount++;
      OType1BLAS.GER(M1, N1, OType1alpha, OType1x, incx1, OType1y, incy1, OType1A, LDA1);
      OType2BLAS.GER(M2, N2, OType2alpha, OType2x, incx2, OType2y, incy2, OType2A, LDA2);
      if(debug)
	{
	  PrintMatrix(OType1A, M1, N1, LDA1, "OType1A_after_operation", matlab);
	  PrintMatrix(OType2A, M2, N2, LDA2, "OType2A_after_operation", matlab);
	}
      if (CompareMatrices(OType1A, M1, N1, LDA1, OType2A, M2, N2, LDA2, TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareMatrices(OType1A, M1, N1, LDA1, OType2A, M2, N2, LDA2, TOL);

      delete [] OType1A;
      delete [] OType1x;
      delete [] OType1y;
      delete [] OType2A;
      delete [] OType2x;
      delete [] OType2y;
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
      M1 = GetRandom(MVMIN, MVMAX);
      N1 = GetRandom(MVMIN, MVMAX);
      P1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      N2 = ConvertType( N1, N2 );
      P2 = ConvertType( P1, P2 );

      if(debug)	{
	  std::cout << "Test #" << TotalTestCount << " -- GEMM -- " << std::endl;
	  std::cout << "TRANSA = " << Teuchos::ETranspChar[TRANSA] << "\t" 
	       << "TRANSB = " << Teuchos::ETranspChar[TRANSB] << std::endl; 
      }
      LDA1 = GetRandom(MVMIN, MVMAX);
      if (Teuchos::ETranspChar[TRANSA] == 'N') {
	while (LDA1 < M1) {  LDA1 = GetRandom(MVMIN, MVMAX); }
        LDA2 = ConvertType( LDA1, LDA2 );
	OType1A = new SType[LDA1 * P1];
	OType2A = new SType[LDA2 * P2];
	for(j1 = 0, j2 = 0; j1 < LDA1 * P1; j1++, j2++)
	{
	  OType1A[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2A[j2] = OType1A[j1];
	}
	if (debug) {
	  PrintMatrix(OType1A, M1, P1, LDA1, "OType1A", matlab);
	  PrintMatrix(OType2A, M2, P2, LDA2, "OType2A", matlab);
	}
      } else {
	while (LDA1 < P1) {  LDA1 = GetRandom(MVMIN, MVMAX); }
        LDA2 = ConvertType( LDA1, LDA2 );
	OType1A = new SType[LDA1 * M1];
	OType2A = new SType[LDA1 * M1];
	for(j1 = 0, j2 = 0; j1 < LDA1 * M1; j1++, j2++)
	{
	  OType1A[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2A[j2] = OType1A[j1];
	}
	if (debug) {
  	  PrintMatrix(OType1A, P1, M1, LDA1, "OType1A", matlab);
	  PrintMatrix(OType2A, P2, M2, LDA2, "OType2A", matlab);
	}
      }

      LDB1 = GetRandom(MVMIN, MVMAX);
      if (Teuchos::ETranspChar[TRANSB] == 'N') {
	while (LDB1 < P1) {  LDB1 = GetRandom(MVMIN, MVMAX); }
        LDB2 = ConvertType( LDB1, LDB2 );
	OType1B = new SType[LDB1 * N1];
	OType2B = new SType[LDB2 * N2];
	for(j1 = 0, j2 = 0; j1 < LDB1 * N1; j1++, j2++)
	{
	  OType1B[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2B[j2] = OType1B[j1];
	}
	if (debug) {
	  PrintMatrix(OType1B, P1, N1, LDB1,"OType1B", matlab);
	  PrintMatrix(OType2B, P2, N2, LDB2,"OType2B", matlab);
	}
      } else { 
	while (LDB1 < N1) {  LDB1 = GetRandom(MVMIN, MVMAX); }
        LDB2 = ConvertType( LDB1, LDB2 );
	OType1B = new SType[LDB1 * P1];
	OType2B = new SType[LDB2 * P2];
	for(j1 = 0, j2 = 0; j1 < LDB1 * P1; j1++, j2++)
	{
	  OType1B[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2B[j2] = OType1B[j1];
	}	
	if (debug) {
	  PrintMatrix(OType1B, N1, P1, LDB1,"OType1B", matlab);
	  PrintMatrix(OType2B, N2, P2, LDB2,"OType2B", matlab);
	}
      }

      LDC1 = GetRandom(MVMIN, MVMAX);
      while (LDC1 < M1) {  LDC1 = GetRandom(MVMIN, MVMAX); }
      LDC2 = ConvertType( LDC1, LDC2 );
      OType1C = new SType[LDC1 * N1];
      OType2C = new SType[LDC2 * N2];
      for(j1 = 0, j2 = 0; j1 < LDC1 * N1; j1++, j2++) {
	  OType1C[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2C[j2] = OType1C[j1];
      }
      if(debug)
	{
	  PrintMatrix(OType1C, M1, N1, LDC1, "OType1C_before_operation", matlab);
	  PrintMatrix(OType2C, M2, N2, LDC2, "OType2C_before_operation", matlab);
	}
	
      OType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      OType1beta = GetRandom(-SCALARMAX, SCALARMAX);
      OType2alpha = OType1alpha;
      OType2beta = OType1beta;

      TotalTestCount++;
      OType1BLAS.GEMM(TRANSA, TRANSB, M1, N1, P1, OType1alpha, OType1A, LDA1, OType1B, LDB1, OType1beta, OType1C, LDC1);
      OType2BLAS.GEMM(TRANSA, TRANSB, M2, N2, P2, OType2alpha, OType2A, LDA2, OType2B, LDB2, OType2beta, OType2C, LDC2);
      if(debug)
	{
	  std::cout << "M1="<<M1 << "\t" << "N1="<<N1 << "\t" << "P1 = " << P1 
                    << "\t" << "LDA1="<<LDA1 << "\t" << "LDB1="<<LDB1 << "\t" << "LDC1=" << LDC1 << std::endl;
	  std::cout << "M2="<<M2 << "\t" << "N2="<<N2 << "\t" << "P2 = " << P2 
                    << "\t" << "LDA2="<<LDA2 << "\t" << "LDB2="<<LDB2 << "\t" << "LDC2=" << LDC2 << std::endl;
	  std::cout << "OType1alpha = " << OType1alpha << std::endl;
	  std::cout << "OType2alpha = " << OType2alpha << std::endl;
	  std::cout << "OType1beta = " << OType1beta << std::endl;
	  std::cout << "OType2beta = " << OType2beta << std::endl;
	  PrintMatrix(OType1C, M1, N1, LDC1, "OType1C_after_operation", matlab);
	  PrintMatrix(OType2C, M2, N2, LDC2, "OType2C_after_operation", matlab);
	}
      if (CompareMatrices(OType1C, M1, N1, LDC1, OType2C, M2, N2, LDC2, TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareMatrices(OType1C, M1, N1, LDC1, OType2C, M2, N2, LDC2, TOL);

      delete [] OType1A;
      delete [] OType1B;
      delete [] OType1C;
      delete [] OType2A;
      delete [] OType2B;
      delete [] OType2C;
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
      M1 = GetRandom(MVMIN, MVMAX);
      N1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      N2 = ConvertType( N1, N2 );
      SIDE = RandomSIDE();
      UPLO = RandomUPLO();

      LDA1 = GetRandom(MVMIN, MVMAX);
      if(Teuchos::ESideChar[SIDE] == 'L') {
	while (LDA1 < M1) { LDA1 = GetRandom(MVMIN, MVMAX); }
        LDA2 = ConvertType( LDA1, LDA2 );
	OType1A = new SType[LDA1 * M1];
	OType2A = new SType[LDA2 * M2];
	for(j1 = 0, j2 = 0; j1 < LDA1 * M1; j1++, j2++) {
	  OType1A[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2A[j2] = OType1A[j1];
	}
      } else {
	while (LDA1 < N1) { LDA1 = GetRandom(MVMIN, MVMAX); }
        LDA2 = ConvertType( LDA1, LDA2 );
	OType1A = new SType[LDA1 * N1];
	OType2A = new SType[LDA2 * N2];
	for(j1 = 0, j2 = 0; j1 < LDA1 * N1; j1++, j2++) {
	  OType1A[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2A[j2] = OType1A[j1];
	}
      }

      LDB1 = GetRandom(MVMIN, MVMAX);
      while (LDB1 < M1) {  LDB1 = GetRandom(MVMIN, MVMAX); }
      LDB2 = ConvertType( LDB1, LDB2 );
      OType1B = new SType[LDB1 * N1];
      OType2B = new SType[LDB2 * N2];
      for(j1 = 0, j2 = 0; j1 < LDB1 * N1; j1++, j2++) {
	  OType1B[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2B[j2] = OType1B[j1];
      }
    
      LDC1 = GetRandom(MVMIN, MVMAX);
      while (LDC1 < M1) {  LDC1 = GetRandom(MVMIN, MVMAX); }
      LDC2 = ConvertType( LDC1, LDC2 );
      OType1C = new SType[LDC1 * N1];
      OType2C = new SType[LDC2 * N2];
      for(j1 = 0, j2 = 0; j1 < LDC1 * N1; j1++, j2++) {
	  OType1C[j1] = GetRandom(-SCALARMAX, SCALARMAX);
	  OType2C[j2] = OType1C[j1];
      }
      
      OType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      OType1beta = GetRandom(-SCALARMAX, SCALARMAX);
      OType2alpha = OType1alpha;
      OType2beta = OType1beta;
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- SYMM -- " << std::endl;
	  std::cout << "SIDE = " << Teuchos::ESideChar[SIDE] << "\t" 
	       << "UPLO = " << Teuchos::EUploChar[UPLO] << std::endl;
	  std::cout << "OType1alpha = " << OType1alpha << std::endl;
	  std::cout << "OType2alpha = " << OType2alpha << std::endl;
	  std::cout << "OType1beta = " << OType1beta << std::endl;
	  std::cout << "OType2beta = " << OType2beta << std::endl;
	  if (Teuchos::ESideChar[SIDE] == 'L') {
	      PrintMatrix(OType1A, M1, M1, LDA1,"OType1A", matlab);
	      PrintMatrix(OType2A, M2, M2, LDA2,"OType2A", matlab);
	  } else {
	    PrintMatrix(OType1A, N1, N1, LDA1, "OType1A", matlab);
	    PrintMatrix(OType2A, N2, N2, LDA2, "OType2A", matlab);
	  }
	  PrintMatrix(OType1B, M1, N1, LDB1,"OType1B", matlab);
	  PrintMatrix(OType1C, M1, N1, LDC1,"OType1C_before_operation", matlab);
	  PrintMatrix(OType2B, M2, N2, LDB2,"OType2B", matlab);
	  PrintMatrix(OType2C, M2, N2, LDC2,"OType2C_before_operation", matlab);
	}
      TotalTestCount++;

      OType1BLAS.SYMM(SIDE, UPLO, M1, N1, OType1alpha, OType1A, LDA1, OType1B, LDB1, OType1beta, OType1C, LDC1);
      OType2BLAS.SYMM(SIDE, UPLO, M2, N2, OType2alpha, OType2A, LDA2, OType2B, LDB2, OType2beta, OType2C, LDC2);
      if(debug)
	{
	  PrintMatrix(OType1C, M1, N1, LDC1,"OType1C_after_operation", matlab);
	  PrintMatrix(OType2C, M2, N2, LDC2,"OType2C_after_operation", matlab);
	}
      if (CompareMatrices(OType1C, M1, N1, LDC1, OType2C, M2, N2, LDC2, TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareMatrices(OType1C, M1, N1, LDC1, OType2C, M2, N2, LDC2, TOL);

      delete [] OType1A;
      delete [] OType1B;
      delete [] OType1C;
      delete [] OType2A;
      delete [] OType2B;
      delete [] OType2C;
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
    N1 = GetRandom(MVMIN, MVMAX);
    K1 = GetRandom(MVMIN, MVMAX);
    while (K1 > N1) { K1 = GetRandom(MVMIN, MVMAX); }
    N2 = ConvertType( N1, N2 );
    K2 = ConvertType( K1, K2 );

    UPLO = RandomUPLO();
    TRANS = RandomTRANS();
#ifdef HAVE_TEUCHOS_COMPLEX
    while (TRANS == Teuchos::CONJ_TRANS) { TRANS = RandomTRANS(); }
#endif

    LDA1 = GetRandom(MVMIN, MVMAX);
    if(Teuchos::ETranspChar[TRANS] == 'N') {
      while (LDA1 < N1) { LDA1 = GetRandom(MVMIN, MVMAX); }
      LDA2 = ConvertType( LDA1, LDA2 );
      OType1A = new SType[LDA1 * K1];
      OType2A = new SType[LDA2 * K2];
      for(j1 = 0, j2 = 0; j1 < LDA1 * K1; j1++, j2++) {
        OType1A[j1] = GetRandom(-SCALARMAX, SCALARMAX);
        OType2A[j2] = OType1A[j1];
      }
    } else {
      while (LDA1 < K1) { LDA1 = GetRandom(MVMIN, MVMAX); }
      LDA2 = ConvertType( LDA1, LDA2 );
      OType1A = new SType[LDA1 * N1];
      OType2A = new SType[LDA2 * N2];
      for(j1 = 0, j2 = 0; j1 < LDA1 * N1; j1++, j2++) {
        OType1A[j1] = GetRandom(-SCALARMAX, SCALARMAX);
        OType2A[j2] = OType1A[j1];
      }
    }

    LDC1 = GetRandom(MVMIN, MVMAX);
    while (LDC1 < N1) {  LDC1 = GetRandom(MVMIN, MVMAX); }
    LDC2 = ConvertType( LDC1, LDC2 );
    OType1C = new SType[LDC1 * N1];
    OType2C = new SType[LDC2 * N2];
    for(j1 = 0, j2 = 0; j1 < N1; j1++, j2++) {

      if(Teuchos::EUploChar[UPLO] == 'U') {
        // The operator is upper triangular, make sure that the entries are
        // only in the upper triangular part of C.
        for(k1 = 0, k2 = 0; k1 < N1; k1++, k2++)
        {
          if(k1 <= j1) {
            OType1C[j1*LDC1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
          } else {
            OType1C[j1*LDC1+k1] = STypezero;
          } 
          OType2C[j2*LDC2+k2] = OType1C[j1*LDC1+k1];
        }
      }
      else {
        for(k1 = 0, k2 = 0; k1 < N1; k1++, k2++)
        {
          if(k1 >= j1) {
            OType1C[j1*LDC1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
          } else {
            OType1C[j1*LDC1+k1] = STypezero;
          }
          OType2C[j2*LDC2+k2] = OType1C[j1*LDC1+k1];
        }
      } 
    }
       
    OType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
    OType1beta = GetRandom(-SCALARMAX, SCALARMAX);
    OType2alpha = OType1alpha;
    OType2beta = OType1beta;
    if(debug)
    {
      std::cout << "Test #" << TotalTestCount << " -- SYRK -- " << std::endl;
      std::cout << "UPLO = " << Teuchos::EUploChar[UPLO] << "\t" 
 	        << "TRANS = " << Teuchos::ETranspChar[TRANS] << std::endl;
      std::cout << "N1="<<N1 << "\t" << "K1 = " << K1 
                << "\t" << "LDA1="<<LDA1 << "\t" << "LDC1=" << LDC1 << std::endl;
      std::cout << "N2="<<N2 << "\t" << "K2 = " << K2 
                << "\t" << "LDA2="<<LDA2 << "\t" << "LDC2=" << LDC2 << std::endl;
      std::cout << "OType1alpha = " << OType1alpha << std::endl;
      std::cout << "OType2alpha = " << OType2alpha << std::endl;
      std::cout << "OType1beta = " << OType1beta << std::endl;
      std::cout << "OType2beta = " << OType2beta << std::endl;
      if (Teuchos::ETranspChar[TRANS] == 'N') {
        PrintMatrix(OType1A, N1, K1, LDA1,"OType1A", matlab);
        PrintMatrix(OType2A, N2, K2, LDA2,"OType2A", matlab);
      } else {
        PrintMatrix(OType1A, K1, N1, LDA1, "OType1A", matlab);
        PrintMatrix(OType2A, K2, N2, LDA2, "OType2A", matlab);
      }
      PrintMatrix(OType1C, N1, N1, LDC1,"OType1C_before_operation", matlab);
      PrintMatrix(OType2C, N2, N2, LDC2,"OType2C_before_operation", matlab);
    }
    TotalTestCount++;

    OType1BLAS.SYRK(UPLO, TRANS, N1, K1, OType1alpha, OType1A, LDA1, OType1beta, OType1C, LDC1);
    OType2BLAS.SYRK(UPLO, TRANS, N2, K2, OType2alpha, OType2A, LDA2, OType2beta, OType2C, LDC2);
    if(debug)
    {
      PrintMatrix(OType1C, N1, N1, LDC1,"OType1C_after_operation", matlab);
      PrintMatrix(OType2C, N2, N2, LDC2,"OType2C_after_operation", matlab);
    }
    if (CompareMatrices(OType1C, N1, N1, LDC1, OType2C, N2, N2, LDC2, TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
    GoodTestSubcount += CompareMatrices(OType1C, N1, N1, LDC1, OType2C, N2, N2, LDC2, TOL);

    delete [] OType1A;
    delete [] OType1C;
    delete [] OType2A;
    delete [] OType2C;
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
      M1 = GetRandom(MVMIN, MVMAX);
      N1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      N2 = ConvertType( N1, N2 );

      LDB1 = GetRandom(MVMIN, MVMAX);
      while (LDB1 < M1) {
	  LDB1 = GetRandom(MVMIN, MVMAX);
      }
      LDB2 = ConvertType( LDB1, LDB2 );

      OType1B = new SType[LDB1 * N1];
      OType2B = new SType[LDB2 * N2];

      SIDE = RandomSIDE();
      UPLO = RandomUPLO();
      TRANSA = RandomTRANS();
      DIAG = RandomDIAG();

      if(Teuchos::ESideChar[SIDE] == 'L')  // The operator is on the left side
	{
          LDA1 = GetRandom(MVMIN, MVMAX);
      	  while (LDA1 < M1) {
	      LDA1 = GetRandom(MVMIN, MVMAX);
       	  }
          LDA2 = ConvertType( LDA1, LDA2 );

	  OType1A = new SType[LDA1 * M1];
	  OType2A = new SType[LDA2 * M2];

	  for(j1 = 0, j2 = 0; j1 < M1; j1++, j2++)
	    {	     
	      if(Teuchos::EUploChar[UPLO] == 'U') {
		// The operator is upper triangular, make sure that the entries are
		// only in the upper triangular part of A and the diagonal is non-zero.
		for(k1 = 0, k2 = 0; k1 < M1; k1++, k2++) 
		{
		    if(k1 < j1) {
	      		OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			OType1A[j1*LDA1+k1] = STypezero;
		    }
	      	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    if(k1 == j1) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
	      		    OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (OType1A[j1*LDA1+k1] == STypezero) {
				OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];	
		    	} else {
	      		    OType1A[j1*LDA1+k1] = STypeone;
	      	    	    OType2A[j2*LDA2+k2] = STypeone;
		    	}			
		    }
		}
	      } else {
		// The operator is lower triangular, make sure that the entries are
		// only in the lower triangular part of A and the diagonal is non-zero.
		for(k1 = 0, k2 = 0; k1 < M1; k1++, k2++) 
		{
		    if(k1 > j1) {
	      		OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			OType1A[j1*LDA1+k1] = STypezero;
		    }
	      	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    if(k1 == j1) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
      			    OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (OType1A[j1*LDA1+k1] == STypezero) {
				OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    	} else {
	      		    OType1A[j1*LDA1+k1] = STypeone;
	      	    	    OType2A[j2*LDA2+k2] = STypeone;
			}
		    }			
		} // end for(k=0 ...		
	      } // end if(UPLO == 'U') ...
	    } // end for(j=0 ...
	} // if(SIDE == 'L') ...
      else // The operator is on the right side
	{
          LDA1 = GetRandom(MVMIN, MVMAX);
      	  while (LDA1 < N1) {
	      LDA1 = GetRandom(MVMIN, MVMAX);
       	  }
          LDA2 = ConvertType( LDA1, LDA2 );

	  OType1A = new SType[LDA1 * N1];
	  OType2A = new SType[LDA2 * N2];

	  for(j1 = 0, j2 = 0; j1 < N1; j1++, j2++)
	    {	     
	      if(Teuchos::EUploChar[UPLO] == 'U') {
		// The operator is upper triangular, make sure that the entries are
		// only in the upper triangular part of A and the diagonal is non-zero.
		for(k1 = 0, k2 = 0; k1 < N1; k1++, k2++) 
		{
		    if(k1 < j1) {
	      		OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			OType1A[j1*LDA1+k1] = STypezero;
		    }
	      	    OType2A[j1*LDA1+k1] = OType1A[j1*LDA1+k1];
		    if(k1 == j1) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
	      		    OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (OType1A[j1*LDA1+k1] == STypezero) {
				OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    	} else {
	      		    OType1A[j1*LDA1+k1] = STypeone;
	      	    	    OType2A[j2*LDA2+k2] = STypeone;
			}
		    }			
		}
	      } else {
		// The operator is lower triangular, make sure that the entries are
		// only in the lower triangular part of A and the diagonal is non-zero.
		for(k1 = 0, k2 = 0; k1 < N1; k1++, k2++) 
		{
		    if(k1 > j1) {
	      		OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			OType1A[j1*LDA1+k1] = STypezero;
		    }
	      	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    if(k1 == j1) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
	      		    OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (OType1A[j1*LDA1+k1] == STypezero) {
				OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    	} else {
	      		    OType1A[j1*LDA1+k1] = STypeone;
	      	    	    OType2A[j2*LDA2+k2] = STypeone;
			}
		    }			
		} // end for(k=0 ...		
	      } // end if(UPLO == 'U') ...
	    } // end for(j=0 ...
	} // end if(SIDE == 'L') ...

      // Fill in the right hand side block B.
      for(j1 = 0, j2 = 0; j1 < N1; j1++, j2++) {
	  for(k1 = 0, k2 = 0; k1 < M1; k1++, k2++) {
	    OType1B[j1*LDB1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
	    OType2B[j2*LDB2+k2] = OType1B[j1*LDB1+k1];
	  }
      }
      OType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      OType2alpha = OType1alpha;
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- TRMM -- " << std::endl;
	  std::cout << "SIDE = " << Teuchos::ESideChar[SIDE] << "\t" 
	       << "UPLO = " << Teuchos::EUploChar[UPLO] << "\t" 
	       << "TRANSA = " << Teuchos::ETranspChar[TRANSA] << "\t" 
	       << "DIAG = " << Teuchos::EDiagChar[DIAG] << std::endl;
	  std::cout << "OType1alpha = " << OType1alpha << std::endl;
	  std::cout << "OType2alpha = " << OType2alpha << std::endl;
          if(Teuchos::ESideChar[SIDE] == 'L') { 
	    PrintMatrix(OType1A, M1, M1, LDA1, "OType1A", matlab);
	    PrintMatrix(OType2A, M2, M2, LDA2, "OType2A", matlab);
	  } else {
	    PrintMatrix(OType1A, N1, N1, LDA1, "OType1A", matlab);
	    PrintMatrix(OType2A, N2, N2, LDA2, "OType2A", matlab);
	  }
	  PrintMatrix(OType1B, M1, N1, LDB1,"OType1B_before_operation", matlab);
	  PrintMatrix(OType2B, M2, N2, LDB2,"OType2B_before_operation", matlab);
	}
      TotalTestCount++;
      OType1BLAS.TRMM(SIDE, UPLO, TRANSA, DIAG, M1, N1, OType1alpha, OType1A, LDA1, OType1B, LDB1);
      OType2BLAS.TRMM(SIDE, UPLO, TRANSA, DIAG, M2, N2, OType2alpha, OType2A, LDA2, OType2B, LDB2);
      if(debug)
	{
	  PrintMatrix(OType1B, M1, N1, LDB1, "OType1B_after_operation", matlab);
	  PrintMatrix(OType2B, M2, N2, LDB2, "OType2B_after_operation", matlab);
	}
      if (CompareMatrices(OType1B, M1, N1, LDB1, OType2B, M2, N2, LDB2, TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareMatrices(OType1B, M1, N1, LDB1, OType2B, M2, N2, LDB2, TOL);
      delete [] OType1A;
      delete [] OType1B;
      delete [] OType2A;
      delete [] OType2B;
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
      M1 = GetRandom(MVMIN, MVMAX);
      N1 = GetRandom(MVMIN, MVMAX);
      M2 = ConvertType( M1, M2 );
      N2 = ConvertType( N1, N2 );

      LDB1 = GetRandom(MVMIN, MVMAX);
      while (LDB1 < M1) {
	  LDB1 = GetRandom(MVMIN, MVMAX);
      }
      LDB2 = ConvertType( LDB1, LDB2 );

      OType1B = new SType[LDB1 * N1];
      OType2B = new SType[LDB2 * N2];

      SIDE = RandomSIDE();
      UPLO = RandomUPLO();
      TRANSA = RandomTRANS();
      // Since the entries are integers, we don't want to use the unit diagonal feature,
      // this creates ill-conditioned, nearly-singular matrices.
      //DIAG = RandomDIAG();  
      DIAG = Teuchos::NON_UNIT_DIAG;

      if(Teuchos::ESideChar[SIDE] == 'L')  // The operator is on the left side
	{
          LDA1 = GetRandom(MVMIN, MVMAX);
      	  while (LDA1 < M1) {
	      LDA1 = GetRandom(MVMIN, MVMAX);
       	  }
          LDA2 = ConvertType( LDA1, LDA2 );

	  OType1A = new SType[LDA1 * M1];
	  OType2A = new SType[LDA2 * M2];

	  for(j1 = 0, j2 = 0; j1 < M1; j1++, j2++)
	    {	     
	      if(Teuchos::EUploChar[UPLO] == 'U') {
		// The operator is upper triangular, make sure that the entries are
		// only in the upper triangular part of A and the diagonal is non-zero.
		for(k1 = 0, k2 = 0; k1 < M1; k1++, k2++) 
		{
		    if(k1 < j1) {
	      		OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			OType1A[j1*LDA1+k1] = STypezero;
		    }
	      	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    if(k1 == j1) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
	      		    OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (OType1A[j1*LDA1+k1] == STypezero) {
				OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];	
		    	} else {
	      		    OType1A[j1*LDA1+k1] = STypeone;
	      	    	    OType2A[j2*LDA2+k2] = STypeone;
		    	}			
		    }
		}
	      } else {
		// The operator is lower triangular, make sure that the entries are
		// only in the lower triangular part of A and the diagonal is non-zero.
		for(k1 = 0, k2 = 0; k1 < M1; k1++, k2++) 
		{
		    if(k1 > j1) {
	      		OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			OType1A[j1*LDA1+k1] = STypezero;
		    }
	      	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    if(k1 == j1) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
      			    OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (OType1A[j1*LDA1+k1] == STypezero) {
				OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    	} else {
	      		    OType1A[j1*LDA1+k1] = STypeone;
	      	    	    OType2A[j2*LDA2+k2] = STypeone;
			}
		    }			
		} // end for(k=0 ...		
	      } // end if(UPLO == 'U') ...
	    } // end for(j=0 ...
	} // if(SIDE == 'L') ...
      else // The operator is on the right side
	{
          LDA1 = GetRandom(MVMIN, MVMAX);
      	  while (LDA1 < N1) {
	      LDA1 = GetRandom(MVMIN, MVMAX);
       	  }
          LDA2 = ConvertType( LDA1, LDA2 );

	  OType1A = new SType[LDA1 * N1];
	  OType2A = new SType[LDA2 * N2];

	  for(j1 = 0, j2 = 0; j1 < N1; j1++, j2++)
	    {	     
	      if(Teuchos::EUploChar[UPLO] == 'U') {
		// The operator is upper triangular, make sure that the entries are
		// only in the upper triangular part of A and the diagonal is non-zero.
		for(k1 = 0, k2 = 0; k1 < N1; k1++, k2++) 
		{
		    if(k1 < j1) {
	      		OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			OType1A[j1*LDA1+k1] = STypezero;
		    }
	      	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    if(k1 == j1) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
	      		    OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (OType1A[j1*LDA1+k1] == STypezero) {
				OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    	} else {
	      		    OType1A[j1*LDA1+k1] = STypeone;
	      	    	    OType2A[j2*LDA2+k2] = STypeone;
			}
		    }			
		}
	      } else {
		// The operator is lower triangular, make sure that the entries are
		// only in the lower triangular part of A and the diagonal is non-zero.
		for(k1 = 0, k2 = 0; k1 < N1; k1++, k2++) 
		{
		    if(k1 > j1) {
	      		OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
		    } else {
			OType1A[j1*LDA1+k1] = STypezero;
		    }
	      	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    if(k1 == j1) {
			if (Teuchos::EDiagChar[DIAG] == 'N') {
	      		    OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    while (OType1A[j1*LDA1+k1] == STypezero) {
				OType1A[j1*LDA1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
			    }
	      	    	    OType2A[j2*LDA2+k2] = OType1A[j1*LDA1+k1];
		    	} else {
	      		    OType1A[j1*LDA1+k1] = STypeone;
	      	    	    OType2A[j2*LDA2+k2] = STypeone;
			}
		    }			
		} // end for(k=0 ...		
	      } // end if(UPLO == 'U') ...
	    } // end for(j=0 ...
	} // end if(SIDE == 'L') ...

      // Fill in the right hand side block B.
      for(j1 = 0, j2 = 0; j1 < N1; j1++, j2++)
	{
	  for(k1 = 0, k2 = 0; k1 < M1; k1++, k2++) 
	    {
	  	OType1B[j1*LDB1+k1] = GetRandom(-SCALARMAX, SCALARMAX);
	  	OType2B[j2*LDB2+k2] = OType1B[j1*LDB1+k1];
	    }
	}
      
      OType1alpha = GetRandom(-SCALARMAX, SCALARMAX);
      OType2alpha = OType1alpha;
      
      if(debug)
	{
	  std::cout << "Test #" << TotalTestCount << " -- TRSM -- " << std::endl;
	  std::cout << "SIDE = " << Teuchos::ESideChar[SIDE] << "\t" 
	       << "UPLO = " << Teuchos::EUploChar[UPLO] << "\t" 
	       << "TRANSA = " << Teuchos::ETranspChar[TRANSA] << "\t" 
	       << "DIAG = " << Teuchos::EDiagChar[DIAG] << std::endl;
	  std::cout << "M1="<<M1 << "\t" << "N1="<<N1 << "\t" << "LDA1="<<LDA1 << "\t" << "LDB1="<<LDB1 << std::endl;
	  std::cout << "M2="<<M2 << "\t" << "N2="<<N2 << "\t" << "LDA2="<<LDA2 << "\t" << "LDB2="<<LDB2 << std::endl;
	  std::cout << "OType1alpha = " << OType1alpha << std::endl;
	  std::cout << "OType2alpha = " << OType2alpha << std::endl;
	  if (Teuchos::ESideChar[SIDE] == 'L') {
	      PrintMatrix(OType1A, M1, M1, LDA1, "OType1A", matlab);
	      PrintMatrix(OType2A, M2, M2, LDA2, "OType2A", matlab);
	  } else {
	      PrintMatrix(OType1A, N1, N1, LDA1, "OType1A", matlab);
	      PrintMatrix(OType2A, N2, N2, LDA2, "OType2A", matlab);
	  }
	  PrintMatrix(OType1B, M1, N1, LDB1, "OType1B_before_operation", matlab);
	  PrintMatrix(OType2B, M2, N2, LDB2, "OType2B_before_operation", matlab);
	}
      TotalTestCount++;

      OType1BLAS.TRSM(SIDE, UPLO, TRANSA, DIAG, M1, N1, OType1alpha, OType1A, LDA1, OType1B, LDB1);
      OType2BLAS.TRSM(SIDE, UPLO, TRANSA, DIAG, M2, N2, OType2alpha, OType2A, LDA2, OType2B, LDB2);
 
      if(debug)
	{
	  PrintMatrix(OType1B, M1, N1, LDB1, "OType1B_after_operation", matlab);
	  PrintMatrix(OType2B, M2, N2, LDB2, "OType2B_after_operation", matlab);
	}

      if (CompareMatrices(OType1B, M1, N1, LDB1, OType2B, M2, N2, LDB2, TOL)==0)
	std::cout << "FAILED TEST!!!!!!" << std::endl;
      GoodTestSubcount += CompareMatrices(OType1B, M1, N1, LDB1, OType2B, M2, N2, LDB2, TOL);

      delete [] OType1A;
      delete [] OType1B;
      delete [] OType2A;
      delete [] OType2B;
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

template<typename T>
std::complex<T> GetRandom( std::complex<T> Low, std::complex<T> High)
{
  T lowMag = Low.real();
  T highMag = High.real();
  T real = (T)(((1.0 * ScalarTraits<int>::random()) / RAND_MAX) * (highMag - lowMag + ScalarTraits<T>::one())) + lowMag;
  T imag = (T)(((1.0 * ScalarTraits<int>::random()) / RAND_MAX) * (highMag - lowMag + ScalarTraits<T>::one())) + lowMag;
  return std::complex<T>( real, imag );
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

template<typename TYPE, typename OTYPE>
void PrintVector(TYPE* Vector, OTYPE Size, std::string Name, bool Matlab)
{
  std::cout << Name << " =" << std::endl;
  OTYPE i;
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

template<typename TYPE, typename OTYPE>
void PrintMatrix(TYPE* Matrix, OTYPE Rows, OTYPE Columns, OTYPE LDM, std::string Name, bool Matlab)
{
  if(!Matlab)
    {
      std::cout << Name << " =" << std::endl;
      OTYPE i, j;
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
      OTYPE i, j;
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

template<typename TYPE>
bool CompareScalars(TYPE Scalar1, TYPE Scalar2, typename ScalarTraits<TYPE>::magnitudeType Tolerance)
{
  typename ScalarTraits<TYPE>::magnitudeType temp = ScalarTraits<TYPE>::magnitude(Scalar2);
  typename ScalarTraits<TYPE>::magnitudeType temp2 = ScalarTraits<TYPE>::magnitude(Scalar1 - Scalar2);
  if (temp != ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::zero()) {
    temp2 /= temp;
  }
  return( temp2 < Tolerance );
}


/*  Function:  CompareVectors
    Purpose:   Compares the difference between two vectors using relative euclidean-norms, i.e. ||v_1-v_2||_2/||v_2||_2
*/
template<typename TYPE, typename OTYPE1, typename OTYPE2>
bool CompareVectors(TYPE* Vector1, OTYPE1 Size1, TYPE* Vector2, OTYPE2 Size2, typename ScalarTraits<TYPE>::magnitudeType Tolerance)
{
  TYPE temp = ScalarTraits<TYPE>::zero();
  typename ScalarTraits<TYPE>::magnitudeType temp2 = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::zero();
  typename ScalarTraits<TYPE>::magnitudeType temp3 = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::zero();
  typename ScalarTraits<TYPE>::magnitudeType sum = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::zero();
  typename ScalarTraits<TYPE>::magnitudeType sum2 = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::zero();
  OTYPE1 i1;
  OTYPE2 i2;
  for(i1 = 0, i2 = 0; i1 < Size1; i1++, i2++)
    {
      sum2 += ScalarTraits<TYPE>::magnitude(ScalarTraits<TYPE>::conjugate(Vector2[i2])*Vector2[i2]);
      temp = Vector1[i1] - Vector2[i2];
      sum += ScalarTraits<TYPE>::magnitude(ScalarTraits<TYPE>::conjugate(temp)*temp);
    }
  temp2 = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::squareroot(sum2);
  if (temp2 != ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::zero())
    temp3 = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::squareroot(sum)/temp2;
  else
    temp3 = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::squareroot(sum);
  if (temp3 > Tolerance )
    return false;
  else
    return true;
}

/*  Function:  CompareMatrices
    Purpose:   Compares the difference between two matrices using relative frobenius-norms, i.e. ||M_1-M_2||_F/||M_2||_F
*/
template<typename TYPE, typename OTYPE1, typename OTYPE2>
bool CompareMatrices(TYPE* Matrix1, OTYPE1 Rows1, OTYPE1 Columns1, OTYPE1 LDM1,
                     TYPE* Matrix2, OTYPE2 Rows2, OTYPE2 Columns2, OTYPE2 LDM2, 
                     typename ScalarTraits<TYPE>::magnitudeType Tolerance)
{
  TYPE temp = ScalarTraits<TYPE>::zero();
  typename ScalarTraits<TYPE>::magnitudeType temp2 = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::zero();
  typename ScalarTraits<TYPE>::magnitudeType temp3 = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::zero();
  typename ScalarTraits<TYPE>::magnitudeType sum = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::zero();
  typename ScalarTraits<TYPE>::magnitudeType sum2 = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::zero();
  OTYPE1 i1, j1;
  OTYPE2 i2, j2;
  for(j1 = 0, j2 = 0; j1 < Columns1; j1++, j2++)
    {
      for(i1 = 0, i2 = 0; i1 < Rows1; i1++, i2++)
	{
	  sum2 = ScalarTraits<TYPE>::magnitude(ScalarTraits<TYPE>::conjugate(Matrix2[j2*LDM2 + i2])*Matrix2[j2*LDM2 + i2]);
	  temp = Matrix1[j1*LDM1 + i1] - Matrix2[j2*LDM2 + i2]; 
	  sum = ScalarTraits<TYPE>::magnitude(ScalarTraits<TYPE>::conjugate(temp)*temp);
	}
    }
  temp2 = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::squareroot(sum2);
  if (temp2 != ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::zero())
    temp3 = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::squareroot(sum)/temp2;
  else
    temp3 = ScalarTraits<typename ScalarTraits<TYPE>::magnitudeType>::squareroot(sum);
  if (temp3 > Tolerance)
    return false;
  else
    return true;
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

