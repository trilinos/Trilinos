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

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "Teuchos_SerialSpdDenseSolver.hpp"
#include "Teuchos_SerialQRDenseSolver.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Version.hpp"

using Teuchos::ScalarTraits;
using Teuchos::SerialDenseMatrix;
using Teuchos::SerialSymDenseMatrix;
using Teuchos::SerialDenseVector;

#define OTYPE int
#ifdef HAVE_TEUCHOS_COMPLEX
#define STYPE std::complex<double>
#else
#define STYPE double
#endif

// SCALARMAX defines the maximum positive value (with a little leeway) generated for matrix and vector elements and scalars:
// random numbers in [-SCALARMAX, SCALARMAX] will be generated.
#ifdef HAVE_TEUCHOS_COMPLEX
#define SCALARMAX  STYPE(10,0)
#else
#define SCALARMAX  STYPE(10)
#endif

template<typename TYPE>
int PrintTestResults(std::string, TYPE, TYPE, bool);

int ReturnCodeCheck(std::string, int, int, bool);

typedef SerialDenseVector<OTYPE, STYPE> DVector;
typedef SerialDenseMatrix<OTYPE, STYPE> DMatrix;
typedef SerialSymDenseMatrix<OTYPE, STYPE> SDMatrix;

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

// Generates random matrices and vectors using GetRandom()
Teuchos::RCP<SDMatrix> GetRandomSpdMatrix(int n);
Teuchos::RCP<DMatrix> GetRandomMatrix(int m, int n);
Teuchos::RCP<DVector> GetRandomVector(int n);

// Compares the difference between two vectors using relative euclidean norms
// Returns 1 if the comparison failed, the relative difference is greater than the tolerance.
int CompareVectors(const SerialDenseVector<OTYPE,STYPE>& Vector1,
		   const SerialDenseVector<OTYPE,STYPE>& Vector2,
		   ScalarTraits<STYPE>::magnitudeType Tolerance );

int main(int argc, char* argv[])
{
  typedef ScalarTraits<STYPE>::magnitudeType MagnitudeType;

  int n=10;
  MagnitudeType tol = 1e-12*ScalarTraits<MagnitudeType>::one();

  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (verbose)
    std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  int numberFailedTests = 0, failedComparison = 0;
  int returnCode = 0;
  std::string testName = "", testType = "";

#ifdef HAVE_TEUCHOS_COMPLEX
  testType = "COMPLEX";
#else
  testType = "REAL";
#endif

  if (verbose) std::cout<<std::endl<<"********** CHECKING TEUCHOS SERIAL SPD DENSE SOLVER - " << testType << "-VALUED **********"<<std::endl<<std::endl;

  // Create dense matrix and vector.
  Teuchos::RCP<SDMatrix> A1 = GetRandomSpdMatrix(n);
  Teuchos::RCP<DVector> x1 = GetRandomVector(n);
  DVector xhat(n), b(n), bt(n);

  Teuchos::RCP<DMatrix> A1full = Teuchos::rcp( new DMatrix(n,n) );
  for (int j=0; j<n; ++j) {
    for (int i=0; i<=j; ++i) {
      (*A1full)(i,j) = (*A1)(i,j);
    }
  }
  for (int j=0; j<n; ++j) {
    for (int i=j+1; i<n; ++i) {
      (*A1full)(i,j) = ScalarTraits<STYPE>::conjugate( (*A1)(i,j) );
    }
  }

  // Compute the right-hand side vector using multiplication.
  returnCode = b.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, ScalarTraits<STYPE>::one() , *A1full, *x1, ScalarTraits<STYPE>::zero());
  testName = "Generating right-hand side vector using A*x, where x is a random vector:";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);

  // Fill the solution vector with zeros.
  xhat.putScalar( ScalarTraits<STYPE>::zero() );

  // Create a serial spd dense solver.
  Teuchos::SerialSpdDenseSolver<OTYPE, STYPE> solver1;

  // Pass in matrix and vectors
  solver1.setMatrix( A1 );
  solver1.setVectors( Teuchos::rcp( &xhat, false ), Teuchos::rcp( &b, false ) );

  // Test1:  Simple factor and solve
  returnCode = solver1.factor();
  testName = "Simple solve: factor() random A:";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);

  // Non-transpose solve
  returnCode = solver1.solve();
  testName = "Simple solve: solve() random A (NO_TRANS):";
  failedComparison = CompareVectors( *x1, xhat, tol );
  if(verbose && failedComparison>0) std::cout << "COMPARISON FAILED : ";
  numberFailedTests += failedComparison;
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);

  // Test2:  Solve with iterative refinement.
#ifdef HAVE_TEUCHOSNUMERICS_EIGEN
  // Iterative refinement not implemented in Eigen
#else
  // Create random linear system
  Teuchos::RCP<SDMatrix> A2 = GetRandomSpdMatrix(n);
  Teuchos::RCP<DVector> x2 = GetRandomVector(n);

  Teuchos::RCP<DMatrix> A2full = Teuchos::rcp( new DMatrix(n,n) );
  for (int j=0; j<n; ++j) {
    for (int i=0; i<=j; ++i) {
      (*A2full)(i,j) = (*A2)(i,j);
    }
  }
  for (int j=0; j<n; ++j) {
    for (int i=j+1; i<n; ++i) {
      (*A2full)(i,j) = ScalarTraits<STYPE>::conjugate( (*A2)(i,j) );
    }
  }

  // Create LHS through multiplication with A2
  xhat.putScalar( ScalarTraits<STYPE>::zero() );
  b.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, ScalarTraits<STYPE>::one() , *A2full, *x2, ScalarTraits<STYPE>::zero());

  // Create a serial dense solver.
  Teuchos::SerialSpdDenseSolver<OTYPE, STYPE> solver2;
  solver2.solveToRefinedSolution( true );

  // Pass in matrix and vectors
  solver2.setMatrix( A2 );
  solver2.setVectors( Teuchos::rcp( &xhat, false ), Teuchos::rcp( &b, false ) );

  // Factor and solve with iterative refinement.
  returnCode = solver2.factor();
  testName = "Solve with iterative refinement: factor() random A:";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);

  // Non-transpose solve
  returnCode = solver2.solve();
  testName = "Solve with iterative refinement: solve() random A (NO_TRANS):";
  failedComparison = CompareVectors( *x2, xhat, tol );
  if(verbose && failedComparison>0) std::cout << "COMPARISON FAILED : ";
  numberFailedTests += failedComparison;
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);

#endif

  // Test3:  Solve with matrix equilibration.

  // Create random linear system
  Teuchos::RCP<SDMatrix> A3 = GetRandomSpdMatrix(n);
  Teuchos::RCP<DVector> x3 = GetRandomVector(n);

  Teuchos::RCP<DMatrix> A3full = Teuchos::rcp( new DMatrix(n,n) );
  for (int j=0; j<n; ++j) {
    for (int i=0; i<=j; ++i) {
      (*A3full)(i,j) = (*A3)(i,j);
    }
  }
  for (int j=0; j<n; ++j) {
    for (int i=j+1; i<n; ++i) {
      (*A3full)(i,j) = ScalarTraits<STYPE>::conjugate( (*A3)(i,j) );
    }
  }

  // Create LHS through multiplication with A3
  xhat.putScalar( ScalarTraits<STYPE>::zero() );
  b.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, ScalarTraits<STYPE>::one() , *A3full, *x3, ScalarTraits<STYPE>::zero());

  Teuchos::RCP<SDMatrix> A3bak = Teuchos::rcp( new SDMatrix( *A3 ) );
  Teuchos::RCP<DVector> b3bak = Teuchos::rcp( new DVector( Teuchos::Copy, b ) );

  // Create a serial dense solver.
  Teuchos::SerialSpdDenseSolver<OTYPE, STYPE> solver3;
  solver3.factorWithEquilibration( true );

  // Pass in matrix and vectors
  solver3.setMatrix( A3 );
  solver3.setVectors( Teuchos::rcp( &xhat, false ), Teuchos::rcp( &b, false ) );

  // Factor and solve with matrix equilibration.
  returnCode = solver3.factor();
  testName = "Solve with matrix equilibration: factor() random A:";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);

  // Non-transpose solve
  returnCode = solver3.solve();
  testName = "Solve with matrix equilibration: solve() random A (NO_TRANS):";
  failedComparison = CompareVectors( *x3, xhat, tol );
  if(verbose && failedComparison>0) std::cout << "COMPARISON FAILED : ";
  numberFailedTests += failedComparison;
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);

  // Factor and solve with matrix equilibration, where only solve is called.
  solver3.setMatrix( A3bak );
  solver3.setVectors( Teuchos::rcp( &xhat, false ), b3bak );
  returnCode = solver3.solve();
  testName = "Solve with matrix equilibration: solve() without factor() random A (NO_TRANS):";
  failedComparison = CompareVectors( *x3, xhat, tol );
  if(verbose && failedComparison>0) std::cout << "COMPARISON FAILED : ";
  numberFailedTests += failedComparison;
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);

  //
  // If a test failed output the number of failed tests.
  //
  if(numberFailedTests > 0)
  {
	    if (verbose) {
		std::cout << "Number of failed tests: " << numberFailedTests << std::endl;
		std::cout << "End Result: TEST FAILED" << std::endl;
		return -1;
	    }
	}
  if(numberFailedTests == 0)
    std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;

}

template<typename TYPE>
int PrintTestResults(std::string testName, TYPE calculatedResult, TYPE expectedResult, bool verbose)
{
  int result;
  if(calculatedResult == expectedResult)
    {
      if(verbose) std::cout << testName << " successful." << std::endl;
      result = 0;
    }
  else
    {
      if(verbose) std::cout << testName << " unsuccessful." << std::endl;
      result = 1;
    }
  return result;
}

int ReturnCodeCheck(std::string testName, int returnCode, int expectedResult, bool verbose)
{
  int result;
  if(expectedResult == 0)
    {
      if(returnCode == 0)
	{
	  if(verbose) std::cout << testName << " test successful." << std::endl;
	  result = 0;
	}
      else
	{
	  if(verbose) std::cout << testName << " test unsuccessful. Return code was " << returnCode << "." << std::endl;
	  result = 1;
	}
    }
  else
    {
      if(returnCode != 0)
	{
	  if(verbose) std::cout << testName << " test successful -- failed as expected." << std::endl;
	  result = 0;
	}
      else
	{
	  if(verbose) std::cout << testName << " test unsuccessful -- did not fail as expected. Return code was " << returnCode << "." << std::endl;
	  result = 1;
	}
    }
  return result;
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

Teuchos::RCP<SDMatrix> GetRandomSpdMatrix(int n)
{

  Teuchos::RCP<SDMatrix> mat = Teuchos::rcp( new SDMatrix(n) );
  Teuchos::RCP<DMatrix> tmp = Teuchos::rcp( new DMatrix(n,n) );
  Teuchos::RCP<DMatrix> tmp2 = Teuchos::rcp( new DMatrix(n,n) );
  Teuchos::RCP<DMatrix> Q = Teuchos::rcp( new DMatrix(n,n) );
  Teuchos::RCP<DMatrix> D = Teuchos::rcp( new DMatrix(n,n) );

  // QR factor a random matrix and get the orthogonal Q
  Teuchos::RCP<DMatrix> A = GetRandomMatrix(n,n);
  Teuchos::SerialQRDenseSolver<OTYPE, STYPE> solver;
  solver.setMatrix( A );
  solver.factor();
  solver.formQ();
  Q = solver.getQ();

  // Get n random positive eigenvalues and put them in a diagonal matrix D
  Teuchos::RCP<DVector> E = GetRandomVector(n);
  for (int i=0; i<n; ++i) {
    (*D)(i,i) = ScalarTraits<STYPE>::magnitude(ScalarTraits<STYPE>::one()) +
                ScalarTraits<STYPE>::magnitude( (*E)(i) ) * ScalarTraits<STYPE>::magnitude( (*E)(i) );
  }

  // Form the spd matrix Q' D Q
  tmp->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, ScalarTraits<STYPE>::one() , *D, *Q, ScalarTraits<STYPE>::zero());
  tmp2->multiply(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, ScalarTraits<STYPE>::one() , *Q, *tmp, ScalarTraits<STYPE>::zero());

  mat->setUpper();
  for (int j=0; j<n; ++j) {
    for (int i=0; i<=j; ++i) {
      (*mat)(i,j) = (*tmp2)(i,j);
    }
  }

  return mat;
}

Teuchos::RCP<DMatrix> GetRandomMatrix(int m, int n)
{
  Teuchos::RCP<DMatrix> newmat = Teuchos::rcp( new DMatrix(m,n) );

  // Fill dense matrix with random entries.
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      (*newmat)(i,j) = GetRandom(-SCALARMAX, SCALARMAX);

  return newmat;
}

Teuchos::RCP<DVector> GetRandomVector(int n)
{
  Teuchos::RCP<DVector> newvec = Teuchos::rcp( new DVector( n ) );

  // Fill dense vector with random entries.
  for (int i=0; i<n; i++)
    (*newvec)(i) = GetRandom(-SCALARMAX, SCALARMAX);

  return newvec;
}

/*  Function:  CompareVectors
    Purpose:   Compares the difference between two vectors using relative euclidean-norms, i.e. ||v_1-v_2||_2/||v_2||_2
*/
int CompareVectors(const SerialDenseVector<OTYPE,STYPE>& Vector1,
		   const SerialDenseVector<OTYPE,STYPE>& Vector2,
		   ScalarTraits<STYPE>::magnitudeType Tolerance )
{
  typedef ScalarTraits<STYPE>::magnitudeType MagnitudeType;

  SerialDenseVector<OTYPE,STYPE> diff( Vector1 );
  diff -= Vector2;

  MagnitudeType norm_diff = diff.normFrobenius();
  MagnitudeType norm_v2 = Vector2.normFrobenius();
  MagnitudeType temp = norm_diff;
  if (norm_v2 != ScalarTraits<MagnitudeType>::zero())
    temp /= norm_v2;

  if (temp > Tolerance)
    return 1;
  else
    return 0;
}
