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

#include "Teuchos_SerialBandDenseMatrix.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "Teuchos_Version.hpp"

#define OTYPE int
#define STYPE std::complex<double>

template<typename TYPE>
int PrintTestResults(std::string, TYPE, TYPE, bool);

int ReturnCodeCheck(std::string, int, int, bool);

typedef Teuchos::SerialBandDenseMatrix<OTYPE, STYPE> BDMatrix;
typedef Teuchos::SerialDenseMatrix<OTYPE, STYPE> DMatrix;
typedef Teuchos::SerialDenseVector<OTYPE, STYPE> DVector;

int main(int argc, char* argv[])
{

  int i;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (verbose)
    std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  int numberFailedTests = 0;
  int returnCode = 0;
  std::string testName = "";

  if (verbose) std::cout<<std::endl<<"********** CHECKING TEUCHOS SERIAL BANDED DENSE MATRIX **********"<<std::endl<<std::endl;

  // default constructor test
  BDMatrix DefConTest;
  if (verbose) std::cout <<"default constructor -- construct empty matrix ";
  if ( DefConTest.values()!=NULL || DefConTest.numCols()!=0 || DefConTest.numRows()!=0 ||DefConTest.stride()!=0 ||DefConTest.empty()!=true ) {
	if (verbose) std::cout << "unsuccessful."<<std::endl;
	numberFailedTests++;
  } else {
	if (verbose) std::cout << "successful."<<std::endl;
  }

  // constructor 1 (matrix w/ dimension but empty)

  BDMatrix Con1Test( 4, 4, 1, 1 );
  if (verbose) std::cout <<"constructor 1 -- empty matrix with given dimensions ";
  if ( Con1Test.numRows()!=4 || Con1Test.numCols()!=4 || Con1Test( 1, 2 )!=0.0 ) {
	if (verbose) std::cout << "unsuccessful."<<std::endl;
	numberFailedTests++;
  } else {
	if (verbose) std::cout << "successful."<<std::endl;
  }

  // constructor 2 (from array) tests

  STYPE a[25];
	      a[4] = 5;  a[8]  = 11;   a[12] = 17;  a[16] = 23;
  a[1] = 0;   a[5] = 6;  a[9]  = 12;   a[13] = 18;  a[17] = 24;
  a[2] = 1;   a[6] = 7;  a[10] = 13;   a[14] = 19;
  a[3] = 2;   a[7] = 8;  a[11] = 14;

  BDMatrix C2T1ER;
  C2T1ER.shape(5, 5, 2, 1);
  C2T1ER(0, 0) = 0;  C2T1ER(0, 1) = 5;
  C2T1ER(1, 0) = 1;  C2T1ER(1, 1) = 6;  C2T1ER(1, 2) = 11;
  C2T1ER(2, 0) = 2;  C2T1ER(2, 1) = 7;  C2T1ER(2, 2) = 12;  C2T1ER(2, 3) = 17;
		     C2T1ER(3, 1) = 8;  C2T1ER(3, 2) = 13;  C2T1ER(3, 3) = 18;  C2T1ER(3, 4) = 23;
					C2T1ER(4, 2) = 14;  C2T1ER(4, 3) = 19;  C2T1ER(4, 4) = 24;

  // Create another lower triangular matrix with a view of 'a'.
  BDMatrix Con2Test1(Teuchos::Copy, a, 4, 5, 5, 2, 1);
  numberFailedTests += PrintTestResults("constructor 2 -- construct matrix from array subrange", Con2Test1, C2T1ER, verbose);

  // constructor 3 (copy constructor)

  BDMatrix Con3TestCopy( C2T1ER );
  if(verbose) std::cout <<"constructor 3 -- copy constructor ";
  if ( Con3TestCopy != C2T1ER ) {
	if (verbose) std::cout << "unsuccessful."<<std::endl;
	numberFailedTests++;
  } else {
	if (verbose) std::cout << "successful."<<std::endl;
  }

  BDMatrix Con3TestCopyTrans( C2T1ER, Teuchos::TRANS );
  if(verbose) std::cout <<"constructor 3 -- copy constructor (transposed) ";
  if ( Con3TestCopyTrans(0, 2) != C2T1ER(2, 0) ) {
	if (verbose) std::cout << "unsuccessful."<<std::endl;
	numberFailedTests++;
  } else {
	if (verbose) std::cout << "successful."<<std::endl;
  }

  // constructor 4 (submatrix)

  BDMatrix Con4TestOrig(Teuchos::Copy, a, 4, 5, 5, 2, 1);
  BDMatrix C4TS;
  C4TS.shape( 3, 3, 2, 1 );
  C4TS(0, 0) = 12;  C4TS(0, 1) = 17;
  C4TS(1, 0) = 13;  C4TS(1, 1) = 18;  C4TS(1, 2) = 23;
  C4TS(2, 0) = 14;  C4TS(2, 1) = 19;  C4TS(2, 2) = 24;

  BDMatrix Con4TestCopy1(Teuchos::Copy, Con4TestOrig, 3, 3, 2);
  numberFailedTests += PrintTestResults("constructor 4 -- submatrix copy", Con4TestCopy1, C4TS, verbose);
  BDMatrix Con4TestCopy2(Teuchos::Copy, Con4TestOrig, 5, 5, 0);
  numberFailedTests += PrintTestResults("constructor 4 -- full matrix copy", Con4TestCopy2, Con4TestOrig, verbose);
  BDMatrix Con4TestView1(Teuchos::View, Con4TestOrig, 5, 5, 0);
  numberFailedTests += PrintTestResults("constructor 4 -- full matrix view", Con4TestView1, Con4TestOrig, verbose);
  BDMatrix Con4TestView2(Teuchos::View, Con4TestOrig, 3, 3, 2);
  numberFailedTests += PrintTestResults("constructor 4 -- submatrix view", Con4TestView2, C4TS, verbose);

  // Norm Tests

  BDMatrix AAA;
  AAA.shape( 5, 5, 2, 1 );
  AAA(0, 0) = 0;  AAA(0, 1) = 5;
  AAA(1, 0) = 1;  AAA(1, 1) = 6;  AAA(1, 2) = 11;
  AAA(2, 0) = 2;  AAA(2, 1) = 7;  AAA(2, 2) = 12;  AAA(2, 3) = 17;
		  AAA(3, 1) = 8;  AAA(3, 2) = 13;  AAA(3, 3) = 18;  AAA(3, 4) = 23;
				  AAA(4, 2) = 14;  AAA(4, 3) = 19;  AAA(4, 4) = 24;

  BDMatrix BBB;
  numberFailedTests += PrintTestResults("normOne of a 5x5", AAA.normOne(), 54.0, verbose);
  numberFailedTests += PrintTestResults("normInf of a 5x5", AAA.normInf(), 62.0, verbose);
  AAA = Teuchos::ScalarTraits<STYPE>::one();
  numberFailedTests += PrintTestResults("normFrobenius of a 5x5", AAA.normFrobenius(), 4.0, verbose);
  numberFailedTests += PrintTestResults("normOne of a 0x0", BBB.normOne(), 0.0, verbose);
  numberFailedTests += PrintTestResults("normInf of a 0x0", BBB.normInf(), 0.0, verbose);
  numberFailedTests += PrintTestResults("normFrobenius of a 0x0", BBB.normFrobenius(), 0.0, verbose);

  //  Set Method Tests.

  BDMatrix CCC( 5, 5, 2, 1 );
  //  Randomize the entries in CCC.
  testName = "random() -- enter random entries into matrix";
  returnCode = CCC.random();
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);
  //  Set the entries of CCC to 1.0.
  testName = "putScalar() -- set every entry of this matrix to 1.0";
  returnCode = CCC.putScalar(Teuchos::ScalarTraits<STYPE>::one());
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);
  //  Check assignment operator.
  BDMatrix CCC2( 5, 5, 2, 1 );
  CCC2.assign( CCC );
  if (verbose) std::cout <<  "assign() -- copy the values of an input matrix ";
  if ( CCC( 3, 4 ) == Teuchos::ScalarTraits<STYPE>::one() ) {
    if (verbose) std::cout<< "successful" <<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful" <<std::endl;
    numberFailedTests++;
  }
  //  Create a view into a submatrix of CCC
  BDMatrix CCCview( Teuchos::View, CCC, 3, 3 );
  BDMatrix CCCtest1( 3, 3, 2, 1 );
  CCCtest1 = CCCview;
  if (verbose) std::cout << "operator= -- small(empty) = large(view) ";
  if (CCCtest1.numRows()==3 && CCCtest1.values()==CCC.values()) {
    if (verbose) std::cout<< "successful" <<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful" <<std::endl;
    numberFailedTests++;
  }
  CCCtest1 = CCC;
  if (verbose) std::cout << "operator= -- small(view) = large(copy) ";
  if (CCCtest1.numRows()==5 && CCCtest1.values()!=CCC.values()) {
    if (verbose) std::cout<< "successful"<<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful"<<std::endl;
    numberFailedTests++;
  }
  BDMatrix CCCtest2( 3, 3, 2, 1 );
  CCCtest2 = Teuchos::ScalarTraits<STYPE>::one();
  CCCtest1 = CCCtest2;
  if (verbose) std::cout << "operator= -- large(copy) = small(copy) ";
  if (CCCtest1.numRows()==3 ) {
    if (verbose) std::cout<< "successful"<<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful"<<std::endl;
    numberFailedTests++;
  }
  CCCtest1 = CCCview;
  if (verbose) std::cout << "operator= -- large(copy) = small(view) ";
  if (CCCtest1.numRows()==3 && CCCtest1.stride()==4) {
    if(verbose) std::cout<<"successful" <<std::endl;
  } else {
    if (verbose) std::cout<<"unsuccessful"<<std::endl;
    numberFailedTests++;
  }
  BDMatrix CCCtest3( CCCview );
  CCCtest1 += CCCtest3;
  if (verbose) std::cout << "operator+= -- add two matrices of the same size, but different leading dimension ";
  if (CCCtest1(1,1)==2.0) {
    if(verbose) std::cout<<"successful" <<std::endl;
  } else {
    if (verbose) std::cout<<"unsuccessful"<<std::endl;
    numberFailedTests++;
  }
  if (verbose) std::cout << "operator+= -- add two matrices of different size (nothing should change) ";
  CCCtest1 += CCC;
  if (CCCtest1(1,1)==2.0) {
    if(verbose) std::cout<<"successful" <<std::endl;
  } else {
    if (verbose) std::cout<<"unsuccessful"<<std::endl;
    numberFailedTests++;
  }

  //  Scale Tests.

  BDMatrix ScalTest( 8, 8, 2, 3 );
  ScalTest = Teuchos::ScalarTraits<STYPE>::one();
  //  Scale the entries by 8, it should be 8.
  //  The matrix is lower triangular, by default, so check a lower triangular entry.
  if (verbose) std::cout << "operator*= -- scale matrix by some number ";
  ScalTest *= 8.0;
  if (ScalTest(5, 7) == 8.0) {
	if (verbose) std::cout<< "successful." <<std::endl;
  } else {
	if (verbose) std::cout<< "unsuccessful." <<std::endl;
	numberFailedTests++;
  }


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
