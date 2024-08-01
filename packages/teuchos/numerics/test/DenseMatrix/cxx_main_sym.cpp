// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_SerialSymDenseMatrix.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "Teuchos_Version.hpp"

#define OTYPE int
#define STYPE std::complex<double>

template<typename TYPE>
int PrintTestResults(std::string, TYPE, TYPE, bool);

int ReturnCodeCheck(std::string, int, int, bool);

typedef Teuchos::SerialSymDenseMatrix<OTYPE, STYPE> SDMatrix;
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

  if (verbose) std::cout<<std::endl<<"********** CHECKING TEUCHOS SERIAL SYMMETRIC DENSE MATRIX **********"<<std::endl<<std::endl;

  // default constructor test
  SDMatrix DefConTest;
  if (verbose) std::cout <<"default constructor -- construct empty matrix ";
  if ( DefConTest.values()!=NULL || DefConTest.numCols()!=0 || DefConTest.numRows()!=0 ||DefConTest.stride()!=0 ||DefConTest.empty()!=true ) {
	if (verbose) std::cout << "unsuccessful."<<std::endl;
	numberFailedTests++;
  } else {
	if (verbose) std::cout << "successful."<<std::endl;
  }

  // constructor 1 (matrix w/ dimension but empty)

  SDMatrix Con1Test( 4 );
  if (verbose) std::cout <<"constructor 1 -- empty matrix with given dimensions ";
  if ( Con1Test.numRows()!=4 || Con1Test.numCols()!=4 || Con1Test( 1, 2 )!=0.0 ) {
	if (verbose) std::cout << "unsuccessful."<<std::endl;
        numberFailedTests++;
  } else {
        if (verbose) std::cout << "successful."<<std::endl;
  }
	
  // constructor 2 (from array) tests

  STYPE a[9];
  for(i = 0; i < 9; i++)
    {
      a[i] = i;
    }
  SDMatrix Con2Test1ExpRes;
  Con2Test1ExpRes.shape(3);
  Con2Test1ExpRes(0, 0) = 0;
  Con2Test1ExpRes(1, 0) = 1;  Con2Test1ExpRes(1, 1) = 4;
  Con2Test1ExpRes(2, 0) = 2;  Con2Test1ExpRes(2, 1) = 5;  Con2Test1ExpRes(2, 2) = 8;

  // Create another lower triangular matrix with a view of 'a'.
  SDMatrix Con2Test1(Teuchos::View, false, a, 3, 3);
  numberFailedTests += PrintTestResults("constructor 2 -- construct matrix from array subrange", Con2Test1, Con2Test1ExpRes, verbose);


  // constructor 3 (copy constructor)

  SDMatrix Con3TestCopy( Con2Test1ExpRes );
  if(verbose) std::cout <<"constructor 3 -- copy constructor ";
  if ( Con3TestCopy != Con2Test1ExpRes ) {
	if (verbose) std::cout << "unsuccessful."<<std::endl;
	numberFailedTests++;
  } else {
	if (verbose) std::cout << "successful."<<std::endl;
  }

  SDMatrix Con3TestCopyTrans( Con2Test1ExpRes );
  Con3TestCopyTrans.setUpper();
  if(verbose) std::cout <<"constructor 3 -- copy constructor (upper active storage) ";
  if ( Con3TestCopyTrans(2, 0) != Con2Test1ExpRes(2, 0) ) {
	if (verbose) std::cout << "unsuccessful."<<std::endl;
	numberFailedTests++;
  } else {
	if (verbose) std::cout << "successful."<<std::endl;
  }

  // constructor 4 (submatrix)

  SDMatrix Con4TestOrig(Teuchos::Copy, false, a, 3, 3);
  SDMatrix Con4TestSubmatrix;
  Con4TestSubmatrix.shape( 2 );
  Con4TestSubmatrix(0, 0) = 4;
  Con4TestSubmatrix(1, 0) = 5; Con4TestSubmatrix(1, 1) = 8;
  SDMatrix Con4TestCopy1(Teuchos::Copy, Con4TestOrig, 2, 1);
  numberFailedTests += PrintTestResults("constructor 4 -- submatrix copy", Con4TestCopy1, Con4TestSubmatrix, verbose);
  SDMatrix Con4TestCopy2(Teuchos::Copy, Con4TestOrig, 3, 0);
  numberFailedTests += PrintTestResults("constructor 4 -- full matrix copy", Con4TestCopy2, Con4TestOrig, verbose);
  SDMatrix Con4TestView1(Teuchos::View, Con4TestOrig, 2, 1);
  numberFailedTests += PrintTestResults("constructor 4 -- full matrix view", Con4TestView1, Con4TestSubmatrix, verbose);
  SDMatrix Con4TestView2(Teuchos::View, Con4TestOrig, 3, 0);
  numberFailedTests += PrintTestResults("constructor 4 -- submatrix view", Con4TestView2, Con4TestOrig, verbose);

  // Norm Tests

  SDMatrix AAA;
  AAA.shape( 3 );
  AAA(0, 0) = 8;
  AAA(1, 0) = 1; AAA(1, 1) = 8;
  AAA(2, 0) = 2; AAA(2, 1) = 3; AAA(2, 2) = 8;
  SDMatrix BBB;
  numberFailedTests += PrintTestResults("normOne of a 3x3", AAA.normOne(), 13.0, verbose);
  numberFailedTests += PrintTestResults("normInf of a 3x3", AAA.normInf(), 13.0, verbose);
  AAA = Teuchos::ScalarTraits<STYPE>::one();
  numberFailedTests += PrintTestResults("normFrobenius of a 3x3", AAA.normFrobenius(), 3.0, verbose);
  numberFailedTests += PrintTestResults("normOne of a 0x0", BBB.normOne(), 0.0, verbose);
  numberFailedTests += PrintTestResults("normInf of a 0x0", BBB.normInf(), 0.0, verbose);
  numberFailedTests += PrintTestResults("normFrobenius of a 0x0", BBB.normFrobenius(), 0.0, verbose);

  // Multiplication Tests

  // Reset values of AAA.
  AAA(0, 0) = 8;
  AAA(1, 0) = 1; AAA(1, 1) = 8;
  AAA(2, 0) = 2; AAA(2, 1) = 3; AAA(2, 2) = 8;

  DMatrix My_Prod( 4, 3 ), My_GenMatrix( 4, 3 );
  My_GenMatrix = Teuchos::ScalarTraits<STYPE>::one();

  // Matrix multiplication ( My_Prod = 1.0*My_GenMatrix*My_Matrix )
  My_Prod.multiply( Teuchos::RIGHT_SIDE, 1.0, AAA, My_GenMatrix, 0.0 );
  numberFailedTests += PrintTestResults("multiply() -- general times symmetric matrix (storage = lower tri)", My_Prod.normOne(), 52.0, verbose);
  AAA.setUpper();
  AAA(2, 1) = 1.0;
  My_Prod.multiply( Teuchos::RIGHT_SIDE, 1.0, AAA, My_GenMatrix, 0.0 );
  numberFailedTests += PrintTestResults("multiply() -- general times symmetric matrix (storage = upper tri)", My_Prod.normOne(), 44.0, verbose);

  //  Set Method Tests.

  SDMatrix CCC( 5 );
  //  Randomize the entries in CCC.
  testName = "random() -- enter random entries into matrix";
  returnCode = CCC.random();
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);
  //  Set the entries of CCC to 1.0.
  testName = "putScalar() -- set every entry of this matrix to 1.0";
  returnCode = CCC.putScalar(Teuchos::ScalarTraits<STYPE>::one());
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);
  //  Check assignment operator.
  SDMatrix CCC2( 5 );
  CCC2.assign( CCC );
  if (verbose) std::cout <<  "assign() -- copy the values of an input matrix ";
  if ( CCC( 3, 4 ) == Teuchos::ScalarTraits<STYPE>::one() ) {
    if (verbose) std::cout<< "successful" <<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful" <<std::endl;
    numberFailedTests++;
  }
  //  Create a swap testing matrix
  SDMatrix CCC2swap( 3 );
  CCC2swap.random();
  SDMatrix copyCCC2(CCC2);
  SDMatrix copyCCC2swap(CCC2swap);
  if (verbose) std::cout <<  "swap() -- swap the values and attributes of two symmetric matrices -- ";
  CCC2swap.swap(CCC2);
  bool op_result = ( (CCC2swap == copyCCC2) && (CCC2 == copyCCC2swap) );
  if (verbose)
    std::cout << (op_result ? "successful" : "failed" )<<std::endl;
  if( !op_result )
    numberFailedTests++;
  // Swap back using other matrix and allow downstream testing to proceed as if without swapping
  CCC2.swap(CCC2swap);

  //  Create a view into a submatrix of CCC
  SDMatrix CCCview( Teuchos::View, CCC, 3 );
  SDMatrix CCCtest1( 2 );
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
  SDMatrix CCCtest2( 2 );
  CCCtest2 = Teuchos::ScalarTraits<STYPE>::one();
  CCCtest1 = CCCtest2;
  if (verbose) std::cout << "operator= -- large(copy) = small(copy) ";
  if (CCCtest1.numRows()==2 ) {
    if (verbose) std::cout<< "successful"<<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful"<<std::endl;
    numberFailedTests++;
  }
  CCCtest1 = CCCview;
  if (verbose) std::cout << "operator= -- large(copy) = small(view) ";
  if (CCCtest1.numRows()==3 && CCCtest1.stride()==5) {
    if(verbose) std::cout<<"successful" <<std::endl;
  } else {
    if (verbose) std::cout<<"unsuccessful"<<std::endl;
    numberFailedTests++;
  }

  SDMatrix CCCtest3( CCCview );
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

  SDMatrix ScalTest( 8 );
  ScalTest = Teuchos::ScalarTraits<STYPE>::one();
  //  Scale the entries by 8, it should be 8.
  //  The matrix is lower triangular, by default, so check a lower triangular entry.
  if (verbose) std::cout << "operator*= -- scale matrix by some number ";
  ScalTest *= 8.0;
  if (ScalTest(7, 1) == 8.0) {
	if (verbose) std::cout<< "successful." <<std::endl;
  } else {
	if (verbose) std::cout<< "unsuccessful." <<std::endl;
	numberFailedTests++;
  }


  //  Matrix Triple-Product Test
  STYPE alpha=0.5*Teuchos::ScalarTraits<STYPE>::one();
  DMatrix W(3,2);
  SDMatrix A1(2), A2(3);
  A1(0,0) = 1.0, A1(1,1) = 2.0;
  A2(0,0) = 1.0, A2(1,1) = 2.0, A2(2,2) = 3.00;
  W = Teuchos::ScalarTraits<STYPE>::one();

  SDMatrix C1upper(3), C1lower(3), C2upper(2), C2lower(2);
  C1upper.setUpper(); C2upper.setUpper();
  C1lower.setLower(); C2lower.setLower();

  // Test all combinations of triple products.

  // These should return a matrix with 1.5 in all entries
  STYPE C1result = 1.5*Teuchos::ScalarTraits<STYPE>::one();
  Teuchos::symMatTripleProduct<OTYPE,STYPE>( Teuchos::NO_TRANS, alpha, A1, W, C1upper );
  Teuchos::symMatTripleProduct<OTYPE,STYPE>( Teuchos::NO_TRANS, alpha, A1, W, C1lower );

  // These should return a matrix with 3 in all entries
  STYPE C2result = 3.0*Teuchos::ScalarTraits<STYPE>::one();
  Teuchos::symMatTripleProduct<OTYPE,STYPE>( Teuchos::TRANS, alpha, A2, W, C2upper );
  Teuchos::symMatTripleProduct<OTYPE,STYPE>( Teuchos::TRANS, alpha, A2, W, C2lower );

  if (verbose) std::cout << "triple product -- compute C = W'*A*W or C = W*A*W' ";
  if (C1upper(2,1)==C1result && C1lower(1,2)==C1result && C2upper(1,0)==C2result && C2lower(0,1)==C2result) {
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
