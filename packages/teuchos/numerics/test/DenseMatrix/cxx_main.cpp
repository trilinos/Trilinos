// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"
#include "Teuchos_Version.hpp"

#define OTYPE int
#define STYPE std::complex<double>

template<typename TYPE>
int PrintTestResults(std::string, TYPE, TYPE, bool);

int ReturnCodeCheck(std::string, int, int, bool);

typedef double Real;
typedef Teuchos::SerialDenseVector<int, std::complex<Real> > DVector;
typedef Teuchos::SerialDenseMatrix<int, std::complex<Real> > DMatrix;
//typedef Teuchos::SerialDenseVector<OTYPE, STYPE> DVector;

int main(int argc, char* argv[])
{

  int i, j;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (verbose)
    std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  int numberFailedTests = 0;
  int returnCode = 0;
  std::string testName = "";



  if (verbose) std::cout<<std::endl<<"********** CHECKING TEUCHOS SERIAL DENSE MATRIX **********"<<std::endl<<std::endl;

  // default constructor test
  DMatrix DefConTest;
  if (verbose) std::cout <<"default constructor -- construct empty matrix ";
  if ( DefConTest.values()!=NULL || DefConTest.numCols()!=0 || DefConTest.numRows()!=0 ||DefConTest.stride()!=0 ||DefConTest.empty()!=true ) {
	if (verbose) std::cout << "unsuccessful."<<std::endl;
	numberFailedTests++;
  } else {
	if (verbose) std::cout << "successful."<<std::endl;
  }

  // constructor 1 (matrix w/ dimension but empty)

  DMatrix Con1Test( 3, 4 );
  if (verbose) std::cout <<"constructor 1 -- empty matrix with given dimensions ";
  if ( Con1Test.numRows()!=3 || Con1Test.numCols()!=4 || Con1Test( 1, 2 )!=0.0 ) {
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
  DMatrix Con2Test1ExpRes;
  Con2Test1ExpRes.shape(2, 3);
  Con2Test1ExpRes(0, 0) = 0;  Con2Test1ExpRes(0, 1) = 2; Con2Test1ExpRes(0, 2) = 4;
  Con2Test1ExpRes(1, 0) = 1;  Con2Test1ExpRes(1, 1) = 3; Con2Test1ExpRes(1, 2) = 5;

  DMatrix Con2Test1(Teuchos::Copy, a, 2, 2, 3);
  numberFailedTests += PrintTestResults("constructor 2 -- construct matrix from array subrange", Con2Test1, Con2Test1ExpRes, verbose);


  // constructor 3 (copy constructor)

  DMatrix Con3TestCopy( Con2Test1ExpRes );
  if(verbose) std::cout <<"constructor 3 -- copy constructor ";
  if ( Con3TestCopy != Con2Test1ExpRes ) {
	if (verbose) std::cout << "unsuccessful."<<std::endl;
	numberFailedTests++;
  } else {
	if (verbose) std::cout << "successful."<<std::endl;
  }

  DMatrix Con3TestCopyTrans( Con2Test1ExpRes, Teuchos::TRANS );
  if(verbose) std::cout <<"constructor 3 -- copy constructor (transposed) ";
  if ( Con3TestCopyTrans(2, 0) != Con2Test1ExpRes(0, 2) ) {
	if (verbose) std::cout << "unsuccessful."<<std::endl;
	numberFailedTests++;
  } else {
	if (verbose) std::cout << "successful."<<std::endl;
  }

  // constructor 4 (submatrix)

  DMatrix Con4TestOrig(Teuchos::Copy, a, 3, 3, 3);
  DMatrix Con4TestSubmatrix;
  Con4TestSubmatrix.shape(2, 2);
  Con4TestSubmatrix(0, 0) = 4; Con4TestSubmatrix(0, 1) = 7;
  Con4TestSubmatrix(1, 0) = 5; Con4TestSubmatrix(1, 1) = 8;
  DMatrix Con4TestCopy1(Teuchos::Copy, Con4TestOrig, 2, 2, 1, 1);
  numberFailedTests += PrintTestResults("constructor 4 -- submatrix copy", Con4TestCopy1, Con4TestSubmatrix, verbose);
  DMatrix Con4TestCopy2(Teuchos::Copy, Con4TestOrig, 3, 3, 0, 0);
  numberFailedTests += PrintTestResults("constructor 4 -- full matrix copy", Con4TestCopy2, Con4TestOrig, verbose);
  DMatrix Con4TestView1(Teuchos::View, Con4TestOrig, 2, 2, 1, 1);
  numberFailedTests += PrintTestResults("constructor 4 -- full matrix view", Con4TestView1, Con4TestSubmatrix, verbose);
  DMatrix Con4TestView2(Teuchos::View, Con4TestOrig, 3, 3, 0, 0);
  numberFailedTests += PrintTestResults("constructor 4 -- submatrix view", Con4TestView2, Con4TestOrig, verbose);

  // Norm Tests

  DMatrix AAA;
  AAA.shape(3, 3);
  AAA(0, 0) = 1; AAA(0, 1) = 2; AAA(0, 2) = 3;
  AAA(1, 0) = 4; AAA(1, 1) = 5; AAA(1, 2) = 6;
  AAA(2, 0) = 7; AAA(2, 1) = 8; AAA(2, 2) = 9;
  DMatrix BBB;
  numberFailedTests += PrintTestResults("normOne of a 3x3", AAA.normOne(), 18.0, verbose);
  numberFailedTests += PrintTestResults("normInf of a 3x3", AAA.normInf(), 24.0, verbose);
  AAA = Teuchos::ScalarTraits<STYPE>::one();
  numberFailedTests += PrintTestResults("normFrobenius of a 3x3", AAA.normFrobenius(), 3.0, verbose);
  numberFailedTests += PrintTestResults("normOne of a 0x0", BBB.normOne(), 0.0, verbose);
  numberFailedTests += PrintTestResults("normInf of a 0x0", BBB.normInf(), 0.0, verbose);
  numberFailedTests += PrintTestResults("normFrobenius of a 0x0", BBB.normFrobenius(), 0.0, verbose);

  // multiply() -- dimensions tests

  DMatrix DimTest0x0A, DimTest0x0B, DimTest2x0, DimTest1x2, DimTest2x1, DimTest2x2A, DimTest2x2B,
    DimTest3x3, DimTest0x2, DimTest0x0Result, DimTest1x1Result, DimTest2x0Result, DimTest1x2Result, DimTest2x1Result, DimTest2x2Result,
    DimTest2x3Result, DimTest0x2Result, DimTest3x3Result;

  DimTest0x2.shape(0, 2);
  DimTest2x0.shape(2, 0);
  DimTest1x2.shape(1, 2);
  DimTest2x1.shape(2, 1);
  DimTest2x2A.shape(2, 2);
  DimTest2x2B.shape(2, 2);
  DimTest3x3.shape(3, 3);
  DimTest0x2Result.shape(0, 2);
  DimTest1x1Result.shape(1, 1);
  DimTest2x0Result.shape(2, 0);
  DimTest1x2Result.shape(1, 2);
  DimTest2x1Result.shape(2, 1);
  DimTest2x2Result.shape(2, 2);
  DimTest2x3Result.shape(2, 3);
  DimTest3x3Result.shape(3, 3);

  returnCode = DimTest2x2Result.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, DimTest2x2A, DimTest2x2B, 1);
  testName = "multiply() -- dimensions -- compatible square matrices";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);
  returnCode = DimTest2x3Result.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, DimTest2x2A, DimTest3x3, 1);
  testName = "multiply() -- dimensions -- incompatible square matrices";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest1x1Result.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, DimTest1x2, DimTest2x1, 1);
  testName = "multiply() -- dimensions -- compatible nonsquare matrices";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);
  returnCode = DimTest2x2Result.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, DimTest2x1, DimTest1x2, 1);
  testName = "multiply() -- dimensions -- compatible nonsquare matrices";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);
  returnCode = DimTest2x1Result.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, DimTest2x1, DimTest2x1, 1);
  testName = "multiply() -- dimensions -- incompatible nonsquare matrices";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest1x2Result.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, DimTest1x2, DimTest1x2, 1);
  testName = "multiply() -- dimensions -- incompatible nonsquare matrices";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest2x2Result.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, DimTest2x0, DimTest2x2A, 1);
  testName = "multiply() -- dimensions -- first operand bad numCols";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest0x2Result.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, DimTest2x2A, DimTest0x2, 1);
  testName = "multiply() -- dimensions -- second operand bad numRows";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest2x2Result.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, DimTest2x2A, DimTest2x0, 1);
  testName = "multiply() -- dimensions -- second operand bad numCols";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);

  // multiply() -- multiplication results tests

  DMatrix MultTest2x2A, MultTest2x2B, MultTest3x3A, MultTest3x3B, MultTest2x2ATimes2x2B,
    MultTest3x3ATimes3x3B, MultTest2x2BTimes2x2A, MultTest3x3BTimes3x3A, MultTest2x2ATimes2x2BExpResult, MultTest2x2BTimes2x2AExpResult,
    MultTest3x3ATimes3x3BExpResult, MultTest3x3BTimes3x3AExpResult, MultTest2x3A, MultTest2x3B, MultTest3x2A, MultTest3x2B,
    MultTest2x3ATimes3x2B, MultTest3x2ATimes2x3B, MultTest2x3BTimes3x2A, MultTest3x2BTimes2x3A, MultTest2x3ATimes3x2BExpResult,
    MultTest3x2ATimes2x3BExpResult, MultTest2x3BTimes3x2AExpResult, MultTest3x2BTimes2x3AExpResult;

  MultTest2x2A.shape(2, 2);
  MultTest2x2B.shape(2, 2);
  MultTest3x3A.shape(3, 3);
  MultTest3x3B.shape(3, 3);
  MultTest2x2ATimes2x2B.shape(2, 2);
  MultTest2x2BTimes2x2A.shape(2, 2);
  MultTest3x3ATimes3x3B.shape(3, 3);
  MultTest3x3BTimes3x3A.shape(3, 3);
  MultTest2x2ATimes2x2BExpResult.shape(2, 2);
  MultTest2x2BTimes2x2AExpResult.shape(2, 2);
  MultTest3x3ATimes3x3BExpResult.shape(3, 3);
  MultTest3x3BTimes3x3AExpResult.shape(3, 3);
  MultTest2x3A.shape(2, 3);
  MultTest2x3B.shape(2, 3);
  MultTest3x2A.shape(3, 2);
  MultTest3x2B.shape(3, 2);
  MultTest2x3ATimes3x2B.shape(2, 2);
  MultTest3x2ATimes2x3B.shape(3, 3);
  MultTest2x3BTimes3x2A.shape(2, 2);
  MultTest3x2BTimes2x3A.shape(3, 3);
  MultTest2x3ATimes3x2BExpResult.shape(2, 2);
  MultTest3x2ATimes2x3BExpResult.shape(3, 3);
  MultTest2x3BTimes3x2AExpResult.shape(2, 2);
  MultTest3x2BTimes2x3AExpResult.shape(3, 3);

  for(i = 0; i < 2; i++)
    {
      for(j = 0; j < 2; j++)
	{
	  MultTest2x2A(i, j) = i + j;
	  MultTest2x2B(i, j) = (i * j) + 1;
	}
    }
  for(i = 0; i < 3; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  MultTest3x3A(i, j) = i + j;
	  MultTest3x3B(i, j) = (i * j) + 1;
	}
    }

  MultTest2x2ATimes2x2BExpResult(0, 0) = 1; MultTest2x2ATimes2x2BExpResult(0, 1) = 2;
  MultTest2x2ATimes2x2BExpResult(1, 0) = 3; MultTest2x2ATimes2x2BExpResult(1, 1) = 5;
  MultTest2x2BTimes2x2AExpResult(0, 0) = 1; MultTest2x2BTimes2x2AExpResult(0, 1) = 3;
  MultTest2x2BTimes2x2AExpResult(1, 0) = 2; MultTest2x2BTimes2x2AExpResult(1, 1) = 5;
  MultTest3x3ATimes3x3BExpResult(0, 0) = 3; MultTest3x3ATimes3x3BExpResult(0, 1) = 8; MultTest3x3ATimes3x3BExpResult(0, 2) = 13;
  MultTest3x3ATimes3x3BExpResult(1, 0) = 6; MultTest3x3ATimes3x3BExpResult(1, 1) = 14; MultTest3x3ATimes3x3BExpResult(1, 2) = 22;
  MultTest3x3ATimes3x3BExpResult(2, 0) = 9; MultTest3x3ATimes3x3BExpResult(2, 1) = 20; MultTest3x3ATimes3x3BExpResult(2, 2) = 31;
  MultTest3x3BTimes3x3AExpResult(0, 0) = 3; MultTest3x3BTimes3x3AExpResult(0, 1) = 6; MultTest3x3BTimes3x3AExpResult(0, 2) = 9;
  MultTest3x3BTimes3x3AExpResult(1, 0) = 8; MultTest3x3BTimes3x3AExpResult(1, 1) = 14; MultTest3x3BTimes3x3AExpResult(1, 2) = 20;
  MultTest3x3BTimes3x3AExpResult(2, 0) = 13; MultTest3x3BTimes3x3AExpResult(2, 1) = 22; MultTest3x3BTimes3x3AExpResult(2, 2) = 31;
  MultTest2x3A(0, 0) = 1; MultTest2x3A(0, 1) = 2; MultTest2x3A(0, 2) = 3;
  MultTest2x3A(1, 0) = 4; MultTest2x3A(1, 1) = 5; MultTest2x3A(1, 2) = 6;
  MultTest3x2A(0, 0) = 1; MultTest3x2A(0, 1) = 2;
  MultTest3x2A(1, 0) = 3; MultTest3x2A(1, 1) = 4;
  MultTest3x2A(2, 0) = 5; MultTest3x2A(2, 1) = 6;
  MultTest2x3B(0, 0) = 0; MultTest2x3B(0, 1) = 2; MultTest2x3B(0, 2) = 4;
  MultTest2x3B(1, 0) = 6; MultTest2x3B(1, 1) = 8; MultTest2x3B(1, 2) = 10;
  MultTest3x2B(0, 0) = 0; MultTest3x2B(0, 1) = 2;
  MultTest3x2B(1, 0) = 4; MultTest3x2B(1, 1) = 6;
  MultTest3x2B(2, 0) = 8; MultTest3x2B(2, 1) = 10;
  MultTest2x3ATimes3x2BExpResult(0, 0) = 32; MultTest2x3ATimes3x2BExpResult(0, 1) = 44;
  MultTest2x3ATimes3x2BExpResult(1, 0) = 68; MultTest2x3ATimes3x2BExpResult(1, 1) = 98;
  MultTest3x2ATimes2x3BExpResult(0, 0) = 12; MultTest3x2ATimes2x3BExpResult(0, 1) = 18; MultTest3x2ATimes2x3BExpResult(0, 2) = 24;
  MultTest3x2ATimes2x3BExpResult(1, 0) = 24; MultTest3x2ATimes2x3BExpResult(1, 1) = 38; MultTest3x2ATimes2x3BExpResult(1, 2) = 52;
  MultTest3x2ATimes2x3BExpResult(2, 0) = 36; MultTest3x2ATimes2x3BExpResult(2, 1) = 58; MultTest3x2ATimes2x3BExpResult(2, 2) = 80;
  MultTest2x3BTimes3x2AExpResult(0, 0) = 26; MultTest2x3BTimes3x2AExpResult(0, 1) = 32;
  MultTest2x3BTimes3x2AExpResult(1, 0) = 80; MultTest2x3BTimes3x2AExpResult(1, 1) = 104;
  MultTest3x2BTimes2x3AExpResult(0, 0) = 8; MultTest3x2BTimes2x3AExpResult(0, 1) = 10; MultTest3x2BTimes2x3AExpResult(0, 2) = 12;
  MultTest3x2BTimes2x3AExpResult(1, 0) = 28; MultTest3x2BTimes2x3AExpResult(1, 1) = 38; MultTest3x2BTimes2x3AExpResult(1, 2) = 48;
  MultTest3x2BTimes2x3AExpResult(2, 0) = 48; MultTest3x2BTimes2x3AExpResult(2, 1) = 66; MultTest3x2BTimes2x3AExpResult(2, 2) = 84;

  MultTest2x2ATimes2x2B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, MultTest2x2A, MultTest2x2B, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 2x2 * 2x2", MultTest2x2ATimes2x2B, MultTest2x2ATimes2x2BExpResult, verbose);
  MultTest2x2BTimes2x2A.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, MultTest2x2B, MultTest2x2A, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 2x2 * 2x2", MultTest2x2BTimes2x2A, MultTest2x2BTimes2x2AExpResult, verbose);
  MultTest3x3ATimes3x3B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, MultTest3x3A, MultTest3x3B, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 3x3 * 3x3", MultTest3x3ATimes3x3B, MultTest3x3ATimes3x3BExpResult, verbose);
  MultTest3x3BTimes3x3A.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, MultTest3x3B, MultTest3x3A, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 3x3 * 3x3", MultTest3x3BTimes3x3A, MultTest3x3BTimes3x3AExpResult, verbose);
  MultTest2x3ATimes3x2B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, MultTest2x3A, MultTest3x2B, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 2x3 * 3x2", MultTest2x3ATimes3x2B, MultTest2x3ATimes3x2BExpResult, verbose);
  MultTest2x3BTimes3x2A.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, MultTest2x3B, MultTest3x2A, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 2x3 * 3x2", MultTest2x3BTimes3x2A, MultTest2x3BTimes3x2AExpResult, verbose);
  MultTest3x2ATimes2x3B.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, MultTest3x2A, MultTest2x3B, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 3x2 * 2x3", MultTest3x2ATimes2x3B, MultTest3x2ATimes2x3BExpResult, verbose);
  MultTest3x2BTimes2x3A.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1, MultTest3x2B, MultTest2x3A, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 3x2 * 2x3", MultTest3x2BTimes2x3A, MultTest3x2BTimes2x3AExpResult, verbose);

  DMatrix MultTestHugeA, MultTestHugeB, MultTestHugeATimesHugeBExpResult,
    MultTestHugeATimesHugeB;

  const int hugeSize = 100;
  MultTestHugeA.shape(hugeSize, hugeSize);
  MultTestHugeB.shape(hugeSize, hugeSize);
  MultTestHugeATimesHugeBExpResult.shape(hugeSize, hugeSize);
  MultTestHugeATimesHugeB.shape(hugeSize, hugeSize);

  for(i = 0; i < hugeSize; i++)
  {
    for(j = 0; j < hugeSize; j++)
    {
      MultTestHugeA(i, j) = j;
      MultTestHugeB(i, j) = i;
      MultTestHugeATimesHugeBExpResult(i, j) = 328350;
    }
  }

  MultTestHugeATimesHugeB.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0,
    MultTestHugeA, MultTestHugeB, 1.0);
  numberFailedTests += PrintTestResults(
    "multiply() -- mult. results -- huge * huge",
    MultTestHugeATimesHugeB, MultTestHugeATimesHugeBExpResult, verbose);

  //
  //  Check scale methods.
  //
  DMatrix ScalTest( 8, 8 );
  ScalTest = Teuchos::ScalarTraits<STYPE>::one();
  //  Scale the entries by 8, it should be 8.
  if (verbose) std::cout << "scale() -- scale matrix by some number ";
  returnCode = ScalTest.scale( 8.0 );
  if (ScalTest(2, 3) == 8.0) {
	if (verbose) std::cout<< "successful." <<std::endl;
  } else {
	if (verbose) std::cout<< "unsuccessful." <<std::endl;
	numberFailedTests++;
  }
  //  Pointwise scale the entries by zero, they all should be zero.
  DMatrix ScalTest2( 8, 8 );
  if (verbose) std::cout << "scale() -- point-wise scale matrix ";
  ScalTest.scale( ScalTest2 );
  if (ScalTest.normOne() == 0.0) {
	if (verbose) std::cout<< "successful." <<std::endl;
  } else {
	if (verbose) std::cout<< "unsuccessful." <<std::endl;
	numberFailedTests++;
  }
  //
  //  Check set methods.
  //
  DMatrix CCC( 5, 5 );
  //  Randomize the entries in CCC.
  testName = "random() -- enter random entries into matrix";
  returnCode = CCC.random();
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);
  //  Set the entries of CCC to 1.0.
  testName = "putScalar() -- set every entry of this matrix to 1.0";
  returnCode = CCC.putScalar(Teuchos::ScalarTraits<STYPE>::one());
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);
  //  Check assignment operator.
  DMatrix CCC2( 5, 5 );
  CCC2.assign( CCC );
  if (verbose) std::cout <<  "assign() -- copy the values of an input matrix ";
  if ( CCC( 3, 4 ) == Teuchos::ScalarTraits<STYPE>::one() ) {
    if (verbose) std::cout<< "successful" <<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful" <<std::endl;
    numberFailedTests++;
  }
  //  Create a swap testing matrix
  DMatrix CCCswap( 2, 3 );
  CCCswap.random();
  DMatrix copyCCC(Teuchos::Copy, CCC);
  DMatrix copyCCCswap(Teuchos::Copy, CCCswap);
  if (verbose) std::cout <<  "swap() -- swap the values and attributes of two matrices -- ";
  CCCswap.swap(CCC);
  bool op_result = ( (CCCswap == copyCCC) && (CCC == copyCCCswap) );
  if (verbose)
    std::cout << (op_result ? "successful" : "failed" )<<std::endl;
  if( !op_result )
    numberFailedTests++;
  // Swap back using other matrix and allow downstream testing to proceed as if without swapping
  CCC.swap(CCCswap);

  //  Create a view into a submatrix of CCC
  DMatrix CCCview( Teuchos::View, CCC, 3, 3 );
  DMatrix CCCtest1( 2, 3 );
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
  DMatrix CCCtest2( 2, 2 );
  CCCtest2 = 3.0;
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

  DMatrix CCCtest3( CCCview );
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
  //
  //  Check overloaded operators.
  //
  op_result;
  MultTestHugeATimesHugeB.reshape(10, 10);
  op_result = (MultTestHugeATimesHugeB == MultTestHugeATimesHugeBExpResult);
  if (verbose) {
	std::cout << "operator== -- results -- small == huge "<< (op_result == false ? "successful" : "failed" )<<std::endl;
  }
  op_result = (MultTestHugeATimesHugeB != MultTestHugeATimesHugeBExpResult);
  if (verbose) {
	std::cout << "operator!= -- results -- small != huge "<< (op_result == true ? "successful" : "failed" )<<std::endl;
	std::cout << std::endl<< MultTestHugeATimesHugeB << std::endl;
  	//These won't work unless boundschecking is enabled.
  	//std::cout << MultTestHugeATimesHugeB(100, 1) << std::endl;
  	//std::cout << MultTestHugeATimesHugeB(1, 100) << std::endl;
  }


  if (verbose) std::cout<<std::endl<<"********** CHECKING TEUCHOS SERIAL DENSE VECTOR **********"<<std::endl<<std::endl;

  DVector DefConTestV;
  if (verbose) std::cout <<"default constructor -- construct empty std::vector ";
  if ( DefConTestV.values()!=NULL || DefConTestV.length()!=0 || DefConTestV.numRows()!=0 ||DefConTestV.stride()!=0 ) {
        if (verbose) std::cout << "unsuccessful."<<std::endl;
        numberFailedTests++;
  } else {
        if (verbose) std::cout << "successful."<<std::endl;
  }

  // constructor 1 (matrix w/ dimension but empty)

  DVector Con1TestV( 3 );
  if (verbose) std::cout <<"constructor 1 -- empty std::vector with given dimensions ";
  if ( Con1TestV.length()!=3 || Con1TestV.numCols()!=1 || Con1TestV( 1 )!=0.0 ) {
	if (verbose) std::cout << "unsuccessful."<<std::endl;
        numberFailedTests++;
  } else {
      	if (verbose) std::cout << "successful."<<std::endl;
  }
	
  // constructor 2 (from array) tests

  DVector Con2Test1V(Teuchos::Copy, a, 4);
  if (verbose) std::cout <<"constructor 2 -- construct std::vector from array subrange ";
  if ( Con2Test1V.numRows()!=4 || Con2Test1V.numCols()!=1 || Con2Test1V[ 2 ]!=2.0 ) {
        if (verbose) std::cout << "unsuccessful."<<std::endl;
        numberFailedTests++;
  } else {
        if (verbose) std::cout << "successful."<<std::endl;
  }

  // constructor 3 (copy constructor)

  DVector Con3TestCopyV( Con2Test1V );
  if(verbose) std::cout <<"constructor 3 -- copy constructor ";
  if ( Con3TestCopyV != Con2Test1V ) {
        if (verbose) std::cout << "unsuccessful."<<std::endl;
        numberFailedTests++;
  } else {
        if (verbose) std::cout << "successful."<<std::endl;
  }

  // non-member helper function (construct vector view of matrix column)

  OTYPE col = Teuchos::OrdinalTraits<OTYPE>::one();
  DVector ColViewTestV = Teuchos::getCol<OTYPE,STYPE>( Teuchos::View, AAA, col );
  if (verbose) std::cout <<"non-method helper function -- construct vector view of second column of matrix ";
  if ( ColViewTestV.normInf() != 1.0 || ColViewTestV.normOne() != 3.0 ) {
        if (verbose) std::cout << "unsuccessful."<<std::endl;
        numberFailedTests++;
  } else {
        if (verbose) std::cout << "successful."<<std::endl;
  }

  // checking norms

  numberFailedTests += PrintTestResults("normOne of a 3x1 std::vector", Con2Test1V.normOne(), 6.0, verbose);
  numberFailedTests += PrintTestResults("normInf of a 3x1 std::vector", Con2Test1V.normInf(), 3.0, verbose);
  Con2Test1V = Teuchos::ScalarTraits<STYPE>::one();
  numberFailedTests += PrintTestResults("normFrobenius of a 3x1 std::vector", Con2Test1V.normFrobenius(), 2.0, verbose);

  // check size/resize

  DVector SizeTestV1;
  SizeTestV1.size( 5 );
  if(verbose) std::cout <<"size() -- test ";
  if (SizeTestV1( 4 )!= 0.0) {
    if (verbose) std::cout << "unsuccessful."<<std::endl;
    numberFailedTests++;
  } else {
    if (verbose) std::cout << "successful."<<std::endl;
  }
  SizeTestV1 = 2.0*Teuchos::ScalarTraits<STYPE>::one();
  SizeTestV1.resize( 10 );
  if(verbose) std::cout <<"resize() -- test small --> large ";
  if (SizeTestV1[ 4 ]!= 2.0 || SizeTestV1[ 8 ]!=0.0 ) {
    if (verbose) std::cout << "unsuccessful."<<std::endl;
    numberFailedTests++;
  } else {
    if (verbose) std::cout << "successful."<<std::endl;
  }
  SizeTestV1.resize( 3 );
  if(verbose) std::cout <<"resize() -- test large --> small ";
  if (SizeTestV1( 2 )!= 2.0) {
    if (verbose) std::cout << "unsuccessful."<<std::endl;
    numberFailedTests++;
  } else {
    if (verbose) std::cout << "successful."<<std::endl;
  }

  DVector OpEqTestV1( 10 ); OpEqTestV1 = 3.0*Teuchos::ScalarTraits<STYPE>::one();
  DVector OpEqTestV2( Teuchos::View, OpEqTestV1.values(), 3 );
  DVector OpEqTestV3( 2 );
  OpEqTestV3 = OpEqTestV2;
  if (verbose) std::cout << "operator= -- small(empty) = large(view) ";
  if (OpEqTestV3.length()==3 && OpEqTestV3.values()==OpEqTestV2.values()) {
    if (verbose) std::cout<< "successful"<<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful"<<std::endl;
    numberFailedTests++;
  }
  OpEqTestV3 = OpEqTestV1;
  if (verbose) std::cout << "operator= -- small(view) = large(copy) ";
  if (OpEqTestV3.length()==10 && OpEqTestV3.values()!=OpEqTestV1.values()) {
    if (verbose) std::cout<< "successful"<<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful"<<std::endl;
    numberFailedTests++;
  }
  OpEqTestV3.size(5);
  OpEqTestV3 = OpEqTestV1;
  if (verbose) std::cout << "operator= -- small(copy) = large(copy) ";
  if (OpEqTestV3.length()==10 && OpEqTestV3.values()!=OpEqTestV1.values() && OpEqTestV3[ 9 ]==3.0) {
    if (verbose) std::cout<< "successful"<<std::endl;
  } else {
    if (verbose) std::cout<< "unsuccessful"<<std::endl;
    numberFailedTests++;
  }

  DVector OpSumTestV1( OpEqTestV2 );
  OpSumTestV1 += OpEqTestV2;
  if (verbose) std::cout << "operator+= -- add two vectors of the same size, but different leading dimension ";
  if (OpSumTestV1( 1 )==6.0) {
    if (verbose) std::cout<<"successful" <<std::endl;
  } else {
    if (verbose) std::cout<<"unsuccessful"<<std::endl;
    numberFailedTests++;
  }
  if (verbose) std::cout << "operator+= -- add two vectors of different size (nothing should change) ";
  OpSumTestV1 += OpEqTestV1;
  if (OpSumTestV1( 1 )==6.0) {
    if (verbose) std::cout<<"successful" <<std::endl;
  } else {
    if (verbose) std::cout<<"unsuccessful"<<std::endl;
    numberFailedTests++;
  }

  DVector OpCompTestV1( 5 );
  OpCompTestV1 = 2.0*Teuchos::ScalarTraits<STYPE>::one();
  if(verbose) std::cout <<"operator== -- test large == small ";
  if (OpCompTestV1 == SizeTestV1) {
    if (verbose) std::cout << "unsuccessful."<<std::endl;
    numberFailedTests++;
  } else {
    if (verbose) std::cout << "successful."<<std::endl;
  }
  if(verbose) std::cout <<"operator!= -- test large != small ";
  if (OpCompTestV1 != SizeTestV1) {
    if (verbose) std::cout << "successful."<<std::endl;
  } else {
    if (verbose) std::cout << "successful."<<std::endl;
    numberFailedTests++;
  }

  DVector ColSetTestV( AAA.numRows() );
  ColSetTestV.putScalar( 2.0 );
  bool ret = Teuchos::setCol<OTYPE,STYPE>( ColSetTestV, col, AAA );
  if (verbose) std::cout <<"non-method helper function -- set second column of matrix with vector ";
  if ( ColViewTestV.normInf() != 2.0 || ColViewTestV.normOne() != 6.0 || ret == false ) {
        if (verbose) std::cout << "unsuccessful."<<std::endl;
        numberFailedTests++;
  } else {
        if (verbose) std::cout << "successful."<<std::endl;
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
