#include <iostream>
#include <string>
#include "Teuchos_DenseMatrix.hpp"

using namespace std;
using namespace Teuchos;

template<typename PTRTYPE>
int PrintTestResults(string, PTRTYPE, PTRTYPE, bool);

int ReturnCodeCheck(string, int, int, bool);

#define USER_OTYPE int
#define USER_STYPE double

int main(int argc, char *argv[])
{

  typedef DenseMatrix<USER_OTYPE, USER_STYPE> DMatrix;

  int i, j, k;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int numberFailedTests = 0;
  int returnCode = 0;
  string testName = "";

  // constructor 1 (from array) tests

  double a[9];
  double b[4];
  for(i = 0; i < 9; i++)
    {
      if(i < 4) b[i] = i;
      a[i] = i;
    }

  DMatrix Con1Test1ExpRes;
  DMatrix Con1Test2ExpRes;
  Con1Test1ExpRes.shape(2, 2);
  Con1Test2ExpRes.shape(3, 3);
  Con1Test1ExpRes(0, 0) = 0;  Con1Test1ExpRes(0, 1) = 2; 
  Con1Test1ExpRes(1, 0) = 1;  Con1Test1ExpRes(1, 1) = 3; 
  
  DMatrix Con1Test1(Copy, a, 2, 2, 2);
  numberFailedTests += PrintTestResults("constructor 1 -- construct from array subrange", Con1Test1, Con1Test1ExpRes, verbose);

  DMatrix Con1Test2(Copy, b, 3, 3, 3);
  numberFailedTests += PrintTestResults("constructor 1 -- construct from array, exceeds bounds", Con1Test2, Con1Test2ExpRes, verbose);

  // constructor 3 (submatrix)

  DMatrix Con3TestOrig(Copy, a, 3, 3, 3);
  DMatrix Con3TestSubmatrix;
  Con3TestSubmatrix.shape(2, 2);
  Con3TestSubmatrix(0, 0) = 4; Con3TestSubmatrix(0, 1) = 7; 
  Con3TestSubmatrix(1, 0) = 5; Con3TestSubmatrix(1, 1) = 8;
  DMatrix Con3TestCopy1(Copy, Con3TestOrig, 2, 2, 1, 1);
  numberFailedTests += PrintTestResults("constructor 3 -- submatrix copy", Con3TestCopy1, Con3TestSubmatrix, verbose);
  DMatrix Con3TestCopy2(Copy, Con3TestOrig, 3, 3, 0, 0);
  numberFailedTests += PrintTestResults("constructor 3 -- full matrix copy", Con3TestCopy2, Con3TestOrig, verbose);
  DMatrix Con3TestView1(View, Con3TestOrig, 2, 2, 1, 1);
  numberFailedTests += PrintTestResults("constructor 3 -- full matrix view", Con3TestView1, Con3TestSubmatrix, verbose);
  DMatrix Con3TestView2(View, Con3TestOrig, 3, 3, 0, 0);
  numberFailedTests += PrintTestResults("constructor 3 -- submatrix view", Con3TestView2, Con3TestOrig, verbose);

  // oneNorm()

  DMatrix AAA;
  AAA.shape(3, 3);
  AAA(0, 0) = 1; AAA(0, 1) = 2; AAA(0, 2) = 3;
  AAA(1, 0) = 4; AAA(1, 1) = 5; AAA(1, 2) = 6;
  AAA(2, 0) = 7; AAA(2, 1) = 8; AAA(2, 2) = 9;
  DMatrix BBB;
  numberFailedTests += PrintTestResults("oneNorm of a 3x3", AAA.oneNorm(), 18.0, verbose);
  numberFailedTests += PrintTestResults("oneNorm of a 0x0", BBB.oneNorm(), 0.0, verbose);

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

  returnCode = DimTest2x2Result.multiply('N', 'N', 1, DimTest2x2A, DimTest2x2B, 1);
  testName = "multiply() -- dimensions -- compatible square matrices";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);
  returnCode = DimTest2x3Result.multiply('N', 'N', 1, DimTest2x2A, DimTest3x3, 1);
  testName = "multiply() -- dimensions -- incompatible square matrices";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest1x1Result.multiply('N', 'N', 1, DimTest1x2, DimTest2x1, 1);
  testName = "multiply() -- dimensions -- compatible nonsquare matrices";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);
  returnCode = DimTest2x2Result.multiply('N', 'N', 1, DimTest2x1, DimTest1x2, 1);
  testName = "multiply() -- dimensions -- compatible nonsquare matrices";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 0, verbose);
  returnCode = DimTest2x1Result.multiply('N', 'N', 1, DimTest2x1, DimTest2x1, 1);
  testName = "multiply() -- dimensions -- incompatible nonsquare matrices";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest1x2Result.multiply('N', 'N', 1, DimTest1x2, DimTest1x2, 1);
  testName = "multiply() -- dimensions -- incompatible nonsquare matrices";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest0x2Result.multiply('N', 'N', 1, DimTest0x2, DimTest2x2A, 1);
  testName = "multiply() -- dimensions -- first operand bad numRows";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest2x2Result.multiply('N', 'N', 1, DimTest2x0, DimTest2x2A, 1);
  testName = "multiply() -- dimensions -- first operand bad numCols";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest0x2Result.multiply('N', 'N', 1, DimTest2x2A, DimTest0x2, 1);
  testName = "multiply() -- dimensions -- second operand bad numRows";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest2x2Result.multiply('N', 'N', 1, DimTest2x2A, DimTest2x0, 1);
  testName = "multiply() -- dimensions -- second operand bad numCols";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest0x0Result.multiply('N', 'N', 1, DimTest0x0A, DimTest0x0B, 1);
  testName = "multiply() -- dimensions -- both operands bad numRows & numCols";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest2x2Result.multiply('N', 'N', 1, DimTest2x0, DimTest0x2, 1);
  testName = "multiply() -- dimensions -- first operand bad numCols, second operand bad numRows";
  numberFailedTests += ReturnCodeCheck(testName, returnCode, 1, verbose);
  returnCode = DimTest0x0Result.multiply('N', 'N', 1, DimTest0x2, DimTest2x0, 1);
  testName = "multiply() -- dimensions -- first operand bad numRows, second operand bad numCols";
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

  MultTest2x2ATimes2x2B.multiply('N', 'N', 1, MultTest2x2A, MultTest2x2B, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 2x2 * 2x2", MultTest2x2ATimes2x2B, MultTest2x2ATimes2x2BExpResult, verbose);
  MultTest2x2BTimes2x2A.multiply('N', 'N', 1, MultTest2x2B, MultTest2x2A, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 2x2 * 2x2", MultTest2x2BTimes2x2A, MultTest2x2BTimes2x2AExpResult, verbose);
  MultTest3x3ATimes3x3B.multiply('N', 'N', 1, MultTest3x3A, MultTest3x3B, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 3x3 * 3x3", MultTest3x3ATimes3x3B, MultTest3x3ATimes3x3BExpResult, verbose);
  MultTest3x3BTimes3x3A.multiply('N', 'N', 1, MultTest3x3B, MultTest3x3A, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 3x3 * 3x3", MultTest3x3BTimes3x3A, MultTest3x3BTimes3x3AExpResult, verbose);
  MultTest2x3ATimes3x2B.multiply('N', 'N', 1, MultTest2x3A, MultTest3x2B, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 2x3 * 3x2", MultTest2x3ATimes3x2B, MultTest2x3ATimes3x2BExpResult, verbose);
  MultTest2x3BTimes3x2A.multiply('N', 'N', 1, MultTest2x3B, MultTest3x2A, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 2x3 * 3x2", MultTest2x3BTimes3x2A, MultTest2x3BTimes3x2AExpResult, verbose);
  MultTest3x2ATimes2x3B.multiply('N', 'N', 1, MultTest3x2A, MultTest2x3B, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 3x2 * 2x3", MultTest3x2ATimes2x3B, MultTest3x2ATimes2x3BExpResult, verbose);
  MultTest3x2BTimes2x3A.multiply('N', 'N', 1, MultTest3x2B, MultTest2x3A, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- 3x2 * 2x3", MultTest3x2BTimes2x3A, MultTest3x2BTimes2x3AExpResult, verbose);

  DMatrix MultTestHugeA, MultTestHugeB, MultTestHugeATimesHugeBExpResult, MultTestHugeATimesHugeB;
  MultTestHugeA.shape(1000, 1000);
  MultTestHugeB.shape(1000, 1000);
  MultTestHugeATimesHugeBExpResult.shape(1000, 1000);
  MultTestHugeATimesHugeB.shape(1000, 1000);

  for(i = 0; i < 1000; i++)
    {
      for(j = 0; j < 1000; j++)
	{
	  MultTestHugeA(i, j) = j;
	  MultTestHugeB(i, j) = i;
	  MultTestHugeATimesHugeBExpResult(i, j) = 332833500;
	}
    }

  MultTestHugeATimesHugeB.multiply('N', 'N', 1, MultTestHugeA, MultTestHugeB, 1);
  numberFailedTests += PrintTestResults("multiply() -- mult. results -- huge * huge", MultTestHugeATimesHugeB, MultTestHugeATimesHugeBExpResult, verbose);

  if(numberFailedTests > 0) cout << "Number of failed tests: " << numberFailedTests << endl;

  return 0;
}

template<typename PTRTYPE>
int PrintTestResults(string testName, PTRTYPE calculatedResult, PTRTYPE expectedResult, bool verbose)
{
  int result;
  if(calculatedResult == expectedResult)
    {
      if(verbose) cout << testName << " successful." << endl;
      result = 0;
    }
  else
    {
      if(verbose) cout << testName << " unsuccessful." << endl;
      result = 1;
    }
  return result;
}

int ReturnCodeCheck(string testName, int returnCode, int expectedResult, bool verbose)
{
  int result;
  if(expectedResult == 0)
    {
      if(returnCode == 0)
	{
	  if(verbose) cout << testName << " test successful." << endl;
	  result = 0;
	}
      else
	{
	  if(verbose) cout << testName << " test unsuccessful. Return code was " << returnCode << "." << endl;
	  result = 1;
	}
    }
  else
    {
      if(returnCode != 0)
	{
	  if(verbose) cout << testName << " test successful -- failed as expected." << endl;
	  result = 0;
	}
      else
	{
	  if(verbose) cout << testName << " test unsuccessful -- did not fail as expected. Return code was " << returnCode << "." << endl;
	  result = 1;
	}
    }
  return result;
}
