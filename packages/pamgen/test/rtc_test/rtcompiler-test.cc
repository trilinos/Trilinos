// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "RTC_FunctionRTC.hh"
#include "RTC_RegistrarRTC.hh"
#include "RTC_TokenizerRTC.hh"

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <cassert>

using namespace std;
using namespace PG_RuntimeCompiler;

void realProgram(double* arrayOne, double* arrayTwo,
                 double& plainValOne, double& plainValTwo, long& plainValThree, char testChar)
{
  int i;
  for (i = 0; i < 9; i = i + 1) {
    if (i == 0) {
      arrayOne[i] = -3.4e10;
      arrayTwo[i] = -arrayTwo[i];
    }
    else if (i == 1) {
      arrayOne[i] = fabs(arrayOne[i-1]);
      arrayTwo[i] = fabs(arrayTwo[i-1]);
    }
    else if (i == 2) {
      arrayOne[i] = sin(arrayOne[i] + pow((double)2,2));
      arrayTwo[i] = log(fabs(-arrayTwo[i]));
    }
    else if (i == 3) {
      arrayOne[i] = !(arrayOne[i] == arrayTwo[i]);
      arrayTwo[i] = !(arrayOne[i]);
    }
    else if (i == 4) {
      arrayOne[i] = sqrt(arrayTwo[i]);
      arrayTwo[i] = sqrt(arrayOne[i]);
    }
    else if (i == 5) {
      char temp = 'B';
      arrayOne[i] = fabs( (double)(temp - testChar));
      arrayTwo[i] = fabs( (double)(testChar - temp));
    }
    else if (i == 6) {
      arrayOne[i] = pow(testChar, arrayOne[i]) * (fabs(10e1));
      arrayTwo[i] = pow((double)testChar, 0) * (fabs(10e1));
    }
    else if (i == 7) {
      arrayOne[i] = sqrt(sin(3.4 + .002) + cos(2.3 / 2));
      arrayTwo[i] = sqrt(sin(3.4 + arrayTwo[(int)fabs((double)i)])
                         + cos(2.3 / 2));
    }
    else {
      int sum = 1;
      for (int a = 0; a < 10; a = a + 2) {
        sum = a + 2;
      }
      arrayOne[i] = sum + arrayOne[i];
      arrayTwo[i] = arrayTwo[i] / sum;
    }
  }

  for (i = 0; i < 9; ++i) {
    plainValOne = plainValOne + arrayOne[i];
    plainValTwo = plainValTwo + arrayTwo[i];
  }

  plainValThree = 123456789012;
}

int main( int argc, char* argv[] )
{
  Tokenizer::test();
  Line::test();

  assert( 1 == 1 );
  assert( 0 != 1 );

  Registrar::setupStandardFunctions();

  Function function(9);

  function.addVar("double[]", "arrayOne"     );
  function.addVar("double[]", "arrayTwo"     );
  function.addVar("double"  , "plainValOne"  );
  function.addVar("double"  , "plainValTwo"  );
  function.addVar("long"    , "plainValThree");
  function.addVar("char"    , "testChar"     );
  function.addVar("double[]", "testCol1"     );
  function.addVar("double[]", "testCol2"     );
  function.addVar("double[]", "testCol3"     );

  double arrayOne[9] = {0.0};
  double arrayOneReal[9] = {0.0};

  function.arrayAddrFill(0, arrayOne, 9);

  double arrayTwo[9] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
  double arrayTwoReal[9] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};

  function.arrayAddrFill(1, arrayTwo, 9);

  double plainValOne     = 0;
  double plainValOneReal = 0;
  function.varAddrFill(2, &plainValOne);

  double plainValTwo     = 5;
  double plainValTwoReal = 5;
  function.varAddrFill(3, &plainValTwo);

  long plainValThree     = 0;
  long plainValThreeReal = 0;
  function.varAddrFill(4, &plainValThree);

  char testChar = 'A';
  function.varValueFill(5, testChar);

  double testCol1[3] = {0.0};

  function.arrayAddrFill(6, testCol1, 3);

  double testCol2[3] = {0.0};

  function.arrayAddrFill(7, testCol2, 3);

  double testCol3[3] = {0.0};

  function.arrayAddrFill(8, testCol3, 3);

  string program = "int i; \n\
    for (i = 0; i < 9; i = i + 1) { \n\
    if (i == 0) { \n\
    arrayOne[i] = -3.4e10;\n\
      arrayTwo[i] = -arrayTwo[i];\n\
    }\n\
    else if (i == 1) \\{\n\
      arrayOne[i] = fabs(arrayOne[i-1]);\n\
      arrayTwo[i] = fabs(arrayTwo[i-1]);\n\
    \\}\n\
    else if (i == 2) \\{\n\
      arrayOne[i] = sin(arrayOne[i] + pow(2,2));\n\
      arrayTwo[i] = log(fabs(-arrayTwo[i]));\n\
    \\}\n\
    else if (i == 3) {\n\
      arrayOne[i] = !(arrayOne[i] == arrayTwo[i]);\n\
      arrayTwo[i] = !(arrayOne[i]);\n\
    }\n\
    else if (i == 4) {\n\
      arrayOne[i] = sqrt(arrayTwo[i]);\n\
      arrayTwo[i] = sqrt(arrayOne[i]);\n\
    }\n\
    else if (i == 5) {\n\
      char temp = 'B';\n\
      arrayOne[i] = fabs(temp - testChar);\n\
      arrayTwo[i] = fabs(testChar - temp);\n\
    }\n\
    else if (i == 6) {\n\
      arrayOne[i] = testChar^(arrayOne[i]) * (fabs(10.0e1));\n\
      arrayTwo[i] = pow(testChar, 0) * (fabs(10.0e1));\n\
    }\n\
    else if (i == 7) {\n\
      arrayOne[i] = sqrt(sin(3.4 + .002) + cos(2.3 / 2));\n\
      arrayTwo[i] = sqrt(sin(3.4 + arrayTwo[fabs(i)]) + cos(2.3 / 2));\n\
    }\n\
    else {\n\
      int sum = 1;\n\
      for (int a = 0; a < 10; a = a + 2) {\n\
        sum = a + 2;\n\
      }\n\
      arrayOne[i] = sum + arrayOne[i];\n\
      arrayTwo[i] = arrayTwo[i]/sum;\n\
    } \n\
  }\n\
  for (i = 0; i < 5; i = i + 1) \\{ \n\
    double temp = rand(); \n\
  }\n\
  for (i = 0; i < 5; i = i + 1) { \n\
    double temp = drand(); \n\
  }\n\
  for (i = 0; i < 9; i = i + 1) {\n\
    plainValOne = plainValOne + arrayOne[i];\n\
    plainValTwo = plainValTwo + arrayTwo[i];\n\
  }\n\
  plainValThree = 123456789012; \n\
  /* tes *(&^&*(^ sdfkjs 32 kjs f  * \n\
   *    * /  *                     */ \n\
  printf(\"One:% Two:% Three:% \", 5-4, 2.0e0, 'c');\n\
  int numRead = 0; \n\
  int index = 0; \n\
  while (numRead != -1) { \n\
    char buffer[256]; \n\
    numRead = readline(\"table.dat\", buffer); \n\
    if (numRead != -1) { \n\
      scanf(buffer, testCol1[index], testCol2[index], testCol3[index]); \n\
    } \n\
    index = index + 1; \n\
  }\n\
";

  function.addBody(program);

  function.execute();
  cout << function.getErrors() << endl;
  assert(function.getErrors() == "");

  cout << "############################################################" << endl;
  cout << "First run done" << endl;
  cout << "############################################################" << endl;

  for (int i = 0; i < 9; i++) {
    arrayOne[i]     = 0.0;
    arrayTwo[i]     = 5.0;
  }

  function.arrayAddrFill(0, arrayOne, 9);

  function.arrayAddrFill(1, arrayTwo, 9);

  plainValOne     = 0;
  function.varAddrFill(2, &plainValOne);

  plainValTwo     = 5;
  function.varAddrFill(3, &plainValTwo);

  plainValThree = 0;
  function.varAddrFill(4, &plainValThree);

  testChar = 'A';
  function.varValueFill(5, testChar);

  for (int i = 0; i < 3; ++i) {
    assert(testCol1[i] == 1.0 + i/10.0);
    assert(testCol2[i] == 2.0 + i/10.0);
    assert(testCol3[i] == 3.0 + i/10.0);

    testCol1[i] = 0.0;
    testCol2[i] = 0.0;
    testCol3[i] = 0.0;
  }

  function.arrayAddrFill(6, testCol1, 3);
  function.arrayAddrFill(7, testCol2, 3);
  function.arrayAddrFill(8, testCol3, 3);

  function.execute();
  cout << function.getErrors() << endl;
  assert(function.getErrors() == "");

  cout << "############################################################" << endl;
  cout << "Second run done" << endl;
  cout << "############################################################" << endl;

  for (int i = 0; i < 3; ++i) {
    assert(testCol1[i] == 1.0 + i/10.0);
    assert(testCol2[i] == 2.0 + i/10.0);
    assert(testCol3[i] == 3.0 + i/10.0);
  }

  realProgram(arrayOneReal, arrayTwoReal, plainValOneReal,
              plainValTwoReal, plainValThreeReal, testChar);

  [[maybe_unused]] const double epsilon = 1.0e-15;

  for (int i = 0; i < 9; ++i) {
    assert ( fabs(arrayOne[i] - arrayOneReal[i]) < epsilon);
    assert ( fabs(arrayTwo[i] - arrayTwoReal[i]) < epsilon);
  }
  assert ( fabs(plainValOne - plainValOneReal) < epsilon );
  assert ( fabs(plainValTwo - plainValTwoReal) < epsilon );

  assert( fabs(plainValThree - plainValThreeReal) < epsilon);

  cout << "\n === rtcunit completed successfully ===\n" << endl;

  return 0;
}
