// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file test_01.cpp
\brief  Unit tests for the FIAT::AffineMap class
\author Created by R. Kirby
*/

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "FIAT_AffineMap.hpp"

using namespace std;
using namespace FIAT;
using Teuchos::SerialDenseMatrix;

int main(int argc, char *argv[]) {

  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                 Unit Test (FIAT::AffineMap)                                 |\n" \
    << "|                                                                             |\n" \
    << "|     1) Construction of affine mappings between simplices                    |\n" \
    << "|     2) Failure modes                                                        |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Robert Kirby (robert.c.kirby@ttu.edu)                  |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n"\
    << "| TEST 1: Map construction, exception testing                                 |\n"\
    << "===============================================================================\n";

  int errorFlag = 0;
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 1;


  try {
    // FIRST ERROR: incompatible numbers of vertices
    // two vertices, 0 and 1 in R^1
    double xs[] = { 0.0 , 1.0 }; 
    SerialDenseMatrix<int,double> X( Teuchos::Copy , xs , 1 , 2 , 1 );
    // three vertices in R^2
    double ys[] = { 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 1.0 };
    SerialDenseMatrix<int,double> Y( Teuchos::Copy , ys , 1 , 3 , 2 );

    AffineMap<double> myMap( Teuchos::rcp( &X , false ) , Teuchos::rcp( &Y , false ) );
  }
  catch (std::invalid_argument err) {
    *outStream << "Expected error 1.) -------------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
  }

  // Check if number of thrown exceptions matches the one we expect (1)
  if (TestForException_getThrowNumber() != endThrowNumber) {
    errorFlag++;
    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
  }

  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 2: correctness of mapping, 1d to 1d cells                               \n"\
    << "===============================================================================\n";

  try{
    double xs[] = {0.0,1.0};
    double ys[] = {-1.0,1.0};
    SerialDenseMatrix<int,double> X( Teuchos::Copy , xs , 1 , 2 , 1 );
    SerialDenseMatrix<int,double> Y( Teuchos::Copy , ys , 1 , 2 , 1 );
    AffineMap<double> myMap( Teuchos::rcp( &X , false ) , Teuchos::rcp( &Y , false ) );
    SerialDenseMatrix<int,double> &a = *(myMap.getMatrix());
    SerialDenseMatrix<int,double> &b = *(myMap.getVector());

    TEST_FOR_EXCEPTION( a.numRows() != 1 || a.numCols() != 1 
                        || b.numRows() != 1 || b.numCols() != 1 
                        || fabs( a(0,0) - 2.0 ) > 1.e-10 
                        || fabs( b(0,0) + 1.0 ) > 1.e-10 ,
                        std::runtime_error , 
                        ">>> ERROR: TEST 2 Failed" );
  }
  catch (std::runtime_error &err){
     *outStream << err.what() << "\n\n";
     errorFlag = -1000;
  };

  try{
    double xs[] = {0.0,0.0,1.0,0.0,0.0,1.0};
    double ys[] = {-1.0,-1.0,1.0,-1.0,-1.0,1.0};
    SerialDenseMatrix<int,double> X( Teuchos::Copy , xs , 1 , 3 , 2 );
    SerialDenseMatrix<int,double> Y( Teuchos::Copy , ys , 1 , 3 , 2 );
    AffineMap<double> myMap( Teuchos::rcp( &X , false ) , Teuchos::rcp( &Y , false ) );
    SerialDenseMatrix<int,double> &a = *(myMap.getMatrix());
    SerialDenseMatrix<int,double> &b = *(myMap.getVector());

    TEST_FOR_EXCEPTION( a.numRows() != 2 || a.numCols() != 2 
                        || b.numRows() != 2 || b.numCols() != 1
                        || fabs( a(0,0) - 2.0 ) > 1.e-10 
                        || fabs( a(0,1) ) > 1.e-10
                        || fabs( a(1,0) ) > 1.e-10
                        || fabs( a(1,1) - 2.0 ) > 1.e-10
                        || fabs( b(0,0) + 1.0 ) > 1.e-10 
                        || fabs( b(1,0) + 1.0 ) > 1.e-10 ,
                        std::runtime_error , 
                        ">>> ERROR: TEST 2 Failed" );
  }
  catch (std::runtime_error &err){
     *outStream << err.what() << "\n\n";
     errorFlag = -1000;
  };



  if (errorFlag != 0) {
    std::cout << "End Result: TEST FAILED\n";
  }
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}

