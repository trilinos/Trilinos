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
// Questions? Contact Robert Kirby (robert.c.kirby@ttu.edu) or
//                    Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file test_02.cpp
\brief  Unit tests for the FIAT::Default class
\author Created by R. Kirby
*/

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "FIAT_AffineMap.hpp"
#include "FIAT_DefaultLine.hpp"

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
    << "|                 Unit Test (FIAT::DefaultLine)                               |\n" \
    << "|                                                                             |\n" \
    << "|     1) Construction of [-1,1] line segment                                  |\n" \
    << "|     2) Failure modes                                                        |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Robert Kirby (robert.c.kirby@ttu.edu)                  |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n"\
    << "| TEST 1: Construction, exception testing                                     |\n"\
    << "===============================================================================\n";

  int errorFlag = 0;
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 3;

  // first error, get number of facets for illegal dimension
  try {
    DefaultLine<double> myLine;
    int numFacets = myLine.getNumFacets( -1 );
    // suppress "unused variable" warning:
    numFacets *= 1;
  }
  catch (std::invalid_argument err) {
    *outStream << "Expected error 1.) -------------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
  }
  // second error, number of facets exceeds spatial dimension
  try {
    DefaultLine<double> myLine;
    int numFacets = myLine.getNumFacets( 2 );
    // suppress "unused variable" warning
    numFacets *= 1;
  }
  catch (std::invalid_argument err) {
    *outStream << "Expected error 2.) -------------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
  }
  // third error: face tangents should barf
  try {
    DefaultLine<double> myLine;
    myLine.getFaceTangents();
  }
  catch (std::invalid_argument err) {
    *outStream << "Expected error 3.) -------------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
  }

  // Check if number of thrown exceptions matches the one we expect
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
    DefaultLine<double> myLine;
    
    // the first vertex should be -1.0
    TEST_FOR_EXCEPTION( fabs( myLine.getVerticesOfSubcomplex(0,0)->operator()(0,0)  + 1.0 ) > 1.e-10 ,
			std::runtime_error , 
			">>> FAILURE (FIAT::DefaultLine test_02) -- wrong value of vertex" );


    // spot check a member of the lattice
    TEST_FOR_EXCEPTION( fabs( myLine.makeLattice( 2 , 0 )->operator()(1,0)  ) > 1.e-10 ,
			std::runtime_error ,
			">>> FAILURE( FIAT::DefaultLine test_02) -- wrong value in lattice" );

    // check the tangent
    TEST_FOR_EXCEPTION( fabs( myLine.getScaledTangents( )->operator()(0,0) - 2.0 ) > 1.e-10 ,
			std::runtime_error ,
			">>> FAILURE( FIAT::DefaultLine test_02) -- wrong scaled tangent" );

    //check a normal
    TEST_FOR_EXCEPTION( fabs( myLine.getNormals( )->operator()(0,0) + 1.0 ) > 1.e-10 ,
			std::runtime_error ,
			">>> FAILURE( FIAT::DefaultLine test_02) -- wrong normal" );

    // area should be 2
    TEST_FOR_EXCEPTION( fabs( myLine.getMeasure( ) - 2.0 ) > 1.e-10 ,
			std::runtime_error ,
			">>> FAILURE( FIAT::DefaultLine test_02) -- wrong area" );

    // area of a vertex should be 1
    TEST_FOR_EXCEPTION( fabs( myLine.getMeasure( 0 , 1 ) - 1.0 ) > 1.e-10 ,
			std::runtime_error ,
			">>> FAILURE( FIAT::DefaultLine test_02) -- wrong vertex area" );

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
