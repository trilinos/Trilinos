// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file test_01.cpp
\brief  Unit tests for the Intrepid::Basis_HGRAD_LINE_Cn_FEM class.
\author Created by P. Bochev, D. Ridzal, K. Peterson.
*/

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_LINE_Cn_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace std;
using namespace Intrepid;

#define INTREPID_TEST_COMMAND( S , throwCounter, nException )                                                              \
{                                                                                                                          \
  ++nException;                                                                                                            \
  try {                                                                                                                    \
    S ;                                                                                                                    \
  }                                                                                                                        \
  catch (std::logic_error err) {                                                                                           \
      ++throwCounter;                                                                                                      \
      *outStream << "Expected Error " << nException << " -------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                                                                    \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";           \
  };                                                                                                                       \
}


int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

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
    << "|               Unit Test (Basis_HGRAD_LINE_Cn_FEM)                           |\n" \
    << "|                                                                             |\n" \
    << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n" \
    << "|     2) Basis values for VALUE, GRAD, CURL, and Dk operators                 |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n"\
    << "| TEST 1: Basis creation, exception testing                                   |\n"\
    << "===============================================================================\n";

  
  // Define basis and error flag
  shards::CellTopology line(shards::getCellTopologyData< shards::Line<> >());   // create cell topology
  const int deg = 5;

  // get the points for the basis
  FieldContainer<double> pts(PointTools::getLatticeSize(line,deg),1);
  PointTools::getLattice<double,FieldContainer<double> >(pts,line,deg);

  Basis_HGRAD_LINE_Cn_FEM<double, FieldContainer<double> > lineBasis(deg, pts);
  int errorFlag = 0;

  // Initialize throw counter for exception testing
  int nException     = 0;
  int throwCounter   = 0;
  
  // Define array containing vertices of the reference Line and a few other points   
  int numIntervals = 100;
  FieldContainer<double> lineNodes(numIntervals+1, 1);
  for (int i=0; i<numIntervals+1; i++) {
    lineNodes(i,0) = -1.0+(2.0*(double)i)/(double)numIntervals;
  }

  // Generic array for the output values; needs to be properly resized depending on the operator type
  FieldContainer<double> vals;

  try{
    // exception #1: DIV cannot be applied to scalar functions
    // resize vals to rank-2 container with dimensions (num. points, num. basis functions)
    vals.resize(lineBasis.getCardinality(), lineNodes.dimension(0) );
    //INTREPID_TEST_COMMAND( lineBasis.getValues(vals, lineNodes, OPERATOR_DIV), throwCounter, nException );

    // Exceptions 1-5: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    // exception #1
    INTREPID_TEST_COMMAND( lineBasis.getDofOrdinal(2,0,0), throwCounter, nException );
    // exception #2
    INTREPID_TEST_COMMAND( lineBasis.getDofOrdinal(1,1,1), throwCounter, nException );
     // exception #3
    INTREPID_TEST_COMMAND( lineBasis.getDofOrdinal(1,0,7), throwCounter, nException );
    // exception #4
    INTREPID_TEST_COMMAND( lineBasis.getDofTag(6), throwCounter, nException );
    // exception #5
    INTREPID_TEST_COMMAND( lineBasis.getDofTag(-1), throwCounter, nException );
    // not an exception
    INTREPID_TEST_COMMAND( lineBasis.getDofTag(5), throwCounter, nException ); --nException;

#ifdef HAVE_INTREPID_DEBUG
    // Exceptions 6-14 test exception handling with incorrectly dimensioned input/output arrays
    // exception #6: input points array must be of rank-2
    FieldContainer<double> badPoints1(4, 5, 3);
    INTREPID_TEST_COMMAND( lineBasis.getValues(vals, badPoints1, OPERATOR_VALUE), throwCounter, nException );

    // exception #7: dimension 1 in the input point array must equal space dimension of the cell
    FieldContainer<double> badPoints2(4, 3);
    INTREPID_TEST_COMMAND( lineBasis.getValues(vals, badPoints2, OPERATOR_VALUE), throwCounter, nException );

    // exception #8: output values must be of rank-2 for OPERATOR_VALUE
    FieldContainer<double> badVals1(4, 3, 1);
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals1, lineNodes, OPERATOR_VALUE), throwCounter, nException );

    // exception #9: output values must be of rank-3 for OPERATOR_GRAD
    FieldContainer<double> badVals2(4, 3);
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals2, lineNodes, OPERATOR_GRAD), throwCounter, nException );

    // exception #10: output values must be of rank-2 for OPERATOR_D1
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals2, lineNodes, OPERATOR_D1), throwCounter, nException );

    // exception #11: incorrect 0th dimension of output array (must equal number of basis functions)
    FieldContainer<double> badVals3(lineBasis.getCardinality() + 1, lineNodes.dimension(0));
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals3, lineNodes, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #12: incorrect 1st dimension of output array (must equal number of points)
    FieldContainer<double> badVals4(lineBasis.getCardinality(), lineNodes.dimension(0) + 1);
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals4, lineNodes, OPERATOR_VALUE), throwCounter, nException );

    // exception #13: incorrect 2nd dimension of output array (must equal spatial dimension)
    FieldContainer<double> badVals5(lineBasis.getCardinality(), lineNodes.dimension(0), 2);
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals5, lineNodes, OPERATOR_GRAD), throwCounter, nException );

    
    // not an exception
    FieldContainer<double> goodVals2(lineBasis.getCardinality(), lineNodes.dimension(0));
    INTREPID_TEST_COMMAND( lineBasis.getValues(goodVals2, lineNodes, OPERATOR_VALUE), throwCounter, nException ); --nException;
#endif

  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  // Check if number of thrown exceptions matches the one we expect
  if (throwCounter != nException) {
    errorFlag++;
    *outStream << std::setw(70) << "FAILURE! Incorrect number of exceptions." << "\n";
  }

   *outStream \
     << "\n"
     << "===============================================================================\n" \
     << "| TEST 2: correctness of basis function values                                |\n" \
     << "===============================================================================\n";
   outStream -> precision(20);


   try {
     // Check Kronecker property for Lagrange polynomials.
     int maxorder = 4;

     for (int ordi=1; ordi <= maxorder; ordi++) {
       FieldContainer<double> pts(PointTools::getLatticeSize(line,ordi),1);
       PointTools::getLattice<double,FieldContainer<double> >(pts,line,ordi);
       Basis_HGRAD_LINE_Cn_FEM<double, FieldContainer<double> > lineBasis(ordi, pts);
       FieldContainer<double> vals(lineBasis.getCardinality(),pts.dimension(0));
       lineBasis.getValues(vals,pts,OPERATOR_VALUE);
       for (int i=0;i<lineBasis.getCardinality();i++) {
  	for (int j=0;j<pts.dimension(0);j++) {
  	  if ( i==j && std::abs( vals(i,j) - 1.0 ) > INTREPID_TOL ) {
  	    errorFlag++;
  	    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
 	    *outStream << " Basis function " << i << " does not have unit value at its node\n";
 	  }
 	  if ( i!=j && std::abs( vals(i,j) ) > INTREPID_TOL ) {
 	    errorFlag++;
 	    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
 	    *outStream << " Basis function " << i << " does not vanish at node " << j << "\n";
 	  }
 	}
       }
     }
   }
  // Catch unexpected errors
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
