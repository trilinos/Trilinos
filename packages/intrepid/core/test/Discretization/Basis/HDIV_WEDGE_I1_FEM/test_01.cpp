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

/** \file   test_01.cpp
    \brief  Unit tests for the Intrepid::HDIV_WEDGE_I1_FEM class.
    \author Created by P. Bochev, D. Ridzal, and K. Peterson.
*/
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HDIV_WEDGE_I1_FEM.hpp"
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
    << "|                 Unit Test (Basis_HDIV_WEDGE_I1_FEM)                         |\n" \
    << "|                                                                             |\n" \
    << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n" \
    << "|     2) Basis values for VALUE and DIV operators                             |\n" \
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
  Basis_HDIV_WEDGE_I1_FEM<double, FieldContainer<double> > wedgeBasis;
  int errorFlag = 0;

  // Initialize throw counter for exception testing
  int nException     = 0;
  int throwCounter   = 0;

  // Define array containing the 6 vertices of the reference WEDGE and 6 other points.
  FieldContainer<double> wedgeNodes(12, 3);
  wedgeNodes(0,0) =  0.0;  wedgeNodes(0,1) =  0.0;  wedgeNodes(0,2) = -1.0;
  wedgeNodes(1,0) =  1.0;  wedgeNodes(1,1) =  0.0;  wedgeNodes(1,2) = -1.0;
  wedgeNodes(2,0) =  0.0;  wedgeNodes(2,1) =  1.0;  wedgeNodes(2,2) = -1.0;
  wedgeNodes(3,0) =  0.0;  wedgeNodes(3,1) =  0.0;  wedgeNodes(3,2) =  1.0;
  wedgeNodes(4,0) =  1.0;  wedgeNodes(4,1) =  0.0;  wedgeNodes(4,2) =  1.0;
  wedgeNodes(5,0) =  0.0;  wedgeNodes(5,1) =  1.0;  wedgeNodes(5,2) =  1.0;

  wedgeNodes(6,0) =  0.25; wedgeNodes(6,1) =  0.5;  wedgeNodes(6,2) = -1.0;
  wedgeNodes(7,0) =  0.5;  wedgeNodes(7,1) =  0.25; wedgeNodes(7,2) =  0.0;
  wedgeNodes(8,0) =  0.25; wedgeNodes(8,1) =  0.25; wedgeNodes(8,2) =  1.0;
  wedgeNodes(9,0) =  0.25; wedgeNodes(9,1) =  0.0;  wedgeNodes(9,2) =  0.75;
  wedgeNodes(10,0)=  0.0;  wedgeNodes(10,1)=  0.5;  wedgeNodes(10,2)= -0.25;
  wedgeNodes(11,0)=  0.5;  wedgeNodes(11,1)=  0.5;  wedgeNodes(11,2)=  0.0;

  // Generic array for the output values; needs to be properly resized depending on the operator type
  FieldContainer<double> vals;

  try{
    // exception #1: GRAD cannot be applied to HDIV functions 
    // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
    vals.resize(wedgeBasis.getCardinality(), wedgeNodes.dimension(0), 3 );
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_GRAD), throwCounter, nException );

    // exception #2: CURL cannot be applied to HDIV functions
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_CURL), throwCounter, nException );
        
    // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    // exception #3
    INTREPID_TEST_COMMAND( wedgeBasis.getDofOrdinal(3,0,0), throwCounter, nException );
    // exception #4
    INTREPID_TEST_COMMAND( wedgeBasis.getDofOrdinal(1,1,1), throwCounter, nException );
    // exception #5
    INTREPID_TEST_COMMAND( wedgeBasis.getDofOrdinal(0,4,1), throwCounter, nException );
    // exception #6
    INTREPID_TEST_COMMAND( wedgeBasis.getDofTag(11), throwCounter, nException );
    // exception #7
    INTREPID_TEST_COMMAND( wedgeBasis.getDofTag(-1), throwCounter, nException );

#ifdef HAVE_INTREPID_DEBUG
    // Exceptions 8-15 test exception handling with incorrectly dimensioned input/output arrays
    // exception #8: input points array must be of rank-2
    FieldContainer<double> badPoints1(4, 5, 3);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(vals, badPoints1, OPERATOR_VALUE), throwCounter, nException );

    // exception #9 dimension 1 in the input point array must equal space dimension of the cell
    FieldContainer<double> badPoints2(4, 2);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(vals, badPoints2, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #10 output values must be of rank-3 for OPERATOR_VALUE
    FieldContainer<double> badVals1(4, 3);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals1, wedgeNodes, OPERATOR_VALUE), throwCounter, nException );
 
    // exception #11 output values must be of rank-2 for OPERATOR_DIV
    FieldContainer<double> badVals2(4, 3, 1);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals2, wedgeNodes, OPERATOR_DIV), throwCounter, nException );
    
    // exception #12 incorrect 0th dimension of output array (must equal number of basis functions)
    FieldContainer<double> badVals3(wedgeBasis.getCardinality() + 1, wedgeNodes.dimension(0), 3);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals3, wedgeNodes, OPERATOR_VALUE), throwCounter, nException );

    // exception #13 incorrect 0th dimension of output array (must equal number of basis functions)
    FieldContainer<double> badVals4(wedgeBasis.getCardinality() + 1, wedgeNodes.dimension(0));
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals4, wedgeNodes, OPERATOR_DIV), throwCounter, nException );

    // exception #14 incorrect 1st dimension of output array (must equal number of points)
    FieldContainer<double> badVals5(wedgeBasis.getCardinality(), wedgeNodes.dimension(0) + 1, 3);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals5, wedgeNodes, OPERATOR_VALUE), throwCounter, nException );

    // exception #15 incorrect 1st dimension of output array (must equal number of points)
    FieldContainer<double> badVals6(wedgeBasis.getCardinality(), wedgeNodes.dimension(0) + 1);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals6, wedgeNodes, OPERATOR_DIV), throwCounter, nException );

    // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
    FieldContainer<double> badVals7(wedgeBasis.getCardinality(), wedgeNodes.dimension(0), 4);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals7, wedgeNodes, OPERATOR_VALUE), throwCounter, nException );
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
    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
  }
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"\
    << "===============================================================================\n";
  
  try{
    std::vector<std::vector<int> > allTags = wedgeBasis.getAllDofTags();
    
    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    for (unsigned i = 0; i < allTags.size(); i++) {
      int bfOrd  = wedgeBasis.getDofOrdinal(allTags[i][0], allTags[i][1], allTags[i][2]);
      
      std::vector<int> myTag = wedgeBasis.getDofTag(bfOrd);
       if( !( (myTag[0] == allTags[i][0]) &&
              (myTag[1] == allTags[i][1]) &&
              (myTag[2] == allTags[i][2]) &&
              (myTag[3] == allTags[i][3]) ) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getDofOrdinal( {" 
          << allTags[i][0] << ", " 
          << allTags[i][1] << ", " 
          << allTags[i][2] << ", " 
          << allTags[i][3] << "}) = " << bfOrd <<" but \n";   
        *outStream << " getDofTag(" << bfOrd << ") = { "
          << myTag[0] << ", " 
          << myTag[1] << ", " 
          << myTag[2] << ", " 
          << myTag[3] << "}\n";        
      }
    }
    
    // Now do the same but loop over basis functions
    for( int bfOrd = 0; bfOrd < wedgeBasis.getCardinality(); bfOrd++) {
      std::vector<int> myTag  = wedgeBasis.getDofTag(bfOrd);
      int myBfOrd = wedgeBasis.getDofOrdinal(myTag[0], myTag[1], myTag[2]);
      if( bfOrd != myBfOrd) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getDofTag(" << bfOrd << ") = { "
          << myTag[0] << ", " 
          << myTag[1] << ", " 
          << myTag[2] << ", " 
          << myTag[3] << "} but getDofOrdinal({" 
          << myTag[0] << ", " 
          << myTag[1] << ", " 
          << myTag[2] << ", " 
          << myTag[3] << "} ) = " << myBfOrd << "\n";
      }
    }
  }
  catch (std::logic_error err){
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 3: correctness of basis function values                                |\n"\
    << "===============================================================================\n";
  
  outStream -> precision(20);
  
  // VALUE: Each row pair gives the 5x3 correct basis set values at an evaluation point
  double basisValues[] = {
    0, -0.500000, 0, 0, 0, 0, -0.500000, 0, 0, 0, 0, -2.00000, 0, 0, 0, \
    0.500000, -0.500000, 0, 0.500000, 0, 0, 0, 0, 0, 0, 0, -2.00000, 0, \
    0, 0, 0, 0, 0, 0, 0.500000, 0, -0.500000, 0.500000, 0, 0, 0, \
    -2.00000, 0, 0, 0, 0, -0.500000, 0, 0, 0, 0, -0.500000, 0, 0, 0, 0, \
    0, 0, 0, 2.00000, 0.500000, -0.500000, 0, 0.500000, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 2.00000, 0, 0, 0, 0, 0.500000, 0, -0.500000, 0.500000, 0, \
    0, 0, 0, 0, 0, 2.00000, 0.125000, -0.250000, 0, 0.125000, 0.250000, \
    0, -0.375000, 0.250000, 0, 0, 0, -2.00000, 0, 0, 0, 0.250000, \
    -0.375000, 0, 0.250000, 0.125000, 0, -0.250000, 0.125000, 0, 0, 0, \
    -1.00000, 0, 0, 1.00000, 0.125000, -0.375000, 0, 0.125000, 0.125000, \
    0, -0.375000, 0.125000, 0, 0, 0, 0, 0, 0, 2.00000, 0.125000, \
    -0.500000, 0, 0.125000, 0, 0, -0.375000, 0, 0, 0, 0, -0.250000, 0, 0, \
    1.75000, 0, -0.250000, 0, 0, 0.250000, 0, -0.500000, 0.250000, 0, 0, \
    0, -1.25000, 0, 0, 0.750000, 0.250000, -0.250000, 0, 0.250000, \
    0.250000, 0, -0.250000, 0.250000, 0, 0, 0, -1.00000, 0, 0, 1.00000};
  
  // DIV: each row pair gives the 5 correct values of the divergence of the 5 basis functions
  double basisDivs[] = {   
    // 6 vertices
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0,
    // 6 other points
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0
  };
  
  try{
        
    // Dimensions for the output arrays:
    int numFields = wedgeBasis.getCardinality();
    int numPoints = wedgeNodes.dimension(0);
    int spaceDim  = wedgeBasis.getBaseCellTopology().getDimension();
    
    // Generic array for values and curls that will be properly sized before each call
    FieldContainer<double> vals;
    
    // Check VALUE of basis functions: resize vals to rank-3 container:
    vals.resize(numFields, numPoints, spaceDim);
    wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_VALUE);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        for (int k = 0; k < spaceDim; k++) {
           int l = k + i * spaceDim + j * spaceDim * numFields;
           if (std::abs(vals(i,j,k) - basisValues[l]) > INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
             *outStream << "}  computed value: " << vals(i,j,k)
               << " but reference value: " << basisValues[l] << "\n";
            }
         }
      }
    }

    
    // Check DIV of basis function: resize vals to rank-2 container
    vals.resize(numFields, numPoints);
    wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_DIV);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
          int l =  i + j * numFields;
           if (std::abs(vals(i,j) - basisDivs[l]) > INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";
             *outStream << "}  computed divergence component: " << vals(i,j)
               << " but reference divergence component: " << basisDivs[l] << "\n";
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
