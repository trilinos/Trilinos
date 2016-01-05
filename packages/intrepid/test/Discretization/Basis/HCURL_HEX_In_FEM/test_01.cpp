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
\brief  Unit tests for the Intrepid::C_HEX_In_FEM class.
\author Created by P. Bochev, D. Ridzal and K. Peterson.
*/
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HCURL_HEX_In_FEM.hpp"
#include "Intrepid_PointTools.hpp"
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
    << "|                 Unit Test (Basis_HCURL_HEX_In_FEM)                          |\n" \
    << "|                                                                             |\n" \
    << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n" \
    << "|     2) Basis values for VALUE and CURL operators                            |\n" \
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
  const int deg = 1;
  shards::CellTopology line(shards::getCellTopologyData< shards::Line<> >()); 
  FieldContainer<double> closedPts(PointTools::getLatticeSize(line,deg),1);
  FieldContainer<double> openPts(PointTools::getLatticeSize(line,deg+1,1),1);
  PointTools::getLattice<double,FieldContainer<double> >(closedPts,line,deg);
  PointTools::getLattice<double,FieldContainer<double> >(openPts,line,deg+1,1);

  Basis_HCURL_HEX_In_FEM<double, FieldContainer<double> > hexBasis(deg,closedPts,openPts);
  int errorFlag = 0;

  // Initialize throw counter for exception testing
  int nException     = 0;
  int throwCounter   = 0;

  // Define array containing the 8 vertices of the reference HEX, its center and 6 face centers
  FieldContainer<double> hexNodes(8, 3);
  hexNodes(0,0) = -1.0; hexNodes(0,1) = -1.0; hexNodes(0,2) = -1.0;
  hexNodes(1,0) = 1.0; hexNodes(1,1) = -1.0; hexNodes(1,2) = -1.0;
  hexNodes(2,0) = -1.0; hexNodes(2,1) = 1.0; hexNodes(2,2) = -1.0;
  hexNodes(3,0) = 1.0; hexNodes(3,1) = 1.0; hexNodes(3,2) = -1.0;
  hexNodes(4,0) = -1.0; hexNodes(4,1) = -1.0; hexNodes(4,2) = 1.0;
  hexNodes(5,0) = 1.0; hexNodes(5,1) = -1.0; hexNodes(5,2) = 1.0;
  hexNodes(6,0) = -1.0; hexNodes(6,1) = 1.0; hexNodes(6,2) = 1.0;
  hexNodes(7,0) = 1.0; hexNodes(7,1) = 1.0; hexNodes(7,2) = 1.0;


  
  // Generic array for the output values; needs to be properly resized depending on the operator type
  FieldContainer<double> vals;

  try{
    // exception #1: GRAD cannot be applied to HCURL functions 
    // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
    vals.resize(hexBasis.getCardinality(), hexNodes.dimension(0), 4 );
    INTREPID_TEST_COMMAND( hexBasis.getValues(vals, hexNodes, OPERATOR_GRAD), throwCounter, nException );

    // exception #2: DIV cannot be applied to HCURL functions
    // resize vals to rank-2 container with dimensions (num. points, num. basis functions)
    vals.resize(hexBasis.getCardinality(), hexNodes.dimension(0) );
    INTREPID_TEST_COMMAND( hexBasis.getValues(vals, hexNodes, OPERATOR_DIV), throwCounter, nException );
        
    // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    // exception #3
    INTREPID_TEST_COMMAND( hexBasis.getDofOrdinal(3,0,0), throwCounter, nException );
    // exception #4
    INTREPID_TEST_COMMAND( hexBasis.getDofOrdinal(1,1,1), throwCounter, nException );
    // exception #5
    INTREPID_TEST_COMMAND( hexBasis.getDofOrdinal(0,4,1), throwCounter, nException );
    // exception #6
    INTREPID_TEST_COMMAND( hexBasis.getDofTag(12), throwCounter, nException );
    // exception #7
    INTREPID_TEST_COMMAND( hexBasis.getDofTag(-1), throwCounter, nException );
    
#ifdef HAVE_INTREPID_DEBUG
    // Exceptions 8-15 test exception handling with incorrectly dimensioned input/output arrays
    // exception #8: input points array must be of rank-2
    FieldContainer<double> badPoints1(4, 5, 3);
    INTREPID_TEST_COMMAND( hexBasis.getValues(vals, badPoints1, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #9 dimension 1 in the input point array must equal space dimension of the cell
    FieldContainer<double> badPoints2(4, 2);
    INTREPID_TEST_COMMAND( hexBasis.getValues(vals, badPoints2, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #10 output values must be of rank-3 for OPERATOR_VALUE
    FieldContainer<double> badVals1(4, 3);
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals1, hexNodes, OPERATOR_VALUE), throwCounter, nException );
 
    // exception #11 output values must be of rank-3 for OPERATOR_CURL
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals1, hexNodes, OPERATOR_CURL), throwCounter, nException );
    
    // exception #12 incorrect 0th dimension of output array (must equal number of basis functions)
    FieldContainer<double> badVals2(hexBasis.getCardinality() + 1, hexNodes.dimension(0), 3);
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals2, hexNodes, OPERATOR_VALUE), throwCounter, nException ) ;
    
    // exception #13 incorrect 1st  dimension of output array (must equal number of points)
    FieldContainer<double> badVals3(hexBasis.getCardinality(), hexNodes.dimension(0) + 1, 3);
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals3, hexNodes, OPERATOR_VALUE), throwCounter, nException ) ;

    // exception #14: incorrect 2nd dimension of output array (must equal the space dimension)
    FieldContainer<double> badVals4(hexBasis.getCardinality(), hexNodes.dimension(0), 4);
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals4, hexNodes, OPERATOR_VALUE), throwCounter, nException ) ;
    
    // exception #15: incorrect 2nd dimension of output array (must equal the space dimension)
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals4, hexNodes, OPERATOR_CURL), throwCounter, nException ) ;
    
    // exception #16: D2 cannot be applied to HCURL functions 
    // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
    vals.resize(hexBasis.getCardinality(), 
                hexNodes.dimension(0),  
                Intrepid::getDkCardinality(OPERATOR_D2, hexBasis.getBaseCellTopology().getDimension()));
    INTREPID_TEST_COMMAND( hexBasis.getValues(vals, hexNodes, OPERATOR_D2), throwCounter, nException );
    
#endif
    
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };
  
  // Check if number of thrown exceptions matches the one we expect 
  // Note Teuchos throw number will not pick up exceptions 3-7 and therefore will not match.
  if (throwCounter != nException) {
    errorFlag++;
    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
  }
//#endif
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"\
    << "===============================================================================\n";
  
  try{
    std::vector<std::vector<int> > allTags = hexBasis.getAllDofTags();
    
    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    for (unsigned i = 0; i < allTags.size(); i++) {
      int bfOrd  = hexBasis.getDofOrdinal(allTags[i][0], allTags[i][1], allTags[i][2]);

//       for (unsigned j=0;j<4;j++) std::cout << allTags[i][j] << " "; std::cout << std::endl;
      
      std::vector<int> myTag = hexBasis.getDofTag(bfOrd);
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
    for( int bfOrd = 0; bfOrd < hexBasis.getCardinality(); bfOrd++) {
      std::vector<int> myTag  = hexBasis.getDofTag(bfOrd);
      int myBfOrd = hexBasis.getDofOrdinal(myTag[0], myTag[1], myTag[2]);
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
  
  // VALUE: Each row pair gives the 12x3 correct basis set values at an evaluation point: (P,F,D) layout
  double basisValues[] = {
    1,0,0, 1,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 1,0,0, 1,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,0,0, 1,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,0,0, 1,0,0,
    0,1,0, 0,0,0, 0,1,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,1,0, 0,0,0, 0,1,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,1,0, 0,0,0, 0,1,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,1,0, 0,0,0, 0,1,0,
    0,0,1, 0,0,0, 0,0,0, 0,0,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0,
    0,0,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0, 0,0,1, 0,0,0, 0,0,0,
    0,0,0, 0,0,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0, 0,0,1, 0,0,0,
    0,0,0, 0,0,0, 0,0,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0, 0,0,1
  };
  
  // CURL
  
  double basisCurls[] = {   
    0,-0.5,0.5, 0,-0.5,0.5, 0,0,0.5, 0,0,0.5, 0,-0.5,0, 0,-0.5,0, 0,0,0, 0,0,0,
    0,0,-0.5, 0,0,-0.5, 0,-0.5,-0.5, 0,-0.5,-0.5, 0,0,0, 0,0,0, 0,-0.5,0, 0,-0.5,0,
    0,0.5,0, 0,0.5,0, 0,0,0, 0,0,0, 0,0.5,0.5, 0,0.5,0.5, 0,0,0.5, 0,0,0.5,
    0,0,0, 0,0,0, 0,0.5,0, 0,0.5,0, 0,0,-0.5, 0,0,-0.5, 0,0.5,-0.5, 0,0.5,-0.5,
    // y-component basis functions
    // first y-component bf is (0,1/4(1-x)(1-z),0)
    // curl is (1/4(1-x),0,-1/4(1-z))
    0.5,0,-0.5, 0,0,-0.5, 0.5,0,-0.5, 0,0,-0.5, 0.5,0,0, 0,0,0, 0.5,0,0, 0,0,0,
    // second y-component bf is (0,1/4(1+x)(1-z),0)
    // curl is (1/4(1+x),0,1/4(1-z))
    0,0,0.5, 0.5,0,0.5, 0,0,0.5, 0.5,0,0.5, 0,0,0, 0.5,0,0, 0,0,0, 0.5,0,0,
    // third y-component bf is (0,1/4(1-x)(1+z),0)
    // curl is (-1/4(1-x),0,-1/4(1+z))
    -0.5,0,0, 0,0,0, -0.5,0,0, 0,0,0, -0.5,0,-0.5, 0,0,-0.5, -0.5,0,-0.5, 0,0,-0.5,
    // fourth y-component bf is (0,1/4(1+x)(1+z),0)
    // curl is (-1/4(1+x),0,1/4(1+z))
    0.0,0,0, -0.5,0,0, 0.0,0,0, -0.5,0,0, 0.0,0,0.5, -0.5,0,0.5, 0.0,0,0.5, -0.5,0,0.5,
    // first z-component bf is (0,0,1/4(1-x)(1-y))
    // curl is (-1/4(1-x),1/4(1-y),0)
    -0.5,0.5,0, 0,0.5,0, -0.5,0,0, 0,0,0, -0.5,0.5,0, 0,0.5,0, -0.5,0,0, 0,0,0,
    // second z-component bf is (0,0,1/4(1+x)(1-y))
    // curl is (-1/4(1+x),1/4(1-y),0)
    0.0,-0.5,0, -0.5,-0.5,0, 0,0,0, -0.5,0,0, 0,-0.5,0, -0.5,-0.5,0, 0,0,0, -0.5,0,0, 
    // third z-component bf is (0,0,1/4(1-x)(1+y))
    // curl is (1/4(1-x),1/4(1+y),0)
    0.5,0,0, 0,0,0, 0.5,0.5,0, 0,0.5,0, 0.5,0,0, 0,0,0, 0.5,0.5,0, 0,0.5,0, 
    // fourth z-component bf is (0,0,1/4(1+x)(1+y))
    // curl is (1/4(1+x),-1/4(1+y),0)
    0,0,0, 0.5,0,0, 0,-0.5,0, 0.5,-0.5,0, 0,0,0, 0.5,0,0, 0,-0.5,0, 0.5,-0.5,0 
  };

  
  try{
        
    // Dimensions for the output arrays:
    int numFields = hexBasis.getCardinality();
    int numPoints = hexNodes.dimension(0);
    int spaceDim  = hexBasis.getBaseCellTopology().getDimension();
    
    // Generic array for values and curls that will be properly sized before each call
    FieldContainer<double> vals;
    
    // Check VALUE of basis functions: resize vals to rank-3 container:
    vals.resize(numFields, numPoints, spaceDim);
    hexBasis.getValues(vals, hexNodes, OPERATOR_VALUE);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        for (int k = 0; k < spaceDim; k++) {
          
          // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
          int l = k + i * spaceDim * numPoints + j * spaceDim;
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
    
    // Check CURL of basis function: resize vals to rank-3 container
    vals.resize(numFields, numPoints, spaceDim);
    hexBasis.getValues(vals, hexNodes, OPERATOR_CURL);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        for (int k = 0; k < spaceDim; k++) {
          // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
           int l = k + i * spaceDim * numPoints + j * spaceDim;
           if (std::abs(vals(i,j,k) - basisCurls[l]) > INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        
             // Output the multi-index of the value where the error is: 
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
             *outStream << "}  computed curl component: " << vals(i,j,k) 
               << " but reference curl component: " << basisCurls[l] << "\n";
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
