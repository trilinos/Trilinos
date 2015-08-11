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
\brief  Unit tests for the Intrepid::HGRAD_TET_C2_FEM class.
\author Created by P. Bochev, D. Ridzal, and K. Peterson.
*/
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_TET_C2_FEM.hpp"
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
    << "|                 Unit Test (Basis_HGRAD_TET_C2_FEM)                          |\n" \
    << "|                                                                             |\n" \
    << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n" \
    << "|     2) Basis values for VALUE, GRAD, and Dk operators                       |\n" \
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
  Basis_HGRAD_TET_C2_FEM<double, FieldContainer<double> > tetBasis;
  int errorFlag = 0;

  // Initialize throw counter for exception testing
  int nException     = 0;
  int throwCounter   = 0;

  // Define array containing the 10 nodes of Tetrahedron<10>: 4 vertices + 6 edge midpoints  
  FieldContainer<double> tetNodes(10, 3);
  tetNodes(0,0) =  0.0;  tetNodes(0,1) =  0.0;  tetNodes(0,2) =  0.0;  
  tetNodes(1,0) =  1.0;  tetNodes(1,1) =  0.0;  tetNodes(1,2) =  0.0;  
  tetNodes(2,0) =  0.0;  tetNodes(2,1) =  1.0;  tetNodes(2,2) =  0.0;
  tetNodes(3,0) =  0.0;  tetNodes(3,1) =  0.0;  tetNodes(3,2) =  1.0;  
  
  tetNodes(4,0) =  0.5;  tetNodes(4,1) =  0.0;  tetNodes(4,2) =  0.0;
  tetNodes(5,0) =  0.5;  tetNodes(5,1) =  0.5;  tetNodes(5,2) =  0.0;  
  tetNodes(6,0) =  0.0;  tetNodes(6,1) =  0.5;  tetNodes(6,2) =  0.0;  
  tetNodes(7,0) =  0.0;  tetNodes(7,1) =  0.0;  tetNodes(7,2) =  0.5;  
  tetNodes(8,0) =  0.5;  tetNodes(8,1) =  0.0;  tetNodes(8,2) =  0.5;  
  tetNodes(9,0) =  0.0;  tetNodes(9,1) =  0.5;  tetNodes(9,2) =  0.5;  


  // Generic array for the output values; needs to be properly resized depending on the operator type
  FieldContainer<double> vals;

  try{
    // exception #1: CURL cannot be applied to scalar functions
    // resize vals to rank-3 container with dimensions (num. points, num. basis functions, arbitrary)
    vals.resize(tetBasis.getCardinality(), tetNodes.dimension(0), 4 );
    INTREPID_TEST_COMMAND( tetBasis.getValues(vals, tetNodes, OPERATOR_CURL), throwCounter, nException );

    // exception #2: DIV cannot be applied to scalar functions
    // resize vals to rank-2 container with dimensions (num. points, num. basis functions)
    vals.resize(tetBasis.getCardinality(), tetNodes.dimension(0) );
    INTREPID_TEST_COMMAND( tetBasis.getValues(vals, tetNodes, OPERATOR_DIV), throwCounter, nException );
        
    // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    // exception #3
    INTREPID_TEST_COMMAND( tetBasis.getDofOrdinal(3,0,0), throwCounter, nException );
    // exception #4
    INTREPID_TEST_COMMAND( tetBasis.getDofOrdinal(1,1,1), throwCounter, nException );
    // exception #5
    INTREPID_TEST_COMMAND( tetBasis.getDofOrdinal(0,4,0), throwCounter, nException );
    // exception #6
    INTREPID_TEST_COMMAND( tetBasis.getDofTag(10), throwCounter, nException );
    // exception #7
    INTREPID_TEST_COMMAND( tetBasis.getDofTag(-1), throwCounter, nException );
    
#ifdef HAVE_INTREPID_DEBUG 
    // Exceptions 8-18 test exception handling with incorrectly dimensioned input/output arrays
    // exception #8: input points array must be of rank-2
    FieldContainer<double> badPoints1(4, 5, 3);
    INTREPID_TEST_COMMAND( tetBasis.getValues(vals, badPoints1, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #9 dimension 1 in the input point array must equal space dimension of the cell
    FieldContainer<double> badPoints2(4, tetBasis.getBaseCellTopology().getDimension() - 1);
    INTREPID_TEST_COMMAND( tetBasis.getValues(vals, badPoints2, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #10 output values must be of rank-2 for OPERATOR_VALUE
    FieldContainer<double> badVals1(4, 3, 1);
    INTREPID_TEST_COMMAND( tetBasis.getValues(badVals1, tetNodes, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #11 output values must be of rank-3 for OPERATOR_GRAD
    FieldContainer<double> badVals2(4, 3);
    INTREPID_TEST_COMMAND( tetBasis.getValues(badVals2, tetNodes, OPERATOR_GRAD), throwCounter, nException );
    
    // exception #12 output values must be of rank-3 for OPERATOR_D1
    INTREPID_TEST_COMMAND( tetBasis.getValues(badVals2, tetNodes, OPERATOR_D1), throwCounter, nException );
    
    // exception #13 output values must be of rank-3 for OPERATOR_D2
    INTREPID_TEST_COMMAND( tetBasis.getValues(badVals2, tetNodes, OPERATOR_D2), throwCounter, nException );
    
    // exception #14 incorrect 0th dimension of output array (must equal number of basis functions)
    FieldContainer<double> badVals3(tetBasis.getCardinality() + 1, tetNodes.dimension(0));
    INTREPID_TEST_COMMAND( tetBasis.getValues(badVals3, tetNodes, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #15 incorrect 1st dimension of output array (must equal number of points)
    FieldContainer<double> badVals4(tetBasis.getCardinality(), tetNodes.dimension(0) + 1);
    INTREPID_TEST_COMMAND( tetBasis.getValues(badVals4, tetNodes, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
    FieldContainer<double> badVals5(tetBasis.getCardinality(), tetNodes.dimension(0), tetBasis.getBaseCellTopology().getDimension() + 1);
    INTREPID_TEST_COMMAND( tetBasis.getValues(badVals5, tetNodes, OPERATOR_GRAD), throwCounter, nException );
    
    // exception #17: incorrect 2nd dimension of output array (must equal D2 cardinality in 2D)
    FieldContainer<double> badVals6(tetBasis.getCardinality(), tetNodes.dimension(0), 40);
    INTREPID_TEST_COMMAND( tetBasis.getValues(badVals6, tetNodes, OPERATOR_D1), throwCounter, nException );
    
    // exception #18: incorrect 2nd dimension of output array (must equal D3 cardinality in 2D)
    INTREPID_TEST_COMMAND( tetBasis.getValues(badVals6, tetNodes, OPERATOR_D2), throwCounter, nException );
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
    std::vector<std::vector<int> > allTags = tetBasis.getAllDofTags();
    
    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    for (unsigned i = 0; i < allTags.size(); i++) {
      int bfOrd  = tetBasis.getDofOrdinal(allTags[i][0], allTags[i][1], allTags[i][2]);
      
      std::vector<int> myTag = tetBasis.getDofTag(bfOrd);
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
    for( int bfOrd = 0; bfOrd < tetBasis.getCardinality(); bfOrd++) {
      std::vector<int> myTag  = tetBasis.getDofTag(bfOrd);
      int myBfOrd = tetBasis.getDofOrdinal(myTag[0], myTag[1], myTag[2]);
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
  
  // VALUE: in (F,P) format
  double basisValues[] = {
    1.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00000, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 1.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00000, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 1.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00000, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    1.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00000, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 1.00000 };
  
  // GRAD and D1: in (F,P,D) format
  double basisGrads[] = {
    -3.00000, -3.00000, -3.00000, 1.00000, 1.00000, 1.00000, 1.00000, \
    1.00000, 1.00000, 1.00000, 1.00000, 1.00000, -1.00000, -1.00000, \
    -1.00000, 1.00000, 1.00000, 1.00000, -1.00000, -1.00000, -1.00000, \
    -1.00000, -1.00000, -1.00000, 1.00000, 1.00000, 1.00000, 1.00000, \
    1.00000, 1.00000, -1.00000, 0, 0, 3.00000, 0, 0, -1.00000, 0, 0, \
    -1.00000, 0, 0, 1.00000, 0, 0, 1.00000, 0, 0, -1.00000, 0, 0, \
    -1.00000, 0, 0, 1.00000, 0, 0, -1.00000, 0, 0, 0, -1.00000, 0, 0, \
    -1.00000, 0, 0, 3.00000, 0, 0, -1.00000, 0, 0, -1.00000, 0, 0, \
    1.00000, 0, 0, 1.00000, 0, 0, -1.00000, 0, 0, -1.00000, 0, 0, \
    1.00000, 0, 0, 0, -1.00000, 0, 0, -1.00000, 0, 0, -1.00000, 0, 0, \
    3.00000, 0, 0, -1.00000, 0, 0, -1.00000, 0, 0, -1.00000, 0, 0, \
    1.00000, 0, 0, 1.00000, 0, 0, 1.00000, 4.00000, 0, 0, -4.00000, \
    -4.00000, -4.00000, 0, 0, 0, 0, 0, 0, 0, -2.00000, -2.00000, \
    -2.00000, -2.00000, -2.00000, 2.00000, 0, 0, 2.00000, 0, 0, -2.00000, \
    -2.00000, -2.00000, 0, 0, 0, 0, 0, 0, 0, 4.00000, 0, 4.00000, 0, 0, \
    0, 0, 0, 0, 2.00000, 0, 2.00000, 2.00000, 0, 2.00000, 0, 0, 0, 0, 0, \
    0, 2.00000, 0, 2.00000, 0, 0, 0, 4.00000, 0, 0, 0, 0, -4.00000, \
    -4.00000, -4.00000, 0, 0, 0, 0, 2.00000, 0, -2.00000, -2.00000, \
    -2.00000, -2.00000, 0, -2.00000, 0, 2.00000, 0, 0, 0, 0, -2.00000, \
    -2.00000, -2.00000, 0, 0, 4.00000, 0, 0, 0, 0, 0, 0, -4.00000, \
    -4.00000, -4.00000, 0, 0, 2.00000, 0, 0, 0, 0, 0, 2.00000, -2.00000, \
    -2.00000, 0, -2.00000, -2.00000, -2.00000, -2.00000, -2.00000, \
    -2.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 4.00000, 0, 0, 0, 0, \
    2.00000, 0, 0, 2.00000, 0, 0, 0, 2.00000, 0, 0, 2.00000, 0, 2.00000, \
    2.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.00000, 0, 4.00000, 0, 0, 0, \
    0, 0, 0, 2.00000, 0, 0, 2.00000, 0, 2.00000, 0, 0, 2.00000, 0, 0, \
    2.00000, 2.00000};
  
  // D2 values in (F,P, Dk) format
  double basisD2[]={
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 4.00000, \
    4.00000, 4.00000, 4.00000, 4.00000, 4.00000, 0, 0, 0, 0, 0, 4.00000, \
    0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, \
    4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, \
    0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, \
    0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, \
    4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, \
    0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, \
    0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 0, 0, 4.00000, \
    0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, \
    4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, \
    0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, \
    0, 0, 4.00000, -8.00000, -4.00000, -4.00000, 0, 0, 0, -8.00000, \
    -4.00000, -4.00000, 0, 0, 0, -8.00000, -4.00000, -4.00000, 0, 0, 0, \
    -8.00000, -4.00000, -4.00000, 0, 0, 0, -8.00000, -4.00000, -4.00000, \
    0, 0, 0, -8.00000, -4.00000, -4.00000, 0, 0, 0, -8.00000, -4.00000, \
    -4.00000, 0, 0, 0, -8.00000, -4.00000, -4.00000, 0, 0, 0, -8.00000, \
    -4.00000, -4.00000, 0, 0, 0, -8.00000, -4.00000, -4.00000, 0, 0, 0, \
    0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, \
    0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, \
    0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, \
    0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, -4.00000, 0, -8.00000, -4.00000, \
    0, 0, -4.00000, 0, -8.00000, -4.00000, 0, 0, -4.00000, 0, -8.00000, \
    -4.00000, 0, 0, -4.00000, 0, -8.00000, -4.00000, 0, 0, -4.00000, 0, \
    -8.00000, -4.00000, 0, 0, -4.00000, 0, -8.00000, -4.00000, 0, 0, \
    -4.00000, 0, -8.00000, -4.00000, 0, 0, -4.00000, 0, -8.00000, \
    -4.00000, 0, 0, -4.00000, 0, -8.00000, -4.00000, 0, 0, -4.00000, 0, \
    -8.00000, -4.00000, 0, 0, 0, -4.00000, 0, -4.00000, -8.00000, 0, 0, \
    -4.00000, 0, -4.00000, -8.00000, 0, 0, -4.00000, 0, -4.00000, \
    -8.00000, 0, 0, -4.00000, 0, -4.00000, -8.00000, 0, 0, -4.00000, 0, \
    -4.00000, -8.00000, 0, 0, -4.00000, 0, -4.00000, -8.00000, 0, 0, \
    -4.00000, 0, -4.00000, -8.00000, 0, 0, -4.00000, 0, -4.00000, \
    -8.00000, 0, 0, -4.00000, 0, -4.00000, -8.00000, 0, 0, -4.00000, 0, \
    -4.00000, -8.00000, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, \
    0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, \
    0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, \
    0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 0, 0, \
    4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, \
    0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, \
    0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, 0, 0, 0, 4.00000, 0, 0, \
    0, 0, 0, 4.00000, 0
  };
  
  
  try{
        
    // Dimensions for the output arrays:
    int numFields = tetBasis.getCardinality();
    int numPoints = tetNodes.dimension(0);
    int spaceDim  = tetBasis.getBaseCellTopology().getDimension();
    
    // Generic array for values, grads, curls, etc. that will be properly sized before each call
    FieldContainer<double> vals;
    
    // Check VALUE of basis functions: resize vals to rank-2 container:
    vals.resize(numFields, numPoints);
    tetBasis.getValues(vals, tetNodes, OPERATOR_VALUE);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
          int l =  i + j * numFields;
           if (std::abs(vals(i,j) - basisValues[l]) > INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";
             *outStream << "}  computed value: " << vals(i,j)
               << " but reference value: " << basisValues[l] << "\n";
         }
      }
    }

    // Check GRAD of basis function: resize vals to rank-3 container
    vals.resize(numFields, numPoints, spaceDim);
    tetBasis.getValues(vals, tetNodes, OPERATOR_GRAD);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        for (int k = 0; k < spaceDim; k++) {
 
          // basisGrads is (F,P,D), compute offset:
          int l = k + j * spaceDim + i * spaceDim * numPoints;
           if (std::abs(vals(i,j,k) - basisGrads[l]) > INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
             *outStream << "}  computed grad component: " << vals(i,j,k)
               << " but reference grad component: " << basisGrads[l] << "\n";
            }
         }
      }
    }

    // Check D1 of basis function (do not resize vals because it has the correct size: D1 = GRAD)
    tetBasis.getValues(vals, tetNodes, OPERATOR_D1);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        for (int k = 0; k < spaceDim; k++) {
          
          // basisGrads is (F,P,D), compute offset:
          int l = k + j * spaceDim + i * spaceDim * numPoints;
           if (std::abs(vals(i,j,k) - basisGrads[l]) > INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
             *outStream << "}  computed D1 component: " << vals(i,j,k)
               << " but reference D1 component: " << basisGrads[l] << "\n";
            }
         }
      }
    }

    // Check D2 of basis function
    int D2cardinality = Intrepid::getDkCardinality(OPERATOR_D2, spaceDim);
    vals.resize(numFields, numPoints, D2cardinality);    
    tetBasis.getValues(vals, tetNodes, OPERATOR_D2);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        for (int k = 0; k < D2cardinality; k++) {
          
          // basisD2 is (F,P,Dk), compute offset:
          int l = k + j * D2cardinality + i * D2cardinality * numPoints;
          if (std::abs(vals(i,j,k) - basisD2[l]) > INTREPID_TOL) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            
            // Output the multi-index of the value where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
            *outStream << "}  computed D2 component: " << vals(i,j,k)
              << " but reference D2 component: " << basisD2[l] << "\n";
          }
        }
      }
    }
    
    
    // Check all higher derivatives - must be zero. 
    for(EOperator op = OPERATOR_D3; op < OPERATOR_MAX; op++) {
      
      // The last dimension is the number of kth derivatives and needs to be resized for every Dk
      int DkCardin  = Intrepid::getDkCardinality(op, spaceDim);
      vals.resize(numFields, numPoints, DkCardin);    

      tetBasis.getValues(vals, tetNodes, op);
      for (int i = 0; i < vals.size(); i++) {
        if (std::abs(vals[i]) > INTREPID_TOL) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          
          // Get the multi-index of the value where the error is and the operator order
          std::vector<int> myIndex;
          vals.getMultiIndex(myIndex,i);
          int ord = Intrepid::getOperatorOrder(op);
          *outStream << " At multi-index { ";
          for(int j = 0; j < vals.rank(); j++) {
            *outStream << myIndex[j] << " ";
          }
          *outStream << "}  computed D"<< ord <<" component: " << vals[i] 
            << " but reference D" << ord << " component:  0 \n";
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
