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
\brief  Unit tests for the Intrepid::HGRAD_QUAD_C2_FEM class.
\author Created by P. Bochev, D. Ridzal, and K. Peterson.
*/
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
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
    << "|                 Unit Test (Basis_HGRAD_QUAD_C2_FEM)                         |\n" \
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
  Basis_HGRAD_QUAD_C2_FEM<double, FieldContainer<double> > quadBasis;
  int errorFlag = 0;

 // Initialize throw counter for exception testing
  int nException     = 0;
  int throwCounter   = 0;

  // Array with the 4 vertices, 4 edge midpoints, center of the reference QUAD and a random point.  
  FieldContainer<double> quadNodes(10, 2);
  quadNodes(0,0) = -1.0;  quadNodes(0,1) = -1.0;
  quadNodes(1,0) =  1.0;  quadNodes(1,1) = -1.0;
  quadNodes(2,0) =  1.0;  quadNodes(2,1) =  1.0;
  quadNodes(3,0) = -1.0;  quadNodes(3,1) =  1.0;
  // edge midpoints
  quadNodes(4,0) =  0.0;  quadNodes(4,1) = -1.0;
  quadNodes(5,0) =  1.0;  quadNodes(5,1) =  0.0;
  quadNodes(6,0) =  0.0;  quadNodes(6,1) =  1.0;
  quadNodes(7,0) = -1.0;  quadNodes(7,1) =  0.0;
  // center & random point
  quadNodes(8,0) =  0.0;  quadNodes(8,1) =  0.0;
  quadNodes(9,0) =1./3.;  quadNodes(9,1) =-3./5.;
  
  
  // Generic array for the output values; needs to be properly resized depending on the operator type
  FieldContainer<double> vals;

  try{
    // exception #1: DIV cannot be applied to scalar functions
    // resize vals to rank-2 container with dimensions (num. points, num. basis functions)
    vals.resize(quadBasis.getCardinality(), quadNodes.dimension(0));
    INTREPID_TEST_COMMAND( quadBasis.getValues(vals, quadNodes, OPERATOR_DIV), throwCounter, nException );
        
    // Exceptions 2-6: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    // exception #2
    INTREPID_TEST_COMMAND( quadBasis.getDofOrdinal(3,0,0), throwCounter, nException );
    // exception #3
    INTREPID_TEST_COMMAND( quadBasis.getDofOrdinal(1,1,1), throwCounter, nException );
    // exception #4
    INTREPID_TEST_COMMAND( quadBasis.getDofOrdinal(0,4,0), throwCounter, nException );
    // exception #5
    INTREPID_TEST_COMMAND( quadBasis.getDofTag(10), throwCounter, nException );
    // exception #6
    INTREPID_TEST_COMMAND( quadBasis.getDofTag(-1), throwCounter, nException );
    
#ifdef HAVE_INTREPID_DEBUG
    // Exceptions 7- test exception handling with incorrectly dimensioned input/output arrays
    // exception #7: input points array must be of rank-2
    FieldContainer<double> badPoints1(4, 5, 3);
    INTREPID_TEST_COMMAND( quadBasis.getValues(vals, badPoints1, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #8 dimension 1 in the input point array must equal space dimension of the cell
    FieldContainer<double> badPoints2(4, quadBasis.getBaseCellTopology().getDimension() + 1);
    INTREPID_TEST_COMMAND( quadBasis.getValues(vals, badPoints2, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #9 output values must be of rank-2 for OPERATOR_VALUE
    FieldContainer<double> badVals1(4, 3, 1);
    INTREPID_TEST_COMMAND( quadBasis.getValues(badVals1, quadNodes, OPERATOR_VALUE), throwCounter, nException );

    // exception #10 output values must be of rank-3 for OPERATOR_GRAD
    FieldContainer<double> badVals2(4, 3);
    INTREPID_TEST_COMMAND( quadBasis.getValues(badVals2, quadNodes, OPERATOR_GRAD), throwCounter, nException );
    
    // exception #11 output values must be of rank-3 for OPERATOR_CURL
    INTREPID_TEST_COMMAND( quadBasis.getValues(badVals2, quadNodes, OPERATOR_CURL), throwCounter, nException );
    
    // exception #12 output values must be of rank-3 for OPERATOR_D2
    INTREPID_TEST_COMMAND( quadBasis.getValues(badVals2, quadNodes, OPERATOR_D2), throwCounter, nException );
    
    // exception #13 incorrect 0th dimension of output array (must equal number of basis functions)
    FieldContainer<double> badVals3(quadBasis.getCardinality() + 1, quadNodes.dimension(0));
    INTREPID_TEST_COMMAND( quadBasis.getValues(badVals3, quadNodes, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #14 incorrect 1st dimension of output array (must equal number of points in quadNodes)
    FieldContainer<double> badVals4(quadBasis.getCardinality(), quadNodes.dimension(0) + 1);
    INTREPID_TEST_COMMAND( quadBasis.getValues(badVals4, quadNodes, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #15: incorrect 2nd dimension of output array (must equal the space dimension)
    FieldContainer<double> badVals5(quadBasis.getCardinality(), quadNodes.dimension(0), quadBasis.getBaseCellTopology().getDimension() + 1);
    INTREPID_TEST_COMMAND( quadBasis.getValues(badVals5, quadNodes, OPERATOR_GRAD), throwCounter, nException );
    
    // exception #16: incorrect 2nd dimension of output array (must equal D2 cardinality in 2D)
    FieldContainer<double> badVals6(quadBasis.getCardinality(), quadNodes.dimension(0), 40);
    INTREPID_TEST_COMMAND( quadBasis.getValues(badVals6, quadNodes, OPERATOR_D2), throwCounter, nException );
    
    // exception #17: incorrect 2nd dimension of output array (must equal D3 cardinality in 2D)
    INTREPID_TEST_COMMAND( quadBasis.getValues(badVals6, quadNodes, OPERATOR_D3), throwCounter, nException );
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
    std::vector<std::vector<int> > allTags = quadBasis.getAllDofTags();
    
    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    for (unsigned i = 0; i < allTags.size(); i++) {
      int bfOrd  = quadBasis.getDofOrdinal(allTags[i][0], allTags[i][1], allTags[i][2]);
      
      std::vector<int> myTag = quadBasis.getDofTag(bfOrd);
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
    for( int bfOrd = 0; bfOrd < quadBasis.getCardinality(); bfOrd++) {
      std::vector<int> myTag  = quadBasis.getDofTag(bfOrd);
      int myBfOrd = quadBasis.getDofOrdinal(myTag[0], myTag[1], myTag[2]);
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
  
  // VALUE: Correct basis values in (F,P) format:
  double basisValues[] = {
    1.000000000000000, 0, 0, 0, 0, 0, 0, 0, 0, -0.05333333333333334, 0, \
    1.000000000000000, 0, 0, 0, 0, 0, 0, 0, 0.1066666666666667, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 0, 0, -0.02666666666666666, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 0, 0.01333333333333333, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 0.4266666666666667, 0, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0.1422222222222222, 0, 0, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, -0.1066666666666667, 0, 0, 0, 0, 0, 0, 0, \
    1.000000000000000, 0, -0.07111111111111112, 0, 0, 0, 0, 0, 0, 0, 0, \
    1.000000000000000, 0.5688888888888890};
  
  // GRAD and D1: Correct gradients and D1 in (F,P,D) format: each row contains 6x2 values of
  double basisGrads[] = {
    -1.500000000000000, -1.500000000000000, 0.5000000000000000, 0, 0, 0, \
    0, 0.5000000000000000, -0.5000000000000000, 0, 0, 0, 0, 0, 0, \
    -0.5000000000000000, 0, 0, -0.08000000000000002, 0.1222222222222222, \
    -0.5000000000000000, 0, 1.500000000000000, -1.500000000000000, 0, \
    0.5000000000000000, 0, 0, 0.5000000000000000, 0, 0, \
    -0.5000000000000000, 0, 0, 0, 0, 0, 0, 0.3999999999999999, \
    -0.2444444444444444, 0, 0, 0, -0.5000000000000000, 1.500000000000000, \
    1.500000000000000, -0.5000000000000000, 0, 0, 0, 0, \
    0.5000000000000000, 0.5000000000000000, 0, 0, 0, 0, 0, \
    -0.09999999999999998, -0.02222222222222221, 0, -0.5000000000000000, \
    0, 0, 0.5000000000000000, 0, -1.500000000000000, 1.500000000000000, \
    0, 0, 0, 0, -0.5000000000000000, 0, 0, 0.5000000000000000, 0, 0, \
    0.02000000000000000, 0.01111111111111112, 2.000000000000000, 0, \
    -2.000000000000000, 0, 0, 0, 0, 0, 0, -1.500000000000000, 0, 0, 0, \
    0.5000000000000000, 0, 0, 0, -0.5000000000000000, \
    -0.3199999999999999, -0.9777777777777779, 0, 0, 0, 2.000000000000000, \
    0, -2.000000000000000, 0, 0, 0, 0, 1.500000000000000, 0, 0, 0, \
    -0.5000000000000000, 0, 0.5000000000000000, 0, 0.5333333333333334, \
    0.2666666666666666, 0, 0, 0, 0, -2.000000000000000, 0, \
    2.000000000000000, 0, 0, -0.5000000000000000, 0, 0, 0, \
    1.500000000000000, 0, 0, 0, 0.5000000000000000, 0.07999999999999997, \
    -0.08888888888888888, 0, 2.000000000000000, 0, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0.5000000000000000, 0, 0, 0, \
    -1.500000000000000, 0, -0.5000000000000000, 0, -0.1066666666666667, \
    -0.1333333333333333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.000000000000000, \
    -2.000000000000000, 0, 0, -2.000000000000000, 2.000000000000000, 0, \
    0, 0, -0.4266666666666667, 1.066666666666667};
  
  // D2: Correct multiset of second order partials in (F,P,Dk) format. D2 cardinality = 3 for 2D 
  double basisD2[] = {
    1.000000000000000, 2.250000000000000, 1.000000000000000, \
    1.000000000000000, -0.7500000000000000, 0, 0, 0.2500000000000000, 0, \
    0, -0.7500000000000000, 1.000000000000000, 1.000000000000000, \
    0.7500000000000000, 0, 0, -0.2500000000000000, 0, 0, \
    -0.2500000000000000, 0, 0, 0.7500000000000000, 1.000000000000000, 0, \
    0.2500000000000000, 0, 0.4800000000000000, 0.1833333333333334, \
    -0.1111111111111111, 1.000000000000000, 0.7500000000000000, 0, \
    1.000000000000000, -2.250000000000000, 1.000000000000000, 0, \
    0.7500000000000000, 1.000000000000000, 0, -0.2500000000000000, 0, \
    1.000000000000000, -0.7500000000000000, 0, 0, -0.7500000000000000, \
    1.000000000000000, 0, 0.2500000000000000, 0, 0, 0.2500000000000000, \
    0, 0, -0.2500000000000000, 0, 0.4800000000000000, \
    -0.9166666666666666, 0.2222222222222222, 0, 0.2500000000000000, 0, 0, \
    -0.7500000000000000, 1.000000000000000, 1.000000000000000, \
    2.250000000000000, 1.000000000000000, 1.000000000000000, \
    -0.7500000000000000, 0, 0, -0.2500000000000000, 0, 0, \
    0.7500000000000000, 1.000000000000000, 1.000000000000000, \
    0.7500000000000000, 0, 0, -0.2500000000000000, 0, 0, \
    0.2500000000000000, 0, -0.1200000000000000, -0.08333333333333331, \
    0.2222222222222222, 0, 0.7500000000000000, 1.000000000000000, 0, \
    -0.2500000000000000, 0, 1.000000000000000, 0.7500000000000000, 0, \
    1.000000000000000, -2.250000000000000, 1.000000000000000, 0, \
    0.2500000000000000, 0, 0, 0.2500000000000000, 0, 1.000000000000000, \
    -0.7500000000000000, 0, 0, -0.7500000000000000, 1.000000000000000, 0, \
    -0.2500000000000000, 0, -0.1200000000000000, 0.01666666666666666, \
    -0.1111111111111111, -2.000000000000000, -3.000000000000000, 0, \
    -2.000000000000000, 3.000000000000000, 0, 0, -1.000000000000000, 0, \
    0, 1.000000000000000, 0, -2.000000000000000, 0, 1.000000000000000, 0, \
    1.000000000000000, 0, 0, 0, 1.000000000000000, 0, -1.000000000000000, \
    0, 0, 0, 1.000000000000000, -0.9600000000000000, 0.7333333333333332, \
    0.8888888888888890, 0, -1.000000000000000, 0, 0, 3.000000000000000, \
    -2.000000000000000, 0, -3.000000000000000, -2.000000000000000, 0, \
    1.000000000000000, 0, 0, 1.000000000000000, 0, 1.000000000000000, 0, \
    -2.000000000000000, 0, -1.000000000000000, 0, 1.000000000000000, 0, \
    0, 1.000000000000000, 0, 0, 0.6400000000000001, 1.000000000000000, \
    -0.4444444444444444, 0, -1.000000000000000, 0, 0, 1.000000000000000, \
    0, -2.000000000000000, -3.000000000000000, 0, -2.000000000000000, \
    3.000000000000000, 0, 0, 0, 1.000000000000000, 0, -1.000000000000000, \
    0, -2.000000000000000, 0, 1.000000000000000, 0, 1.000000000000000, 0, \
    0, 0, 1.000000000000000, 0.2400000000000000, 0.06666666666666665, \
    0.8888888888888890, 0, -3.000000000000000, -2.000000000000000, 0, \
    1.000000000000000, 0, 0, -1.000000000000000, 0, 0, 3.000000000000000, \
    -2.000000000000000, 0, -1.000000000000000, 0, 1.000000000000000, 0, \
    0, 0, 1.000000000000000, 0, 1.000000000000000, 0, -2.000000000000000, \
    1.000000000000000, 0, 0, 0.6400000000000001, -0.2000000000000001, \
    0.2222222222222222, 0, 4.000000000000000, 0, 0, -4.000000000000000, \
    0, 0, 4.000000000000000, 0, 0, -4.000000000000000, 0, 0, 0, \
    -2.000000000000000, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, -2.000000000000000, 0, 0, -2.000000000000000, 0, \
    -2.000000000000000, -1.280000000000000, -0.7999999999999998, \
    -1.777777777777778};
  
  //D3: Correct multiset of second order partials in (F,P,Dk) format. D3 cardinality = 4 for 2D
  double basisD3[] = {
    0, -1.500000000000000, -1.500000000000000, 0, 0, -1.500000000000000, \
    0.5000000000000000, 0, 0, 0.5000000000000000, 0.5000000000000000, 0, \
    0, 0.5000000000000000, -1.500000000000000, 0, 0, -1.500000000000000, \
    -0.5000000000000000, 0, 0, -0.5000000000000000, 0.5000000000000000, \
    0, 0, 0.5000000000000000, -0.5000000000000000, 0, 0, \
    -0.5000000000000000, -1.500000000000000, 0, 0, -0.5000000000000000, \
    -0.5000000000000000, 0, 0, -1.100000000000000, -0.1666666666666667, \
    0, 0, -1.500000000000000, -0.5000000000000000, 0, 0, \
    -1.500000000000000, 1.500000000000000, 0, 0, 0.5000000000000000, \
    1.500000000000000, 0, 0, 0.5000000000000000, -0.5000000000000000, 0, \
    0, -1.500000000000000, 0.5000000000000000, 0, 0, -0.5000000000000000, \
    1.500000000000000, 0, 0, 0.5000000000000000, 0.5000000000000000, 0, \
    0, -0.5000000000000000, -0.5000000000000000, 0, 0, \
    -0.5000000000000000, 0.5000000000000000, 0, 0, -1.100000000000000, \
    0.8333333333333333, 0, 0, -0.5000000000000000, -0.5000000000000000, \
    0, 0, -0.5000000000000000, 1.500000000000000, 0, 0, \
    1.500000000000000, 1.500000000000000, 0, 0, 1.500000000000000, \
    -0.5000000000000000, 0, 0, -0.5000000000000000, 0.5000000000000000, \
    0, 0, 0.5000000000000000, 1.500000000000000, 0, 0, 1.500000000000000, \
    0.5000000000000000, 0, 0, 0.5000000000000000, -0.5000000000000000, 0, \
    0, 0.5000000000000000, 0.5000000000000000, 0, 0, \
    -0.09999999999999998, 0.8333333333333333, 0, 0, -0.5000000000000000, \
    -1.500000000000000, 0, 0, -0.5000000000000000, 0.5000000000000000, 0, \
    0, 1.500000000000000, 0.5000000000000000, 0, 0, 1.500000000000000, \
    -1.500000000000000, 0, 0, -0.5000000000000000, -0.5000000000000000, \
    0, 0, 0.5000000000000000, 0.5000000000000000, 0, 0, \
    1.500000000000000, -0.5000000000000000, 0, 0, 0.5000000000000000, \
    -1.500000000000000, 0, 0, 0.5000000000000000, -0.5000000000000000, 0, \
    0, -0.09999999999999998, -0.1666666666666667, 0, 0, \
    3.000000000000000, 2.000000000000000, 0, 0, 3.000000000000000, \
    -2.000000000000000, 0, 0, -1.000000000000000, -2.000000000000000, 0, \
    0, -1.000000000000000, 2.000000000000000, 0, 0, 3.000000000000000, 0, \
    0, 0, 1.000000000000000, -2.000000000000000, 0, 0, \
    -1.000000000000000, 0, 0, 0, 1.000000000000000, 2.000000000000000, 0, \
    0, 1.000000000000000, 0, 0, 0, 2.200000000000000, \
    -0.6666666666666665, 0, 0, 2.000000000000000, 1.000000000000000, 0, \
    0, 2.000000000000000, -3.000000000000000, 0, 0, -2.000000000000000, \
    -3.000000000000000, 0, 0, -2.000000000000000, 1.000000000000000, 0, \
    0, 2.000000000000000, -1.000000000000000, 0, 0, 0, \
    -3.000000000000000, 0, 0, -2.000000000000000, -1.000000000000000, 0, \
    0, 0, 1.000000000000000, 0, 0, 0, -1.000000000000000, 0, 0, \
    1.200000000000000, -1.666666666666667, 0, 0, 1.000000000000000, \
    2.000000000000000, 0, 0, 1.000000000000000, -2.000000000000000, 0, 0, \
    -3.000000000000000, -2.000000000000000, 0, 0, -3.000000000000000, \
    2.000000000000000, 0, 0, 1.000000000000000, 0, 0, 0, \
    -1.000000000000000, -2.000000000000000, 0, 0, -3.000000000000000, 0, \
    0, 0, -1.000000000000000, 2.000000000000000, 0, 0, \
    -1.000000000000000, 0, 0, 0, 0.2000000000000000, -0.6666666666666665, \
    0, 0, 2.000000000000000, 3.000000000000000, 0, 0, 2.000000000000000, \
    -1.000000000000000, 0, 0, -2.000000000000000, -1.000000000000000, 0, \
    0, -2.000000000000000, 3.000000000000000, 0, 0, 2.000000000000000, \
    1.000000000000000, 0, 0, 0, -1.000000000000000, 0, 0, \
    -2.000000000000000, 1.000000000000000, 0, 0, 0, 3.000000000000000, 0, \
    0, 0, 1.000000000000000, 0, 0, 1.200000000000000, 0.3333333333333334, \
    0, 0, -4.000000000000000, -4.000000000000000, 0, 0, \
    -4.000000000000000, 4.000000000000000, 0, 0, 4.000000000000000, \
    4.000000000000000, 0, 0, 4.000000000000000, -4.000000000000000, 0, 0, \
    -4.000000000000000, 0, 0, 0, 0, 4.000000000000000, 0, 0, \
    4.000000000000000, 0, 0, 0, 0, -4.000000000000000, 0, 0, 0, 0, 0, 0, \
    -2.400000000000000, 1.333333333333333, 0};
  //D4: Correct multiset of second order partials in (F,P,Dk) format. D4 cardinality = 5 for 2D
  double basisD4[] = {
    0, 0, 1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    1.000000000000000, 0, 0, 0, 0, 1.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    -2.000000000000000, 0, 0, 0, 0, -2.000000000000000, 0, 0, 0, 0, \
    4.000000000000000, 0, 0, 0, 0, 4.000000000000000, 0, 0, 0, 0, \
    4.000000000000000, 0, 0, 0, 0, 4.000000000000000, 0, 0, 0, 0, \
    4.000000000000000, 0, 0, 0, 0, 4.000000000000000, 0, 0, 0, 0, \
    4.000000000000000, 0, 0, 0, 0, 4.000000000000000, 0, 0, 0, 0, \
    4.000000000000000, 0, 0, 0, 0, 4.000000000000000, 0, 0};
  
  try{
        
    // Dimensions for the output arrays:
    int numFields = quadBasis.getCardinality();
    int numPoints = quadNodes.dimension(0);
    int spaceDim  = quadBasis.getBaseCellTopology().getDimension();
    
    // Generic array for values, grads, curls, etc. that will be properly sized before each call
    FieldContainer<double> vals;
    
    // Check VALUE of basis functions: resize vals to rank-2 container:
    vals.resize(numFields, numPoints);
    quadBasis.getValues(vals, quadNodes, OPERATOR_VALUE);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        
        // Compute offset for (F,P) container
        int l =  j + i * numPoints;
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
    quadBasis.getValues(vals, quadNodes, OPERATOR_GRAD);
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
    quadBasis.getValues(vals, quadNodes, OPERATOR_D1);
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


    // Check CURL of basis function: resize vals just for illustration! 
    vals.resize(numFields, numPoints, spaceDim);
    quadBasis.getValues(vals, quadNodes, OPERATOR_CURL);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        // We will use "rotated" basisGrads to check CURL: get offsets to extract (u_y, -u_x)
        int curl_0 = 1 + j * spaceDim + i * spaceDim * numPoints;               // position of y-derivative
        int curl_1 = 0 + j * spaceDim + i * spaceDim * numPoints;               // position of x-derivative
        
        double curl_value_0 = basisGrads[curl_0];
        double curl_value_1 =-basisGrads[curl_1];
        if (std::abs(vals(i,j,0) - curl_value_0) > INTREPID_TOL) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          // Output the multi-index of the value where the error is:
          *outStream << " At multi-index { ";
          *outStream << i << " ";*outStream << j << " ";*outStream << 0 << " ";
          *outStream << "}  computed curl component: " << vals(i,j,0)
            << " but reference curl component: " << curl_value_0 << "\n";
        }
        if (std::abs(vals(i,j,1) - curl_value_1) > INTREPID_TOL) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          // Output the multi-index of the value where the error is:
          *outStream << " At multi-index { ";
          *outStream << i << " ";*outStream << j << " ";*outStream << 1 << " ";
          *outStream << "}  computed curl component: " << vals(i,j,1)
            << " but reference curl component: " << curl_value_1 << "\n";
        }
      }
    }
    
    // Check D2 of basis function
    int D2cardinality = Intrepid::getDkCardinality(OPERATOR_D2, spaceDim);
    vals.resize(numFields, numPoints, D2cardinality);    
    quadBasis.getValues(vals, quadNodes, OPERATOR_D2);
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

    
    // Check D3 of basis function
    int D3cardinality = Intrepid::getDkCardinality(OPERATOR_D3, spaceDim);
    vals.resize(numFields, numPoints, D3cardinality);    
    quadBasis.getValues(vals, quadNodes, OPERATOR_D3);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        for (int k = 0; k < D3cardinality; k++) {
          
          // basisD3 is (F,P,Dk), compute offset:
          int l = k + j * D3cardinality + i * D3cardinality * numPoints;
          if (std::abs(vals(i,j,k) - basisD3[l]) > INTREPID_TOL) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            
            // Output the multi-index of the value where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
            *outStream << "}  computed D3 component: " << vals(i,j,k)
              << " but reference D3 component: " << basisD2[l] << "\n";
          }
        }
      }
    }

    
    // Check D4 of basis function
    int D4cardinality = Intrepid::getDkCardinality(OPERATOR_D4, spaceDim);
    vals.resize(numFields, numPoints, D4cardinality);    
    quadBasis.getValues(vals, quadNodes, OPERATOR_D4);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        for (int k = 0; k < D4cardinality; k++) {
          
          // basisD4 is (F,P,Dk), compute offset:
          int l = k + j * D4cardinality + i * D4cardinality * numPoints;
          if (std::abs(vals(i,j,k) - basisD4[l]) > INTREPID_TOL) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            
            // Output the multi-index of the value where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
            *outStream << "}  computed D4 component: " << vals(i,j,k)
              << " but reference D4 component: " << basisD2[l] << "\n";
          }
        }
      }
    }
    

    // Check all higher derivatives - must be zero. 
    for(EOperator op = OPERATOR_D5; op < OPERATOR_MAX; op++) {
      
      // The last dimension is the number of kth derivatives and needs to be resized for every Dk
      int DkCardin  = Intrepid::getDkCardinality(op, spaceDim);
      vals.resize(numFields, numPoints, DkCardin);    

      quadBasis.getValues(vals, quadNodes, op);
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


  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 4: correctness of DoF locations                                        |\n"\
    << "===============================================================================\n";

  try{
    Teuchos::RCP<Basis<double, FieldContainer<double> > > basis =
      Teuchos::rcp(new Basis_HGRAD_QUAD_C2_FEM<double, FieldContainer<double> >);
    Teuchos::RCP<DofCoordsInterface<FieldContainer<double> > > coord_iface =
      Teuchos::rcp_dynamic_cast<DofCoordsInterface<FieldContainer<double> > >(basis);

    FieldContainer<double> cvals;
    FieldContainer<double> bvals(basis->getCardinality(), basis->getCardinality());

    // Check exceptions.
#ifdef HAVE_INTREPID_DEBUG
    cvals.resize(1,2,3);
    INTREPID_TEST_COMMAND( coord_iface->getDofCoords(cvals), throwCounter, nException );
    cvals.resize(8,2);
    INTREPID_TEST_COMMAND( coord_iface->getDofCoords(cvals), throwCounter, nException );
    cvals.resize(9,1);
    INTREPID_TEST_COMMAND( coord_iface->getDofCoords(cvals), throwCounter, nException );
#endif
    cvals.resize(9,2);
    INTREPID_TEST_COMMAND( coord_iface->getDofCoords(cvals), throwCounter, nException ); nException--;
    // Check if number of thrown exceptions matches the one we expect
    if (throwCounter != nException) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
    }

    // Check mathematical correctness.
    basis->getValues(bvals, cvals, OPERATOR_VALUE);
    char buffer[120];
    for (int i=0; i<bvals.dimension(0); i++) {
      for (int j=0; j<bvals.dimension(1); j++) {
        if ((i != j) && (std::abs(bvals(i,j) - 0.0) > INTREPID_TOL)) {
          errorFlag++;
          sprintf(buffer, "\nValue of basis function %d at (%6.4e, %6.4e) is %6.4e but should be %6.4e!\n", i, cvals(i,0), cvals(i,1), bvals(i,j), 0.0);
          *outStream << buffer;
        }
        else if ((i == j) && (std::abs(bvals(i,j) - 1.0) > INTREPID_TOL)) {
          errorFlag++;
          sprintf(buffer, "\nValue of basis function %d at (%6.4e, %6.4e) is %6.4e but should be %6.4e!\n", i, cvals(i,0), cvals(i,1), bvals(i,j), 1.0);
          *outStream << buffer;
        }
      }
    }

  }
  catch (std::logic_error err){
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
