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

/**   \file   test_01.cpp
      \brief  Unit tests for the Intrepid::HGRAD_WEDGE_C2_FEM class.
      \author Created by P. Bochev, D. Ridzal, and K. Peterson.
*/
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_WEDGE_C2_FEM.hpp"
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
    << "|                 Unit Test (Basis_HGRAD_WEDGE_C2_FEM)                        |\n" \
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
  Basis_HGRAD_WEDGE_C2_FEM<double, FieldContainer<double> > wedgeBasis;
  int errorFlag = 0;

  // Initialize throw counter for exception testing
  int nException     = 0;
  int throwCounter   = 0;

  // Nodes of Wedde<18>: vertices, edge midpoints, Quadrilateral face centers 
  FieldContainer<double> wedgeNodes(18, 3);
  wedgeNodes(0,0) =  0.0;  wedgeNodes(0,1) =  0.0;  wedgeNodes(0,2) = -1.0;  
  wedgeNodes(1,0) =  1.0;  wedgeNodes(1,1) =  0.0;  wedgeNodes(1,2) = -1.0;  
  wedgeNodes(2,0) =  0.0;  wedgeNodes(2,1) =  1.0;  wedgeNodes(2,2) = -1.0;
  wedgeNodes(3,0) =  0.0;  wedgeNodes(3,1) =  0.0;  wedgeNodes(3,2) =  1.0;  
  wedgeNodes(4,0) =  1.0;  wedgeNodes(4,1) =  0.0;  wedgeNodes(4,2) =  1.0;  
  wedgeNodes(5,0) =  0.0;  wedgeNodes(5,1) =  1.0;  wedgeNodes(5,2) =  1.0;
    
  wedgeNodes(6,0) =  0.5;  wedgeNodes(6,1) =  0.0;  wedgeNodes(6,2) = -1.0;  
  wedgeNodes(7,0) =  0.5;  wedgeNodes(7,1) =  0.5;  wedgeNodes(7,2) = -1.0;  
  wedgeNodes(8,0) =  0.0;  wedgeNodes(8,1) =  0.5;  wedgeNodes(8,2) = -1.0;
  wedgeNodes(9,0) =  0.0;  wedgeNodes(9,1) =  0.0;  wedgeNodes(9,2) =  0.0;
  wedgeNodes(10,0)=  1.0;  wedgeNodes(10,1)=  0.0;  wedgeNodes(10,2)=  0.0;  
  wedgeNodes(11,0)=  0.0;  wedgeNodes(11,1)=  1.0;  wedgeNodes(11,2)=  0.0;  
 
  wedgeNodes(12,0)=  0.5;  wedgeNodes(12,1)=  0.0;  wedgeNodes(12,2)=  1.0;  
  wedgeNodes(13,0)=  0.5;  wedgeNodes(13,1)=  0.5;  wedgeNodes(13,2)=  1.0;  
  wedgeNodes(14,0)=  0.0;  wedgeNodes(14,1)=  0.5;  wedgeNodes(14,2)=  1.0;  
  wedgeNodes(15,0)=  0.5;  wedgeNodes(15,1)=  0.0;  wedgeNodes(15,2)=  0.0;  
  wedgeNodes(16,0)=  0.5;  wedgeNodes(16,1)=  0.5;  wedgeNodes(16,2)=  0.0;  
  wedgeNodes(17,0)=  0.0;  wedgeNodes(17,1)=  0.5;  wedgeNodes(17,2)=  0.0;  


  // Generic array for the output values; needs to be properly resized depending on the operator type
  FieldContainer<double> vals;

  try{
    // exception #1: CURL cannot be applied to scalar functions
    // resize vals to rank-3 container with dimensions (num. points, num. basis functions)
    vals.resize(wedgeBasis.getCardinality(), wedgeNodes.dimension(0), 3 );
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_DIV), throwCounter, nException );

    // exception #2: DIV cannot be applied to scalar functions
    // resize vals to rank-2 container with dimensions (num. points, num. basis functions)
    vals.resize(wedgeBasis.getCardinality(), wedgeNodes.dimension(0) );
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_DIV), throwCounter, nException );
        
    // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    // exception #3
    INTREPID_TEST_COMMAND( wedgeBasis.getDofOrdinal(3,0,0), throwCounter, nException );
    // exception #4
    INTREPID_TEST_COMMAND( wedgeBasis.getDofOrdinal(1,1,1), throwCounter, nException );
    // exception #5
    INTREPID_TEST_COMMAND( wedgeBasis.getDofOrdinal(0,9,0), throwCounter, nException );
    // exception #6
    INTREPID_TEST_COMMAND( wedgeBasis.getDofTag(18), throwCounter, nException );
    // exception #7
    INTREPID_TEST_COMMAND( wedgeBasis.getDofTag(-1), throwCounter, nException );
    
#ifdef HAVE_INTREPID_DEBUG
    // Exceptions 8-18 test exception handling with incorrectly dimensioned input/output arrays
    // exception #8: input points array must be of rank-2
    FieldContainer<double> badPoints1(4, 5, 3);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(vals, badPoints1, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #9 dimension 1 in the input point array must equal space dimension of the cell
    FieldContainer<double> badPoints2(4, wedgeBasis.getBaseCellTopology().getDimension() + 1);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(vals, badPoints2, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #10 output values must be of rank-2 for OPERATOR_VALUE
    FieldContainer<double> badVals1(4, 3, 1);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals1, wedgeNodes, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #11 output values must be of rank-3 for OPERATOR_GRAD
    FieldContainer<double> badVals2(4, 3);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals2, wedgeNodes, OPERATOR_GRAD), throwCounter, nException );
    
    // exception #12 output values must be of rank-3 for OPERATOR_D1
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals2, wedgeNodes, OPERATOR_D1), throwCounter, nException );
    
    // exception #13 output values must be of rank-3 for OPERATOR_D2
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals2, wedgeNodes, OPERATOR_D2), throwCounter, nException );
    
    // exception #14 incorrect 0th dimension of output array (must equal number of basis functions)
    FieldContainer<double> badVals3(wedgeBasis.getCardinality() + 1, wedgeNodes.dimension(0));
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals3, wedgeNodes, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #15 incorrect 1st dimension of output array (must equal number of points)
    FieldContainer<double> badVals4(wedgeBasis.getCardinality(), wedgeNodes.dimension(0) + 1);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals4, wedgeNodes, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
    FieldContainer<double> badVals5(wedgeBasis.getCardinality(), wedgeNodes.dimension(0), wedgeBasis.getBaseCellTopology().getDimension() - 1);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals5, wedgeNodes, OPERATOR_GRAD), throwCounter, nException );
    
    // exception #17: incorrect 2nd dimension of output array (must equal D2 cardinality in 3D)
    FieldContainer<double> badVals6(wedgeBasis.getCardinality(), wedgeNodes.dimension(0), 40);
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals6, wedgeNodes, OPERATOR_D2), throwCounter, nException );
    
    // exception #18: incorrect 2nd dimension of output array (must equal D3 cardinality in 3D)
    INTREPID_TEST_COMMAND( wedgeBasis.getValues(badVals6, wedgeNodes, OPERATOR_D3), throwCounter, nException );
#endif

  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };
  
  // Check if number of thrown exceptions matches the one we expect - 18
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
  
  // VALUE: correct basis function values in (F,P) format
  double basisValues[] = {
    1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00  };
  
  // GRAD, D1, D2, D3 and D4 test values are stored in files due to their large size
  std::string     fileName;
  std::ifstream   dataFile;
  
  // GRAD and D1 values are stored in (F,P,D) format in a data file. Read file and do the test
  std::vector<double> basisGrads;           // Flat array for the gradient values.
  
  fileName = "./testdata/WEDGE_C2_GradVals.dat";
  dataFile.open(fileName.c_str());
  TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                      ">>> ERROR (HGRAD_WEDGE_C2/test01): could not open GRAD values data file, test aborted.");
  while (!dataFile.eof() ){
    double temp;
    string line;                            // string for one line of input file
    std::getline(dataFile, line);           // get next line from file
    stringstream data_line(line);           // convert to stringstream
    while(data_line >> temp){               // extract value from line
      basisGrads.push_back(temp);           // push into vector
    }
  }
  // It turns out that just closing and then opening the ifstream variable does not reset it
  // and subsequent open() command fails. One fix is to explicitely clear the ifstream, or
  // scope the variables.
  dataFile.close();
  dataFile.clear();
  
  
  //D2: flat array with the values of D2 applied to basis functions. Multi-index is (F,P,D2cardinality)
  std::vector<double> basisD2;
  
  fileName = "./testdata/WEDGE_C2_D2Vals.dat";
  dataFile.open(fileName.c_str());
  TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                      ">>> ERROR (HGRAD_WEDGE_C2/test01): could not open D2 values data file, test aborted.");
  
  while (!dataFile.eof() ){
    double temp;
    string line;                            // string for one line of input file
    std::getline(dataFile, line);           // get next line from file
    stringstream data_line(line);           // convert to stringstream
    while(data_line >> temp){               // extract value from line
      basisD2.push_back(temp);              // push into vector
    }
  }
  dataFile.close();
  dataFile.clear();
  
  
  //D3: flat array with the values of D3 applied to basis functions. Multi-index is (F,P,D3cardinality)
  std::vector<double> basisD3;
  
  fileName = "./testdata/WEDGE_C2_D3Vals.dat";
  dataFile.open(fileName.c_str());
  TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                      ">>> ERROR (HGRAD_WEDGE_C2/test01): could not open D3 values data file, test aborted.");
  
  while (!dataFile.eof() ){
    double temp;
    string line;                            // string for one line of input file
    std::getline(dataFile, line);           // get next line from file
    stringstream data_line(line);           // convert to stringstream
    while(data_line >> temp){               // extract value from line
      basisD3.push_back(temp);              // push into vector
    }
  }
  dataFile.close();
  dataFile.clear();
  
  
  //D4: flat array with the values of D3 applied to basis functions. Multi-index is (F,P,D4cardinality)
  std::vector<double> basisD4;
  
  fileName = "./testdata/WEDGE_C2_D4Vals.dat";
  dataFile.open(fileName.c_str());
  TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                      ">>> ERROR (HGRAD_WEDGE_C2/test01): could not open D4 values data file, test aborted.");
  
  while (!dataFile.eof() ){
    double temp;
    string line;                            // string for one line of input file
    std::getline(dataFile, line);           // get next line from file
    stringstream data_line(line);           // convert to stringstream
    while(data_line >> temp){               // extract value from line
      basisD4.push_back(temp);              // push into vector
    }
  }
  dataFile.close();
  dataFile.clear();

  
  try{
        
    // Dimensions for the output arrays:
    int numFields = wedgeBasis.getCardinality();
    int numPoints = wedgeNodes.dimension(0);
    int spaceDim  = wedgeBasis.getBaseCellTopology().getDimension();
    
    // Generic array for values, grads, curls, etc. that will be properly sized before each call
    FieldContainer<double> vals;
    
    // Check VALUE of basis functions: resize vals to rank-2 container:
    vals.resize(numFields, numPoints);
    wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_VALUE);
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
    wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_GRAD);
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
    wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_D1);
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
    wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_D2);
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
    wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_D3);
    
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
              << " but reference D3 component: " << basisD3[l] << "\n";
          }
        }
      }
    }
    
    
    // Check D4 of basis function
    int D4cardinality = Intrepid::getDkCardinality(OPERATOR_D4, spaceDim);
    vals.resize(numFields, numPoints, D4cardinality);    
    wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_D4);
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

      wedgeBasis.getValues(vals, wedgeNodes, op);
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
