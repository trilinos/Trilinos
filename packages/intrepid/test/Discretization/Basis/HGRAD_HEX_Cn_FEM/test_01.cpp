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
\brief  Unit tests for the Intrepid::G_HEX_C2_FEM class.
\author Created by P. Bochev, D. Ridzal, and K. Peterson.
*/
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_HEX_Cn_FEM.hpp"
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
    << "|                 Unit Test (Basis_HGRAD_HEX_C2_FEM)                          |\n" \
    << "|                                                                             |\n" \
    << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n" \
    << "|     2) Basis values for VALUE, GRAD, and Dk operators                       |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Robert Kirby  (robert.c.kirby@ttu.edu),                |\n" \
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
  // get points
  const int deg=2;
  shards::CellTopology line(shards::getCellTopologyData< shards::Line<> >()); 
  FieldContainer<double> pts(PointTools::getLatticeSize(line,deg),1);
  PointTools::getLattice<double,FieldContainer<double> >(pts,line,deg);

  Basis_HGRAD_HEX_Cn_FEM<double, FieldContainer<double> > hexBasis(deg,deg,deg,pts,pts,pts);

  int errorFlag = 0;

  // Initialize throw counter for exception testing
  int nException     = 0;
  int throwCounter   = 0;

  // Define arrayS containing the 27 nodes of hexahedron<27> topology
  FieldContainer<double> hexNodes(27, 3);
  // do it lexicographically as a lattice
  hexNodes(0, 0) = -1.0;   hexNodes(0, 1) = -1.0;  hexNodes(0, 2) = -1.0;
  hexNodes(1, 0) =  0.0;   hexNodes(1, 1) = -1.0;  hexNodes(1, 2) = -1.0;
  hexNodes(2, 0) =  1.0;   hexNodes(2, 1) = -1.0;  hexNodes(2, 2) = -1.0;
  hexNodes(3, 0) = -1.0;   hexNodes(3, 1) =  0.0;  hexNodes(3, 2) = -1.0;
  hexNodes(4, 0) =  0.0;   hexNodes(4, 1) =  0.0;  hexNodes(4, 2) = -1.0;
  hexNodes(5, 0) =  1.0;   hexNodes(5, 1) =  0.0;  hexNodes(5, 2) = -1.0;
  hexNodes(6, 0) = -1.0;   hexNodes(6, 1) =  1.0;  hexNodes(6, 2) = -1.0;
  hexNodes(7, 0) = 0.0;    hexNodes(7, 1) =  1.0;  hexNodes(7, 2) = -1.0;
  hexNodes(8, 0) = 1.0;    hexNodes(8, 1) =  1.0;  hexNodes(8, 2) = -1.0;
  hexNodes(9, 0) = -1.0;   hexNodes(9, 1) = -1.0;  hexNodes(9, 2) = 0.0;
  hexNodes(10, 0) =  0.0;   hexNodes(10, 1) = -1.0;  hexNodes(10, 2) = 0.0;
  hexNodes(11, 0) =  1.0;   hexNodes(11, 1) = -1.0;  hexNodes(11, 2) = 0.0;
  hexNodes(12, 0) = -1.0;   hexNodes(12, 1) =  0.0;  hexNodes(12, 2) = 0.0;
  hexNodes(13, 0) =  0.0;   hexNodes(13, 1) =  0.0;  hexNodes(13, 2) = 0.0;
  hexNodes(14, 0) =  1.0;   hexNodes(14, 1) =  0.0;  hexNodes(14, 2) = 0.0;
  hexNodes(15, 0) = -1.0;   hexNodes(15, 1) =  1.0;  hexNodes(15, 2) = 0.0;
  hexNodes(16, 0) = 0.0;    hexNodes(16, 1) =  1.0;  hexNodes(16, 2) = 0.0;
  hexNodes(17, 0) = 1.0;    hexNodes(17, 1) =  1.0;  hexNodes(17, 2) = 0.0;
  hexNodes(18, 0) = -1.0;   hexNodes(18, 1) = -1.0;  hexNodes(18, 2) = 1.0;
  hexNodes(19, 0) =  0.0;   hexNodes(19, 1) = -1.0;  hexNodes(19, 2) = 1.0;
  hexNodes(20, 0) =  1.0;   hexNodes(20, 1) = -1.0;  hexNodes(20, 2) = 1.0;
  hexNodes(21, 0) = -1.0;   hexNodes(21, 1) =  0.0;  hexNodes(21, 2) = 1.0;
  hexNodes(22, 0) =  0.0;   hexNodes(22, 1) =  0.0;  hexNodes(22, 2) = 1.0;
  hexNodes(23, 0) =  1.0;   hexNodes(23, 1) =  0.0;  hexNodes(23, 2) = 1.0;
  hexNodes(24, 0) = -1.0;   hexNodes(24, 1) =  1.0;  hexNodes(24, 2) = 1.0;
  hexNodes(25, 0) = 0.0;    hexNodes(25, 1) =  1.0;  hexNodes(25, 2) = 1.0;
  hexNodes(26, 0) = 1.0;    hexNodes(26, 1) =  1.0;  hexNodes(26, 2) = 1.0;

 
  // Generic array for the output values; needs to be properly resized depending on the operator type
  FieldContainer<double> vals;

  try{
    // exception #1: CURL cannot be applied to scalar functions in 3D
    // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
    vals.resize(hexBasis.getCardinality(), hexNodes.dimension(0), 4 );
    INTREPID_TEST_COMMAND( hexBasis.getValues(vals, hexNodes, OPERATOR_CURL), throwCounter, nException );

    // exception #2: DIV cannot be applied to scalar functions in 3D
    // resize vals to rank-2 container with dimensions (num. basis functions, num. points)
    vals.resize(hexBasis.getCardinality(), hexNodes.dimension(0) );
    INTREPID_TEST_COMMAND( hexBasis.getValues(vals, hexNodes, OPERATOR_DIV), throwCounter, nException );
        
    // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    // exception #3
    INTREPID_TEST_COMMAND( hexBasis.getDofOrdinal(3,10,0), throwCounter, nException );    
    // exception #4
    INTREPID_TEST_COMMAND( hexBasis.getDofOrdinal(1,2,1), throwCounter, nException );
    // exception #5
    INTREPID_TEST_COMMAND( hexBasis.getDofOrdinal(0,4,1), throwCounter, nException );
    // exception #6
    INTREPID_TEST_COMMAND( hexBasis.getDofTag(28), throwCounter, nException );
    // exception #7
    INTREPID_TEST_COMMAND( hexBasis.getDofTag(-1), throwCounter, nException );

#ifdef HAVE_INTREPID_DEBUG 
    // Exceptions 8-18 test exception handling with incorrectly dimensioned input/output arrays
    // exception #8: input points array must be of rank-2
    FieldContainer<double> badPoints1(4, 5, 3);
    INTREPID_TEST_COMMAND( hexBasis.getValues(vals, badPoints1, OPERATOR_VALUE), throwCounter, nException );

    // exception #9 dimension 1 in the input point array must equal space dimension of the cell
    FieldContainer<double> badPoints2(4, hexBasis.getBaseCellTopology().getDimension() - 1);
    INTREPID_TEST_COMMAND( hexBasis.getValues(vals, badPoints2, OPERATOR_VALUE), throwCounter, nException );

    // exception #10 output values must be of rank-2 for OPERATOR_VALUE
    FieldContainer<double> badVals1(4, 3, 1);
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals1, hexNodes, OPERATOR_VALUE), throwCounter, nException );
 
    // exception #11 output values must be of rank-3 for OPERATOR_GRAD
    FieldContainer<double> badVals2(4, 3);
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals2, hexNodes, OPERATOR_GRAD), throwCounter, nException );
    
    // exception #12 output values must be of rank-3 for OPERATOR_D1
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals2, hexNodes, OPERATOR_D1), throwCounter, nException );
    
    // exception #13 output values must be of rank-3 for OPERATOR_D2
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals2, hexNodes, OPERATOR_D2), throwCounter, nException );
    
    // exception #14 incorrect 0th dimension of output array (must equal number of basis functions)
    FieldContainer<double> badVals3(hexBasis.getCardinality() + 1, hexNodes.dimension(0));
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals3, hexNodes, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #15 incorrect 1st dimension of output array (must equal number of points)
    FieldContainer<double> badVals4(hexBasis.getCardinality(), hexNodes.dimension(0) + 1);
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals4, hexNodes, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
    FieldContainer<double> badVals5(hexBasis.getCardinality(), hexNodes.dimension(0), hexBasis.getBaseCellTopology().getDimension() - 1);
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals5, hexNodes, OPERATOR_GRAD), throwCounter, nException );
    
    // exception #17: incorrect 2nd dimension of output array (must equal D2 cardinality in 3D)
    FieldContainer<double> badVals6(hexBasis.getCardinality(), hexNodes.dimension(0), 40);
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals6, hexNodes, OPERATOR_D2), throwCounter, nException );
    
    // exception #18: incorrect 2nd dimension of output array (must equal D3 cardinality in 3D)
    FieldContainer<double> badVals7(hexBasis.getCardinality(), hexNodes.dimension(0), 50);
    INTREPID_TEST_COMMAND( hexBasis.getValues(badVals7, hexNodes, OPERATOR_D3), throwCounter, nException );
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
  
  // VALUE: Each row gives the 27 correct basis set values at an evaluation point
  double basisValues[] = {
    1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, 0, \
    0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000, 0, \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.000  };
  
  
  // GRAD, D1, D2, D3 and D4 test values are stored in files due to their large size
  std::string     fileName;
  std::ifstream   dataFile;

  // GRAD and D1 values are stored in (F,P,D) format in a data file. Read file and do the test
  std::vector<double> basisGrads;           // Flat array for the gradient values.
  
  fileName = "./testdata/HEX_C2_GradVals.dat";
  dataFile.open(fileName.c_str());
  TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                      ">>> ERROR (HGRAD_HEX_C2/test01): could not open GRAD values data file, test aborted.");
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
  fileName = "./testdata/HEX_C2_D2Vals.dat";  
  dataFile.open(fileName.c_str());
  TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                      ">>> ERROR (HGRAD_HEX_C2/test01): could not open D2 values data file, test aborted.");
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
  
  fileName = "./testdata/HEX_C2_D3Vals.dat";  
  dataFile.open(fileName.c_str());
  TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                      ">>> ERROR (HGRAD_HEX_C2/test01): could not open D3 values data file, test aborted.");
  
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
 
  
  //D4: flat array with the values of D applied to basis functions. Multi-index is (F,P,D4cardinality)
  std::vector<double> basisD4;
  
  fileName = "./testdata/HEX_C2_D4Vals.dat";
  dataFile.open(fileName.c_str());
  TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                      ">>> ERROR (HGRAD_HEX_C2/test01): could not open D4 values data file, test aborted.");
  
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
    int numFields = hexBasis.getCardinality();
    int numPoints = hexNodes.dimension(0);
    int spaceDim  = hexBasis.getBaseCellTopology().getDimension();
    
    // Generic array for values, grads, curls, etc. that will be properly sized before each call
    FieldContainer<double> vals;
    
    // Check VALUE of basis functions: resize vals to rank-2 container:
    vals.resize(numFields, numPoints);
    hexBasis.getValues(vals, hexNodes, OPERATOR_VALUE);
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
    hexBasis.getValues(vals, hexNodes, OPERATOR_GRAD);
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
    hexBasis.getValues(vals, hexNodes, OPERATOR_D1);
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
    hexBasis.getValues(vals, hexNodes, OPERATOR_D2);
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
    hexBasis.getValues(vals, hexNodes, OPERATOR_D3);
     
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
    hexBasis.getValues(vals, hexNodes, OPERATOR_D4);
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
    


  // Check D7 to D10 - must be zero. This basis does not cover D5 and D6
    for(EOperator op = OPERATOR_D7; op < OPERATOR_MAX; op++) {
      
      // The last dimension is the number of kth derivatives and needs to be resized for every Dk
      int DkCardin  = Intrepid::getDkCardinality(op, spaceDim);
      vals.resize(numFields, numPoints, DkCardin);    

      hexBasis.getValues(vals, hexNodes, op);
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
