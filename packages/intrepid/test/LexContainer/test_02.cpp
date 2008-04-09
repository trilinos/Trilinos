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


/** \file
\brief  Unit test of LexContainer class
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_LexContainer.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

using namespace Intrepid;

int main(int argc, char *argv[]) {
  
  // This little trick lets us print to cout  only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original cout .
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream  \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                           Unit Test LexContainer                            |\n" \
    << "|                                                                             |\n" \
    << "|     1) Testing exception handling                                           |\n" \
    << "|       requires intrepid to be configured with --enable-intrepid-debug       |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  
  // Test initializations
  int errorFlag  = 0;
  
  // Define variables to create and use LexContainers
  Teuchos::Array<int> indexRange;
  Teuchos::Array<int> multiIndex;
  
  // Initialize indexRange for rank-4 multi-index value
  indexRange.resize(4);
  indexRange[0] = 5;
  indexRange[1] = 3;
  indexRange[2] = 2;
  indexRange[3] = 7;
  
  // Create a LexContainer
  LexContainer<double> myContainer(indexRange);
  
  // These tests should only run if intrepid was configured with --enable-intrepid-debug option
  // Each test is designed to cause an exception. The total number of all caught exceptions should
  // be the same as the number of tests in this section.
#ifdef HAVE_INTREPID_DEBUG
  
  *outStream << "\n" \
    << "===============================================================================\n"\
    << "| TEST 1: Catching exceptions                                                 |\n"\
    << "===============================================================================\n\n";
  
  int numTestException = 7;
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + numTestException;
  
  try{ // Outer try block contains all tests for exception
    
    
    try{ // catch exception (1)
      
      //  Trying to  get enumeration using multi-index with the wrong rank:
      *outStream << "\n" \
      << "===============================================================================\n"\
      << "  Trying to  get enumeration using multi-index with the wrong rank: \n";
      multiIndex.resize(5);
      multiIndex[0] = 3; 
      multiIndex[1] = 1;
      multiIndex[2] = 2;
      multiIndex[3] = 2;
      multiIndex[4] = 6;
      myContainer.getEnumeration(multiIndex);
    }
    catch (std::logic_error err) {
      *outStream  << err.what() << "\n";
    };
    
    
    
    try{ // catch exception (2)
      
      // Trying to get enumeration using multi-index that is out of bounds: 3rd index is 4, must be <2
      *outStream << "\n" \
      << "===============================================================================\n"\
      << " Trying to get enumeration using multi-index that is out of bounds: \n";
      multiIndex.resize(4);
      multiIndex[0] = 3; 
      multiIndex[1] = 1;
      multiIndex[2] = 4;
      multiIndex[3] = 2;
      myContainer.getEnumeration(multiIndex);
    }
    catch (std::logic_error err) {
      *outStream  << err.what() << "\n";
    };
    
 
    
    try{  // catch exception (3)
      
      // Trying to set values from array whose size is less than LexContainer's size:
      *outStream << "\n" \
      << "===============================================================================\n"\
      <<  " Trying to set values from array whose size is less than LexContainer's size: \n";
      
      // Change one of the values of the indexRange to a lesser value: original value was 5 
      indexRange[0] = 4;
      
      // Define Teuchos::Array to store values with dimension equal to the number of multi-indexed values
      Teuchos::Array<double> dataTeuchosArray(4*3*2*7);
      
      // Fill with data
      int counter = 0;
      for(int i=0; i < indexRange[0]; i++){
        for(int j=0; j < indexRange[1]; j++){
          for(int k=0; k < indexRange[2]; k++){
            for(int l = 0; l < indexRange[3]; l++){
              dataTeuchosArray[counter] = (double)counter;
              counter++;
            }
          }
        }
      }
      
      // Now try to stuff this data into LexContainer
      myContainer.setValues(dataTeuchosArray);
    }
    catch (std::logic_error err) {
      *outStream  << err.what() << "\n";
    };
    
    
    
    try{  // catch exception (4)
      
      // Trying to set values from array whose size is greater than LexContainer's size:
      *outStream << "\n" \
      << "===============================================================================\n"\
      <<  " Trying to set values from array whose size is greater than LexContainer's size: \n";
      
      // Change one of the values of the indexRange to a lesser value: restore indexRange[0] to the 
      // value used to construct the LexArray and change indexRange[2] to a greater value
      indexRange[0] = 5;
      indexRange[2] = 3;
      
      // Define Teuchos::Array to store values with dimension equal to the number of multi-indexed values
      Teuchos::Array<double> dataTeuchosArray(5*3*3*7);
      
      // Fill with data
      int counter = 0;
      for(int i=0; i < indexRange[0]; i++){
        for(int j=0; j < indexRange[1]; j++){
          for(int k=0; k < indexRange[2]; k++){
            for(int l = 0; l < indexRange[3]; l++){
              dataTeuchosArray[counter] = (double)counter;
              counter++;
            }
          }
        }
      }
      
      // Now try to stuff this data into LexContainer
      myContainer.setValues(dataTeuchosArray);
    }
    catch (std::logic_error err) {
      *outStream  << err.what() << "\n";
    };
    
    
    
    try{ // catch exception (5)
      
      // Trying to use [] with enumeration that is out of range:
      *outStream << "\n" \
      << "===============================================================================\n"\
      << " Trying to use [] with enumeration that is out of range: \n";
      myContainer[1000];
    }
    catch (std::logic_error err) {
      *outStream  << err.what() << "\n";
    }
    
    
    
    try{ // catch exception (6)
      
      // Trying to get multi-index from enumeration that is out of bounds:
      *outStream << "\n" \
      << "===============================================================================\n"\
      << " Trying to get multi-index from enumeration that is out of bounds: \n";
      myContainer.getMultiIndex(multiIndex,10000);
    }
    catch(std::logic_error err) {
     *outStream << err.what() << "\n";
    }
    
  
    
    try{ // catch exception (7)
      
      //Trying to self-assign LexContainer
      *outStream << "\n" \
      << "===============================================================================\n"\
      << " Trying to self-assign LexContainer \n";
      myContainer = myContainer;
    }
    catch(std::logic_error err) {
      *outStream << err.what() << "\n"; 
    }

    // Check if number of caught exceptions matches the expected number
    if (TestForException_getThrowNumber() != endThrowNumber) {
      errorFlag++;
    }
  } // outer try block
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream  << err.what() << "\n";
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  }
#endif
  
  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return errorFlag;
}
