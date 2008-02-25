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
  
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
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
    << "|                           Unit Test LexContainer                            |\n" \
    << "|                                                                             |\n" \
    << "|     1) Value accesss by multi-index and address, setting values             |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  
  int errorFlag  = 0;
  
  // Define variables needed to construct and use rank-5 LexContainer:
  Teuchos::Array<int> indexRange(5);
  Teuchos::Array<int> multiIndex(5);
    
  // Set upper index ranges
  indexRange[0] = 3;
  indexRange[1] = 7;
  indexRange[2] = 11;
  indexRange[3] = 5;
  indexRange[4] = 6;
  
  // Define Teuchos::Array with dimension equal to LexContainer capacity (determined from indexRange)
  int dataSize = indexRange[0];
  for(int r = 1; r < (int)indexRange.size(); r++){
    dataSize *= indexRange[r];
  }
  Teuchos::Array<double> dataTeuchosArray(dataSize);
  
  // Fill with data
  int counter = 0;
  for(int i=0; i < indexRange[0]; i++){
    for(int j=0; j < indexRange[1]; j++){
      for(int k=0; k < indexRange[2]; k++){
        for(int l = 0; l < indexRange[3]; l++){
          for(int m = 0; m < indexRange[4]; m++){
            dataTeuchosArray[counter] = (double)counter;
            counter++;
          }
        }
      }
    }
  }
  
  try{
    
    // Create LexContainer from data array and index array and show it
    LexContainer<double> myNewContainer(indexRange, dataTeuchosArray);
    *outStream << myNewContainer;
    
    // Initialize some variables
    int containerSize = myNewContainer.getSize();
    int containerRank = myNewContainer.getRank();
    Teuchos::Array<int> multiIndex(containerRank);

    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 1: Traverse LexContainer by address                                    |\n"\
      << "===============================================================================\n";
    
    // Loop over container by address
    for(int address = 0; address < containerSize; address++) {
      
      // Convert address to multi-index
      myNewContainer.getMultiIndex(multiIndex, address);
      
      // Check if access by address gives the same value as access by multi-index
      if( myNewContainer[address] != myNewContainer.getValue(multiIndex) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 2: Traverse LexContainer by multi-index & compare with data array      |\n"\
      << "===============================================================================\n";
    
    // Loop over container by multi-index
    for(int i=0; i < indexRange[0]; i++){
      multiIndex[0] = i;
      for(int j=0; j < indexRange[1]; j++){
        multiIndex[1] = j;
        for(int k=0; k < indexRange[2]; k++){
          multiIndex[2] = k; 
          for(int l = 0; l < indexRange[3]; l++){
            multiIndex[3] = l;
            for(int m = 0; m < indexRange[4]; m++){
              multiIndex[4] = m;
              
              // Convert multi-index to address
              int address = myNewContainer.getAddress(multiIndex);
              
              // Check if access by address gives the same value as access by multi-index
              if( myNewContainer[address] != myNewContainer.getValue(multiIndex)) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              }
              
              // Check if container data matches values in dataTeuchosArray 
              if(dataTeuchosArray[address] != myNewContainer.getValue(multiIndex)) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              }
            }
          }
        }
      }
    }
    
  } //try 
  
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
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
