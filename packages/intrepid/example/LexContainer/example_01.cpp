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
\brief  Illustrates use of the LexContainer class.
\author Created by P. Bochev and D. Ridzal
*/

#include "Intrepid_LexContainer.hpp"


using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {
  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                  Example use of the LexContainer class                      |\n" \
  << "|                                                                             |\n" \
  << "|    1) Creating and filling LexContainer objects                             |\n" \
  << "|    2) Accessing elements in LexContainer objects                            |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n\n";
  
  // Define variables to create and use LexContainers
  Teuchos::Array<int> indexRange;
  Teuchos::Array<int> multiIndex;
  
  cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 1: rank 2 multi-index: {u(p,i) | 0 <= p < 5; 0 <= i < 3 }           |\n"\
    << "===============================================================================\n\n";
  // This rank 2 multi-indexed value can be used to store values of 3D vector field
  // evaluated at 5 points, or values of gradient of scalar field evaluated at the points
  
  // Resize indexRange and multiIndex for rank-2 multi-indexed value 
  indexRange.resize(2);
  multiIndex.resize(2);
  
  // Load upper ranges for the two indices in the multi-indexed value
  indexRange[0] = 5;
  indexRange[1] = 3;
    
  // Create LexContainer that can hold the rank-2 multi-indexed value
  LexContainer<double> myContainer(indexRange);
  
  // Fill with some data: leftmost index changes last, rightmost index changes first!  
  for(int p = 0; p < indexRange[0]; p++){
    multiIndex[0] = p;
    
    for(int i = 0; i < indexRange[1]; i++){
      multiIndex[1] = i;
      
      // Load value with multi-index {p,i} to container
      myContainer.setValue((double)(i+p), multiIndex);
    }
  }
  
  // Show container contents
  cout << myContainer;
  
  cout << "\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 2: rank 3 multi-index: {u(p,i,j) | 0 <=p< 5; 0 <= i<2, 0<=j<3       |\n"\
  << "===============================================================================\n\n";
  // This rank-3 value can be used to store subset of second partial derivatives  values
  // of a scalar function at p points.
  
  // Resize indexRange and multiIndex for rank-3 multi-indexed value 
  indexRange.resize(3);
  multiIndex.resize(3);
  
  // Define upper ranges for the three indices in the multi-indexed value
  indexRange[0] = 5;
  indexRange[1] = 2;
  indexRange[2] = 3;
  
  // Reset the existing container to accept rank-3 value with the specified index ranges
  myContainer.resize(indexRange);
  
  // Fill with some data
  for(int p = 0; p < indexRange[0]; p++){
    multiIndex[0] = p;
    for(int i = 0; i < indexRange[1]; i++){
      multiIndex[1] = i;
      for(int j = 0; j < indexRange[2]; j++){
        multiIndex[2] = j;
        
        // Load value with multi-index {p,i} to container
        myContainer.setValue((double)(p+i+j), multiIndex);
      }
    }
  }
  
  // Display contents
  cout << myContainer;
  
  // Access values by multi-index
  multiIndex[0] = 3; 
  multiIndex[1] = 1;
  multiIndex[2] = 2;
  
  cout << "Access by multi-index: myContainer{" \
    << multiIndex[0] << multiIndex[1] << multiIndex[2] << "}="\
    << myContainer.getValue(multiIndex) <<"\n";
  
  // Access same value by enumeration
  int enumeration = myContainer.getEnumeration(multiIndex);
  cout << "Access by enumeration:     myContainer[" << enumeration << "]=" << myContainer[enumeration] <<"\n";
    
  cout << "\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 3: setting rank-4 values to existing LexContainer from data array   |\n"\
  << "===============================================================================\n\n";
  
  // Initialize indexRange for rank-4 multi-index value
  indexRange.resize(4);
  indexRange[0] = 2;
  indexRange[1] = 3;
  indexRange[2] = 2;
  indexRange[3] = 4;
  
  // Resize myContainer using the new indexRange
  myContainer.resize(indexRange);
  
  // Define array to store values with dimension equal to the number of multi-indexed values
  double dataArray[2*3*2*4];
  
  // Fill with data
  int counter = 0;
  for(int i=0; i < indexRange[0]; i++){
    for(int j=0; j < indexRange[1]; j++){
      for(int k=0; k < indexRange[2]; k++){
        for(int l=0; l < indexRange[3]; l++){
          dataArray[counter] = (double)counter;
          counter++;
        }
      }
    }
  }
  
  // Fill myContainer with values from the array and show it. data size = counter!
  myContainer.setValues(dataArray,counter);
  cout << myContainer;
  
  cout << "\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 4: making rank-5 LexContainer from data array and index range array |\n"\
  << "===============================================================================\n\n";
  // Initialize indexRange for rank-5 multi-index value
  indexRange.resize(5);
  indexRange[0] = 5;
  indexRange[1] = 2;
  indexRange[2] = 3;
  indexRange[3] = 4;
  indexRange[4] = 6;
  
  // Define Teuchos::Array to store values with dimension equal to the number of multi-indexed values
  Teuchos::Array<double> dataTeuchosArray(5*2*3*4*6);
  
  // Fill with data
  counter = 0;
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
  
  // Create LexContainer from data array and index array and show it
  LexContainer<double> myNewContainer(indexRange, dataTeuchosArray);
  cout << myNewContainer;

  // Access values by multi-index
  multiIndex.resize(myNewContainer.getRank());
  multiIndex[0] = 3; 
  multiIndex[1] = 1;
  multiIndex[2] = 2;
  multiIndex[3] = 2;
  multiIndex[4] = 5;
  
  cout << "Access by multi-index: myNewContainer{" \
    << multiIndex[0] << multiIndex[1] << multiIndex[2] << multiIndex[3] << multiIndex[4] << "}="\
    << myNewContainer.getValue(multiIndex) <<"\n";
  
  // Access same value by enumeration
  enumeration = myNewContainer.getEnumeration(multiIndex);
  cout << "Access by enumeration: myNewContainer[" << enumeration << "]=" << myNewContainer[enumeration] <<"\n";

  cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 5: making trivial LexContainers and storing a single zero           |\n"\
    << "===============================================================================\n\n";
  
  // Make trivial container by resetting the index range to zero rank (no indices) and then
  // using resize method
  indexRange.resize(0);
  myContainer.resize(indexRange);
  std::cout << myContainer;
  
  // Make trivial container by using empty method:
  myNewContainer.empty();
  std::cout << myNewContainer;
  
  // Now use storeZero() to reset the container to hold a single zero
  myNewContainer.storeZero();
  std::cout << myNewContainer;
  
  
   return 0;
}
