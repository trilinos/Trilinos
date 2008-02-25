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
  << "|    3) using LexContainer in DEBUG mode                                      |\n" \
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
  myContainer.resetContainer(indexRange);
  
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
  
  // Access same value by address
  int address = myContainer.getAddress(multiIndex);
  cout << "Access by address:     myContainer[" << address << "]=" << myContainer[address] <<"\n";
    
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
  myContainer.resetContainer(indexRange);
  
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
  
  // Access same value by address
  address = myNewContainer.getAddress(multiIndex);
  cout << "Access by address:     myNewContainer[" << address << "]=" << myNewContainer[address] <<"\n";

  cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 5: Debug mode                                                       |\n"\
    << "===============================================================================\n\n";
  
  // Trying to  get address using multi-index with the wrong rank (myContainer's rank is 4, 
  // whereas multiIndex has rank 5)
  cout \
    << "===============================================================================\n"\
    << " Trying to  get address using multi-index with the wrong rank: \n\n";
  try{
    myContainer.getAddress(multiIndex);
  }
  catch(std::logic_error err){
    cout << err.what() << "\n"; 
  }

  // Trying to get address using multi-index that is out of bounds (the last value should be < 6)
  cout \
    << "===============================================================================\n"\
    << " Trying to get address using multi-index that is out of bounds: \n\n";
  try{
    multiIndex[0] = 3; 
    multiIndex[1] = 1;
    multiIndex[2] = 2;
    multiIndex[3] = 2;
    multiIndex[4] = 6;
    myNewContainer.getAddress(multiIndex);
  }
  catch(std::logic_error err){
    cout << err.what() << "\n\n"; 
  }
  
  // Trying to set values from array whose size does not match LexContainer size
  cout \
    << "===============================================================================\n"\
    << " Trying to set values from array whose size does not match LexContainer size: \n\n";
  try{
    myContainer.setValues(dataTeuchosArray);
  }
  catch(std::logic_error err){
    cout << err.what() << "\n";
  }
  
  // Trying to use [] with address that is out of range (size of myNewContainer is 720)
  cout \
    << "===============================================================================\n"\
    << " Trying to use [] with address that is out of range: \n\n";
  try{
    myNewContainer[1000];
  }
  catch(std::logic_error err){
    cout << err.what() << "\n\n"; 
  }

  // Trying to create LexContainer using incompatible data array and indexRange. In this example
  // dataTeuchosArray corresponds to indexRange = {5,2,3,4,6} but we change the indexRange to one
  // that does not match the data. Note that if we permute the values in indexRange it will be
  // compatible with the data because it will specify the same size for the container. However, 
  // index bound permutation reshapes the container!
  cout \
    << "===============================================================================\n"\
    << " Trying to create LexContainer using incompatible data array and indexRange: \n\n";
  try{
    indexRange[0] = 5;
    indexRange[1] = 2;
    indexRange[2] = 3;
    indexRange[3] = 8;
    indexRange[4] = 6;
    
    LexContainer<double> myOtherContainer(indexRange, dataTeuchosArray);
  }
  catch(std::logic_error err){
    cout << err.what() << endl;
  }
  
   return 0;
}
