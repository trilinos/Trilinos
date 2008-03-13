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
  << "|    1) using LexContainer in DEBUG mode:                                     |\n" \
  << "|       requires intrepid to be configured with --enable-intrepid-debug       |\n" \
  << "|       See /test/LexContainer/test_02.cpp for more examples                  |\n" \
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
  
  // Initialize indexRange for rank-4 multi-index value
  indexRange.resize(4);
  indexRange[0] = 5;
  indexRange[1] = 3;
  indexRange[2] = 2;
  indexRange[3] = 7;
  
  // Create a LexContainer
  LexContainer<double> myContainer(indexRange);
  
  cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 1: Debug mode                                                       |\n"\
    << "===============================================================================\n\n";
  
  // Trying to  get enumeration using multi-index with the wrong rank (myContainer's rank is 4, 
  // whereas multiIndex has rank 5)
  cout \
    << "===============================================================================\n"\
    << " Trying to  get enumeration using multi-index with the wrong rank: \n\n";
  try{
    multiIndex.resize(5);
    multiIndex[0] = 3; 
    multiIndex[1] = 1;
    multiIndex[2] = 2;
    multiIndex[3] = 2;
    multiIndex[4] = 6;
    myContainer.getEnumeration(multiIndex);
  }
  catch(std::logic_error err){
    cout << err.what() << "\n"; 
  }
  
  // Trying to get enumeration using multi-index that is out of bounds: 3rd index is 4, must be <2 
  cout \
    << "===============================================================================\n"\
    << " Trying to get enumeration using multi-index that is out of bounds: \n\n";
  try{
    multiIndex.resize(4);
    multiIndex[0] = 3; 
    multiIndex[1] = 1;
    multiIndex[2] = 4;
    multiIndex[3] = 2;
    myContainer.getEnumeration(multiIndex);
  }
  catch(std::logic_error err){
    cout << err.what() << "\n\n"; 
  }
  
  // Trying to set values from array whose size is less than LexContainer size
  cout \
    << "===============================================================================\n"\
    << " Trying to set values from array whose size is less than LexContainer's size: \n\n";

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
  try{
    myContainer.setValues(dataTeuchosArray);
  }
  catch(std::logic_error err){
    cout << err.what() << "\n";
  }
  
  // Trying to set values from array whose size is greater than LexContainer's size
  cout \
    << "===============================================================================\n"\
    << " Trying to set values from array whose size is greater than LexContainer's size: \n\n";
  // Change one of the values of the indexRange to a lesser value: restore indexRange[0] to the 
  // value used to construct the LexArray and change indexRange[2] to a greater value
  indexRange[0] = 5;
  indexRange[2] = 3;
  
  // Define Teuchos::Array to store values with dimension equal to the number of multi-indexed values
  dataTeuchosArray.resize(5*3*3*7);
  
  // Fill with data
  counter = 0;
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
  try{
    myContainer.setValues(dataTeuchosArray);
  }
  catch(std::logic_error err){
    cout << err.what() << "\n";
  }
  
  
  // Trying to use [] with enumeration that is out of range (size of myContainer is 210)
  cout \
    << "===============================================================================\n"\
    << " Trying to use [] with enumeration that is out of range: \n\n";
  try{
    myContainer[1000];
  }
  catch(std::logic_error err){
    cout << err.what() << "\n\n"; 
  }
  
  // Trying to create LexContainer using incompatible data array and indexRange. In this example
  // dataTeuchosArray corresponds to indexRange = {5,3,3,7} but we change the indexRange to one
  // that does not match the data. Note that if we permute the values in indexRange it will be
  // compatible with the data because it will specify the same size for the container. However, 
  // index bound permutation reshapes the container!
  cout \
    << "===============================================================================\n"\
    << " Trying to create LexContainer using incompatible data array and indexRange: \n\n";
  try{
    indexRange[0] = 5;
    indexRange[1] = 3;
    indexRange[2] = 3;
    indexRange[3] = 8;
    
    LexContainer<double> myOtherContainer(indexRange, dataTeuchosArray);
  }
  catch(std::logic_error err){
    cout << err.what() << endl;
  }
  
  return 0;
}
