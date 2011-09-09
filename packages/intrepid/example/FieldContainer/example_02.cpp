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

/** \file
\brief  Illustrates use of the FieldContainer class.
\author Created by P. Bochev and D. Ridzal
*/

#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                  Example use of the FieldContainer class                    |\n" \
  << "|                                                                             |\n" \
  << "|    1) using FieldContainer in DEBUG mode:                                   |\n" \
  << "|       requires intrepid to be configured with --enable-intrepid-debug       |\n" \
  << "|       See /test/FieldContainer/test_02.cpp for more examples                |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n\n";
  
  // Define variables to create and use FieldContainers
  Teuchos::Array<int> dimensions;
  Teuchos::Array<int> multiIndex;
  
  // Initialize dimensions for rank-4 multi-index value
  dimensions.resize(4);
  dimensions[0] = 5;
  dimensions[1] = 3;
  dimensions[2] = 2;
  dimensions[3] = 7;
  
  // Create a FieldContainer
  FieldContainer<double> myContainer(dimensions);
  
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
  
  // Trying to set values from array whose size is less than FieldContainer size
  cout \
    << "===============================================================================\n"\
    << " Trying to set values from array whose size is less than FieldContainer's size: \n\n";

  // Change one of the values of the dimensions to a lesser value: original value was 5 
  dimensions[0] = 4;
  
  // Define Teuchos::Array to store values with dimension equal to the number of multi-indexed values
  Teuchos::Array<double> dataTeuchosArray(4*3*2*7);
  
  // Fill with data
  int counter = 0;
  for(int i=0; i < dimensions[0]; i++){
    for(int j=0; j < dimensions[1]; j++){
      for(int k=0; k < dimensions[2]; k++){
        for(int l = 0; l < dimensions[3]; l++){
          dataTeuchosArray[counter] = (double)counter;
          counter++;
        }
      }
    }
  }
  
  // Now try to stuff this data into FieldContainer
  try{
    myContainer.setValues(dataTeuchosArray);
  }
  catch(std::logic_error err){
    cout << err.what() << "\n";
  }
  
  // Trying to set values from array whose size is greater than FieldContainer's size
  cout \
    << "===============================================================================\n"\
    << " Trying to set values from array whose size is greater than FieldContainer's size: \n\n";
  // Change one of the values of the dimensions to a lesser value: restore dimensions[0] to the 
  // value used to construct the LexArray and change dimensions[2] to a greater value
  dimensions[0] = 5;
  dimensions[2] = 3;
  
  // Define Teuchos::Array to store values with dimension equal to the number of multi-indexed values
  dataTeuchosArray.resize(5*3*3*7);
  
  // Fill with data
  counter = 0;
  for(int i=0; i < dimensions[0]; i++){
    for(int j=0; j < dimensions[1]; j++){
      for(int k=0; k < dimensions[2]; k++){
        for(int l = 0; l < dimensions[3]; l++){
          dataTeuchosArray[counter] = (double)counter;
          counter++;
        }
      }
    }
  }
  
  // Now try to stuff this data into FieldContainer
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
  
  // Trying to create FieldContainer using incompatible data array and dimensions. In this example
  // dataTeuchosArray corresponds to dimensions = {5,3,3,7} but we change the dimensions to one
  // that does not match the data. Note that if we permute the values in dimensions it will be
  // compatible with the data because it will specify the same size for the container. However, 
  // index bound permutation reshapes the container!
  cout \
    << "===============================================================================\n"\
    << " Trying to create FieldContainer using incompatible data array and dimensions: \n\n";
  try{
    dimensions[0] = 5;
    dimensions[1] = 3;
    dimensions[2] = 3;
    dimensions[3] = 8;
    
    FieldContainer<double> myOtherContainer(dimensions, dataTeuchosArray);
  }
  catch(std::logic_error err){
    cout << err.what() << endl;
  }
  
  return 0;
}
