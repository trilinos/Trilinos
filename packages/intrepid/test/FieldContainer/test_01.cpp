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
\brief  Unit test of FieldContainer class
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_FieldContainer.hpp"
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
    << "|                           Unit Test FieldContainer                          |\n" \
    << "|                                                                             |\n" \
    << "|     1) Value accesss by multi-index and enumeration, setting values         |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  
  int errorFlag  = 0;
  
  // Define variables to construct and use FieldContainer:
  Teuchos::Array<double> dataArray;
  Teuchos::Array<int> dimension;
  Teuchos::Array<int> multiIndex;
  int dim0, dim1, dim2, dim3, dim4;
  int containerSize;
  int containerRank;
  int dataSize;
  int counter;
  

  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 1: Rank-1 FieldContainer                                               |\n"\
    << "===============================================================================\n";
  
  // Adjust for rank-1 containers
  dimension.resize(1);
  multiIndex.resize(1);
  dim0 = 267; 
  
  // Set dimensions
  dimension[0] = dim0;
  
  // Define Teuchos::Array with dimension equal to FieldContainer capacity (determined from dimension)
  dataSize = dimension[0];
  dataArray.resize(dataSize);
  
  // Fill with data
  counter = 0;
  for(int i=0; i < dimension[0]; i++){
    dataArray[counter] = (double)counter;
    counter++;
  }
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 1-a: Compare constructors with array and list of dimensions            |\n"\
    << "===============================================================================\n";
  try{
    
    // Using ctor with array of dimensions and array of data
    FieldContainer<double> rank1Container(dimension, dataArray);
    
    // Using ctor with list of dimensions and no data
    FieldContainer<double> rank1ContainerAlt(dim0);
    
    // Initialize variables
    containerSize = rank1Container.getSize();
    containerRank = rank1Container.getRank();
    multiIndex.resize(containerRank);
    
    if( rank1Container.getSize() != rank1ContainerAlt.getSize() ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Size using ctor with array of dimensions = " << rank1Container.getSize() << "\n";
      *outStream << " Size using ctor with list of dimensions  = " << rank1ContainerAlt.getSize() << "\n";
    }
    
    if( rank1Container.getRank() != rank1ContainerAlt.getRank() ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Rank using ctor with array of dimensions = " << rank1Container.getRank() << "\n";
      *outStream << " Rank using ctor with list of dimensions  = " << rank1ContainerAlt.getRank() << "\n";
    }
    
    for(int dim = 0; dim < rank1Container.getRank(); dim ++ ) {
      if( rank1Container.getDimension(dim) != rank1ContainerAlt.getDimension(dim) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Dimension " << dim << " using ctor with array of dimensions = " << rank1Container.getDimension(dim) << "\n";
        *outStream << " Dimension " << dim << " using ctor with list of dimensions  = " << rank1ContainerAlt.getDimension(dim) << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 1-b: Access by enumeration, multi-index array and multi-index list     |\n"\
      << "===============================================================================\n";
    
    // Loop over container by enumeration
    for(int enumeration = 0; enumeration < containerSize; enumeration++) {
      int i0;
      
      // Convert enumeration to multi-index array and multi-index list and compare values
      rank1Container.getMultiIndex(multiIndex, enumeration);
      rank1Container.getMultiIndex(i0, enumeration);
      if( (multiIndex[0] != i0) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Multi-index as array = [" 
          << multiIndex[0] << "]\n";
        *outStream << " Multi-index as list  = (" << i0 << ")\n";
      }      
      
      // Check if access by enumeration gives the same value as access by multi-index array
      if( rank1Container[enumeration] != rank1Container.getValue(multiIndex) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by enumeration       = " << rank1Container[enumeration] << "\n";
        *outStream << " Value by multi-index array = " << rank1Container.getValue(multiIndex) << "\n";
      }
      
      // Check if access by multi-index list gives the same value as access by multi-index array
      if( rank1Container(i0) != rank1Container.getValue(multiIndex) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by multi-index list  = " << rank1Container(i0) << "\n";
        *outStream << " Value by multi-index array = " << rank1Container.getValue(multiIndex) << "\n";
      }
      
      // Check if access by multi-index list gives the same value as access by [] (only for rank-1!)
      if( rank1Container(i0) != rank1Container[i0] ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by multi-index list  = " << rank1Container(i0) << "\n";
        *outStream << " Value by overloaded []     = " << rank1Container[i0] << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 1-c: Access by multi-index array and list & compare with data array    |\n"\
      << "===============================================================================\n";
    
    // Loop over container by multi-index
    for(int i=0; i < dimension[0]; i++){
      multiIndex[0] = i;
      
      // Method that takes array of multi-indices
      int enumeration = rank1Container.getEnumeration(multiIndex);
      
      // Compare with method that takes a list of multi-indices
      if( enumeration != rank1Container.getEnumeration(i) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Enumeration from multi-index array  = " << enumeration << "\n";
        *outStream << " Enumeration from multi-index list   = " << rank1Container.getEnumeration(i) << "\n";
      }
      
      // Check if access by multi-index array matches values in original dataArray 
      if(dataArray[enumeration] != rank1Container.getValue(multiIndex)) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by multi-index array = " << rank1Container.getValue(multiIndex) << "\n";
        *outStream << " Value from data array      = " << dataArray[enumeration] << "\n";
      }
      
      // Check if access by multi-index list matches values in original dataArray
      if(dataArray[enumeration] != rank1Container(i)) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by multi-index list = " << rank1Container(i) << "\n";
        *outStream << " Value from data array     = " << dataArray[enumeration] << "\n";
      }
      
      // Check if access by multi-index array matches access by multi-index list
      if( rank1Container(i) != rank1Container.getValue(multiIndex)) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by multi-index array = " << rank1Container.getValue(multiIndex) << "\n";
        *outStream << " Value by multi-index list  = " << rank1Container(i) << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 1-d: Store zeroes and empty container                                  |\n"\
      << "===============================================================================\n";
    
    rank1Container.storeZero();
    double sum = 0.0;
    for (int i=0; i<rank1Container.getSize(); i++) {
      sum += rank1Container[i];
    }
    if( (sum != 0.0) ) {
      errorFlag++;
      *outStream << " Container size = " << rank1Container.getSize() << "\n";
      *outStream << " Container rank = " << rank1Container.getRank() << "\n";      
    }
    
    rank1Container.clear();
    if( !(rank1Container.getSize() == 0 && rank1Container.getRank() == 0)) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Container size = " << rank1Container.getSize() << "\n";
      *outStream << " Container rank = " << rank1Container.getRank() << "\n";
    }
  } //try 
  
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  };
  
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 2: Rank-2 FieldContainer                                               |\n"\
    << "===============================================================================\n";
  
  // Adjust for rank-2 containers
  dimension.resize(2);
  multiIndex.resize(2);
  dim0 = 1; 
  dim1 = 55; 
  
  // Set dimensions
  dimension[0] = dim0;
  dimension[1] = dim1;
  
  // Define Teuchos::Array with dimension equal to FieldContainer capacity (determined from dimension)
  dataSize = dimension[0];
  for(int r = 1; r < (int)dimension.size(); r++){
    dataSize *= dimension[r];
  }
  dataArray.resize(dataSize);
  
  // Fill with data
  counter = 0;
  for(int i=0; i < dimension[0]; i++){
    for(int j=0; j < dimension[1]; j++){
      dataArray[counter] = (double)counter;
      counter++;
    }
  }
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 2-a: Compare constructors with array and list of dimensions            |\n"\
    << "===============================================================================\n";
  try{
    
    // Using ctor with array of dimensions and array of data
    FieldContainer<double> rank2Container(dimension, dataArray);
    
    // Using ctor with list of dimensions and no data
    FieldContainer<double> rank2ContainerAlt(dim0, dim1);
    
    // Initialize variables
    containerSize = rank2Container.getSize();
    containerRank = rank2Container.getRank();
    multiIndex.resize(containerRank);
    
    if( rank2Container.getSize() != rank2ContainerAlt.getSize() ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Size using ctor with array of dimensions = " << rank2Container.getSize() << "\n";
      *outStream << " Size using ctor with list of dimensions  = " << rank2ContainerAlt.getSize() << "\n";
    }
    
    if( rank2Container.getRank() != rank2ContainerAlt.getRank() ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Rank using ctor with array of dimensions = " << rank2Container.getRank() << "\n";
      *outStream << " Rank using ctor with list of dimensions  = " << rank2ContainerAlt.getRank() << "\n";
    }
    
    for(int dim = 0; dim < rank2Container.getRank(); dim ++ ) {
      if( rank2Container.getDimension(dim) != rank2ContainerAlt.getDimension(dim) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Dimension " << dim << " using ctor with array of dimensions = " << rank2Container.getDimension(dim) << "\n";
        *outStream << " Dimension " << dim << " using ctor with list of dimensions  = " << rank2ContainerAlt.getDimension(dim) << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 2-b: Access by enumeration, multi-index array and multi-index list     |\n"\
      << "===============================================================================\n";
    
    // Loop over container by enumeration
    for(int enumeration = 0; enumeration < containerSize; enumeration++) {
      int i0,i1;
      
      // Convert enumeration to multi-index array and multi-index list and compare values
      rank2Container.getMultiIndex(multiIndex, enumeration);
      rank2Container.getMultiIndex(i0, i1, enumeration);
      if( (multiIndex[0] != i0) ||
          (multiIndex[1] != i1) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Multi-index as array = [" 
          << multiIndex[0] << multiIndex[1] << "]\n";
        *outStream << " Multi-index as list  = (" << i0 << "," << i1  << ")\n";
      }      
      
      // Check if access by enumeration gives the same value as access by multi-index array
      if( rank2Container[enumeration] != rank2Container.getValue(multiIndex) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by enumeration       = " << rank2Container[enumeration] << "\n";
        *outStream << " Value by multi-index array = " << rank2Container.getValue(multiIndex) << "\n";
      }
      
      // Check if access by multi-index list gives the same value as access by multi-index array
      if( rank2Container(i0,i1) != rank2Container.getValue(multiIndex) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by multi-index list  = " << rank2Container(i0,i1) << "\n";
        *outStream << " Value by multi-index array = " << rank2Container.getValue(multiIndex) << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 2-c: Access by multi-index array and list & compare with data array    |\n"\
      << "===============================================================================\n";
    
    // Loop over container by multi-index
    for(int i=0; i < dimension[0]; i++){
      multiIndex[0] = i;
      for(int j=0; j < dimension[1]; j++){
        multiIndex[1] = j;
        
        // Method that takes array of multi-indices
        int enumeration = rank2Container.getEnumeration(multiIndex);
        
        // Compare with method that takes a list of multi-indices
        if( enumeration != rank2Container.getEnumeration(i,j) ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Enumeration from multi-index array  = " << enumeration << "\n";
          *outStream << " Enumeration from multi-index list   = " << rank2Container.getEnumeration(i,j) << "\n";
        }
        
        // Check if access by multi-index array matches values in original dataArray 
        if(dataArray[enumeration] != rank2Container.getValue(multiIndex)) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Value by multi-index array = " << rank2Container.getValue(multiIndex) << "\n";
          *outStream << " Value from data array      = " << dataArray[enumeration] << "\n";
        }
        
        // Check if access by multi-index list matches values in original dataArray
        if(dataArray[enumeration] != rank2Container(i,j)) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Value by multi-index list = " << rank2Container(i,j) << "\n";
          *outStream << " Value from data array     = " << dataArray[enumeration] << "\n";
        }
        
        // Check if access by multi-index array matches access by multi-index list
        if( rank2Container(i,j) != rank2Container.getValue(multiIndex)) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Value by multi-index array = " << rank2Container.getValue(multiIndex) << "\n";
          *outStream << " Value by multi-index list  = " << rank2Container(i,j) << "\n";
        }
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 2-d: Store zeroes and empty container                                  |\n"\
      << "===============================================================================\n";
    
    rank2Container.storeZero();
    double sum = 0.0;
    for (int i=0; i<rank2Container.getSize(); i++) {
      sum += rank2Container[i];
    }
    if( (sum != 0.0) ) {
      errorFlag++;
      *outStream << " Container size = " << rank2Container.getSize() << "\n";
      *outStream << " Container rank = " << rank2Container.getRank() << "\n";      
    }
    
    rank2Container.clear();
    if( !(rank2Container.getSize() == 0 && rank2Container.getRank() == 0)) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Container size = " << rank2Container.getSize() << "\n";
      *outStream << " Container rank = " << rank2Container.getRank() << "\n";
    }
  } //try 
  
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  };
  
  
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 3: Rank-3 FieldContainer                                               |\n"\
    << "===============================================================================\n";
  
  // Adjust for rank-3 containers
  dimension.resize(3);
  multiIndex.resize(3);
  dim0 = 17; 
  dim1 = 55; 
  dim2 = 10; 
  
  // Set dimensions
  dimension[0] = dim0;
  dimension[1] = dim1;
  dimension[2] = dim2;
  
  // Define Teuchos::Array with dimension equal to FieldContainer capacity (determined from dimension)
  dataSize = dimension[0];
  for(int r = 1; r < (int)dimension.size(); r++){
    dataSize *= dimension[r];
  }
  dataArray.resize(dataSize);
  
  // Fill with data
  counter = 0;
  for(int i=0; i < dimension[0]; i++){
    for(int j=0; j < dimension[1]; j++){
      for(int k=0; k < dimension[2]; k++){
        dataArray[counter] = (double)counter;
        counter++;
      }
    }
  }
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 3-a: Compare constructors with array and list of dimensions            |\n"\
    << "===============================================================================\n";
  try{
    
    // Using ctor with array of dimensions and array of data
    FieldContainer<double> rank3Container(dimension, dataArray);
    
    // Using ctor with list of dimensions and no data
    FieldContainer<double> rank3ContainerAlt(dim0, dim1, dim2);
    
    // Initialize variables
    containerSize = rank3Container.getSize();
    containerRank = rank3Container.getRank();
    multiIndex.resize(containerRank);
    
    if( rank3Container.getSize() != rank3ContainerAlt.getSize() ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Size using ctor with array of dimensions = " << rank3Container.getSize() << "\n";
      *outStream << " Size using ctor with list of dimensions  = " << rank3ContainerAlt.getSize() << "\n";
    }
    
    if( rank3Container.getRank() != rank3ContainerAlt.getRank() ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Rank using ctor with array of dimensions = " << rank3Container.getRank() << "\n";
      *outStream << " Rank using ctor with list of dimensions  = " << rank3ContainerAlt.getRank() << "\n";
    }
    
    for(int dim = 0; dim < rank3Container.getRank(); dim ++ ) {
      if( rank3Container.getDimension(dim) != rank3ContainerAlt.getDimension(dim) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Dimension " << dim << " using ctor with array of dimensions = " << rank3Container.getDimension(dim) << "\n";
        *outStream << " Dimension " << dim << " using ctor with list of dimensions  = " << rank3ContainerAlt.getDimension(dim) << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 3-b: Access by enumeration, multi-index array and multi-index list     |\n"\
      << "===============================================================================\n";
    
    // Loop over container by enumeration
    for(int enumeration = 0; enumeration < containerSize; enumeration++) {
      int i0,i1,i2;
      
      // Convert enumeration to multi-index array and multi-index list and compare values
      rank3Container.getMultiIndex(multiIndex, enumeration);
      rank3Container.getMultiIndex(i0, i1, i2, enumeration);
      if( (multiIndex[0] != i0) ||
          (multiIndex[1] != i1) ||
          (multiIndex[2] != i2) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Multi-index as array = [" 
          << multiIndex[0] << multiIndex[1] << multiIndex[2]  << "]\n";
        *outStream << " Multi-index as list  = (" << i0 << "," << i1 << "," << i2 << ")\n";
      }      
      
      // Check if access by enumeration gives the same value as access by multi-index array
      if( rank3Container[enumeration] != rank3Container.getValue(multiIndex) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by enumeration       = " << rank3Container[enumeration] << "\n";
        *outStream << " Value by multi-index array = " << rank3Container.getValue(multiIndex) << "\n";
      }
      
      // Check if access by multi-index list gives the same value as access by multi-index array
      if( rank3Container(i0,i1,i2) != rank3Container.getValue(multiIndex) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by multi-index list  = " << rank3Container(i0,i1,i2) << "\n";
        *outStream << " Value by multi-index array = " << rank3Container.getValue(multiIndex) << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 3-c: Access by multi-index array and list & compare with data array    |\n"\
      << "===============================================================================\n";
    
    // Loop over container by multi-index
    for(int i=0; i < dimension[0]; i++){
      multiIndex[0] = i;
      for(int j=0; j < dimension[1]; j++){
        multiIndex[1] = j;
        for(int k=0; k < dimension[2]; k++){
          multiIndex[2] = k; 
          
          // Method that takes array of multi-indices
          int enumeration = rank3Container.getEnumeration(multiIndex);
          
          // Compare with method that takes a list of multi-indices
          if( enumeration != rank3Container.getEnumeration(i,j,k) ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Enumeration from multi-index array  = " << enumeration << "\n";
            *outStream << " Enumeration from multi-index list   = " << rank3Container.getEnumeration(i,j,k) << "\n";
          }
          
          // Check if access by multi-index array matches values in original dataArray 
          if(dataArray[enumeration] != rank3Container.getValue(multiIndex)) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Value by multi-index array = " << rank3Container.getValue(multiIndex) << "\n";
            *outStream << " Value from data array      = " << dataArray[enumeration] << "\n";
          }
          
          // Check if access by multi-index list matches values in original dataArray
          if(dataArray[enumeration] != rank3Container(i,j,k)) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Value by multi-index list = " << rank3Container(i,j,k) << "\n";
            *outStream << " Value from data array     = " << dataArray[enumeration] << "\n";
          }
          
          // Check if access by multi-index array matches access by multi-index list
          if( rank3Container(i,j,k) != rank3Container.getValue(multiIndex)) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Value by multi-index array = " << rank3Container.getValue(multiIndex) << "\n";
            *outStream << " Value by multi-index list  = " << rank3Container(i,j,k) << "\n";
          }
        }
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 3-d: Store zeroes and empty container                                  |\n"\
      << "===============================================================================\n";
    
    rank3Container.storeZero();
    double sum = 0.0;
    for (int i=0; i<rank3Container.getSize(); i++) {
      sum += rank3Container[i];
    }
    if( (sum != 0.0) ) {
      errorFlag++;
      *outStream << " Container size = " << rank3Container.getSize() << "\n";
      *outStream << " Container rank = " << rank3Container.getRank() << "\n";      
    }
    
    rank3Container.clear();
    if( !(rank3Container.getSize() == 0 && rank3Container.getRank() == 0)) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Container size = " << rank3Container.getSize() << "\n";
      *outStream << " Container rank = " << rank3Container.getRank() << "\n";
    }
  } //try 
  
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  };
  

  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 4: Rank-4 FieldContainer                                               |\n"\
    << "===============================================================================\n";
  
  // Adjust for rank-4 containers
  dimension.resize(4);
  multiIndex.resize(4);
  dim0 = 27; 
  dim1 = 4; 
  dim2 = 11; 
  dim3 = 6; 
    
  // Set dimensions
  dimension[0] = dim0;
  dimension[1] = dim1;
  dimension[2] = dim2;
  dimension[3] = dim3;
  
  // Define Teuchos::Array with dimension equal to FieldContainer capacity (determined from dimension)
  dataSize = dimension[0];
  for(int r = 1; r < (int)dimension.size(); r++){
    dataSize *= dimension[r];
  }
  dataArray.resize(dataSize);
  
  // Fill with data
  counter = 0;
  for(int i=0; i < dimension[0]; i++){
    for(int j=0; j < dimension[1]; j++){
      for(int k=0; k < dimension[2]; k++){
        for(int l = 0; l < dimension[3]; l++){
          dataArray[counter] = (double)counter;
          counter++;
        }
      }
    }
  }
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 4-a: Compare constructors with array and list of dimensions            |\n"\
    << "===============================================================================\n";
  try{
    
    // Using ctor with array of dimensions and array of data
    FieldContainer<double> rank4Container(dimension, dataArray);
    
    // Using ctor with list of dimensions and no data
    FieldContainer<double> rank4ContainerAlt(dim0, dim1, dim2, dim3);
    
    // Initialize variables
    containerSize = rank4Container.getSize();
    containerRank = rank4Container.getRank();
    multiIndex.resize(containerRank);
    
    if( rank4Container.getSize() != rank4ContainerAlt.getSize() ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Size using ctor with array of dimensions = " << rank4Container.getSize() << "\n";
      *outStream << " Size using ctor with list of dimensions  = " << rank4ContainerAlt.getSize() << "\n";
    }
    
    if( rank4Container.getRank() != rank4ContainerAlt.getRank() ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Rank using ctor with array of dimensions = " << rank4Container.getRank() << "\n";
      *outStream << " Rank using ctor with list of dimensions  = " << rank4ContainerAlt.getRank() << "\n";
    }
    
    for(int dim = 0; dim < rank4Container.getRank(); dim ++ ) {
      if( rank4Container.getDimension(dim) != rank4ContainerAlt.getDimension(dim) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Dimension " << dim << " using ctor with array of dimensions = " << rank4Container.getDimension(dim) << "\n";
        *outStream << " Dimension " << dim << " using ctor with list of dimensions  = " << rank4ContainerAlt.getDimension(dim) << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 4-b: Access by enumeration, multi-index array and multi-index list     |\n"\
      << "===============================================================================\n";
    
    // Loop over container by enumeration
    for(int enumeration = 0; enumeration < containerSize; enumeration++) {
      int i0,i1,i2,i3;
      
      // Convert enumeration to multi-index array and multi-index list and compare values
      rank4Container.getMultiIndex(multiIndex, enumeration);
      rank4Container.getMultiIndex(i0, i1, i2, i3, enumeration);
      if( (multiIndex[0] != i0) ||
          (multiIndex[1] != i1) ||
          (multiIndex[2] != i2) ||
          (multiIndex[3] != i3) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Multi-index as array = [" 
          << multiIndex[0] << multiIndex[1] << multiIndex[2] << multiIndex[3]  << "]\n";
        *outStream << " Multi-index as list  = (" << i0 << "," << i1 << "," << i2 << "," << i3  << ")\n";
      }      
      
      // Check if access by enumeration gives the same value as access by multi-index array
      if( rank4Container[enumeration] != rank4Container.getValue(multiIndex) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by enumeration       = " << rank4Container[enumeration] << "\n";
        *outStream << " Value by multi-index array = " << rank4Container.getValue(multiIndex) << "\n";
      }
      
      // Check if access by multi-index list gives the same value as access by multi-index array
      if( rank4Container(i0,i1,i2,i3) != rank4Container.getValue(multiIndex) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by multi-index list  = " << rank4Container(i0,i1,i2,i3) << "\n";
        *outStream << " Value by multi-index array = " << rank4Container.getValue(multiIndex) << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 4-c: Access by multi-index array and list & compare with data array    |\n"\
      << "===============================================================================\n";
    
    // Loop over container by multi-index
    for(int i=0; i < dimension[0]; i++){
      multiIndex[0] = i;
      for(int j=0; j < dimension[1]; j++){
        multiIndex[1] = j;
        for(int k=0; k < dimension[2]; k++){
          multiIndex[2] = k; 
          for(int l = 0; l < dimension[3]; l++){
            multiIndex[3] = l;
            
            // Method that takes array of multi-indices
            int enumeration = rank4Container.getEnumeration(multiIndex);
            
            // Compare with method that takes a list of multi-indices
            if( enumeration != rank4Container.getEnumeration(i,j,k,l) ) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << " Enumeration from multi-index array  = " << enumeration << "\n";
              *outStream << " Enumeration from multi-index list   = " << rank4Container.getEnumeration(i,j,k,l) << "\n";
            }
            
            // Check if access by multi-index array matches values in original dataArray 
            if(dataArray[enumeration] != rank4Container.getValue(multiIndex)) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << " Value by multi-index array = " << rank4Container.getValue(multiIndex) << "\n";
              *outStream << " Value from data array      = " << dataArray[enumeration] << "\n";
            }
            
            // Check if access by multi-index list matches values in original dataArray
            if(dataArray[enumeration] != rank4Container(i,j,k,l)) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << " Value by multi-index list = " << rank4Container(i,j,k,l) << "\n";
              *outStream << " Value from data array     = " << dataArray[enumeration] << "\n";
            }
            
            // Check if access by multi-index array matches access by multi-index list
            if( rank4Container(i,j,k,l) != rank4Container.getValue(multiIndex)) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << " Value by multi-index array = " << rank4Container.getValue(multiIndex) << "\n";
              *outStream << " Value by multi-index list  = " << rank4Container(i,j,k,l) << "\n";
            }
          }
        }
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 4-d: Store zeroes and empty container                                  |\n"\
      << "===============================================================================\n";
    
    rank4Container.storeZero();
    double sum = 0.0;
    for (int i=0; i<rank4Container.getSize(); i++) {
      sum += rank4Container[i];
    }
    if( (sum != 0.0) ) {
      errorFlag++;
      *outStream << " Container size = " << rank4Container.getSize() << "\n";
      *outStream << " Container rank = " << rank4Container.getRank() << "\n";      
    }
    
    rank4Container.clear();
    if( !(rank4Container.getSize() == 0 && rank4Container.getRank() == 0)) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Container size = " << rank4Container.getSize() << "\n";
      *outStream << " Container rank = " << rank4Container.getRank() << "\n";
    }
  } //try 
  
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  };
  
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 5: Rank-5 FieldContainer                                               |\n"\
    << "===============================================================================\n";
  
  // Adjust for rank-5 containers
  dimension.resize(5);
  multiIndex.resize(5);
  dim0 = 3; 
  dim1 = 7; 
  dim2 = 11; 
  dim3 = 5; 
  dim4 = 6;
  
  // Set dimensions
  dimension[0] = dim0;
  dimension[1] = dim1;
  dimension[2] = dim2;
  dimension[3] = dim3;
  dimension[4] = dim4;
  
  // Define Teuchos::Array with dimension equal to FieldContainer capacity (determined from dimension)
  dataSize = dimension[0];
  for(int r = 1; r < (int)dimension.size(); r++){
    dataSize *= dimension[r];
  }
  dataArray.resize(dataSize);
  
  // Fill with data
  counter = 0;
  for(int i=0; i < dimension[0]; i++){
    for(int j=0; j < dimension[1]; j++){
      for(int k=0; k < dimension[2]; k++){
        for(int l = 0; l < dimension[3]; l++){
          for(int m = 0; m < dimension[4]; m++){
            dataArray[counter] = (double)counter;
            counter++;
          }
        }
      }
    }
  }
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 5-a: Compare constructors with array and list of dimensions            |\n"\
    << "===============================================================================\n";
  try{
    
    // Using ctor with array of dimensions and array of data
    FieldContainer<double> rank5Container(dimension, dataArray);
    
    // Using ctor with list of dimensions and no data
    FieldContainer<double> rank5ContainerAlt(dim0, dim1, dim2, dim3, dim4);
    
    // Initialize variables
    containerSize = rank5Container.getSize();
    containerRank = rank5Container.getRank();
    multiIndex.resize(containerRank);
    
    if( rank5Container.getSize() != rank5ContainerAlt.getSize() ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Size using ctor with array of dimensions = " << rank5Container.getSize() << "\n";
      *outStream << " Size using ctor with list of dimensions  = " << rank5ContainerAlt.getSize() << "\n";
    }
    
    if( rank5Container.getRank() != rank5ContainerAlt.getRank() ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Rank using ctor with array of dimensions = " << rank5Container.getRank() << "\n";
      *outStream << " Rank using ctor with list of dimensions  = " << rank5ContainerAlt.getRank() << "\n";
    }
    
    for(int dim = 0; dim < rank5Container.getRank(); dim ++ ) {
      if( rank5Container.getDimension(dim) != rank5ContainerAlt.getDimension(dim) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Dimension " << dim << " using ctor with array of dimensions = " << rank5Container.getDimension(dim) << "\n";
        *outStream << " Dimension " << dim << " using ctor with list of dimensions  = " << rank5ContainerAlt.getDimension(dim) << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 5-b: Access by enumeration, multi-index array and multi-index list     |\n"\
      << "===============================================================================\n";
    
    // Loop over container by enumeration
    for(int enumeration = 0; enumeration < containerSize; enumeration++) {
      int i0,i1,i2,i3,i4;
      
      // Convert enumeration to multi-index array and multi-index list and compare values
      rank5Container.getMultiIndex(multiIndex, enumeration);
      rank5Container.getMultiIndex(i0, i1, i2, i3, i4, enumeration);
      if( (multiIndex[0] != i0) ||
          (multiIndex[1] != i1) ||
          (multiIndex[2] != i2) ||
          (multiIndex[3] != i3) ||
          (multiIndex[4] != i4) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Multi-index as array = [" 
          << multiIndex[0] << multiIndex[1] << multiIndex[2] << multiIndex[3] << multiIndex[4] << "]\n";
        *outStream << " Multi-index as list  = (" << i0 << "," << i1 << "," << i2 << "," << i3 << "," << i4 << ")\n";
      }      
      
      // Check if access by enumeration gives the same value as access by multi-index array
      if( rank5Container[enumeration] != rank5Container.getValue(multiIndex) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by enumeration       = " << rank5Container[enumeration] << "\n";
        *outStream << " Value by multi-index array = " << rank5Container.getValue(multiIndex) << "\n";
      }
      
      // Check if access by multi-index list gives the same value as access by multi-index array
      if( rank5Container(i0,i1,i2,i3,i4) != rank5Container.getValue(multiIndex) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Value by multi-index list  = " << rank5Container(i0,i1,i2,i3,i4) << "\n";
        *outStream << " Value by multi-index array = " << rank5Container.getValue(multiIndex) << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 5-c: Access by multi-index array and list & compare with data array    |\n"\
      << "===============================================================================\n";
    
    // Loop over container by multi-index
    for(int i=0; i < dimension[0]; i++){
      multiIndex[0] = i;
      for(int j=0; j < dimension[1]; j++){
        multiIndex[1] = j;
        for(int k=0; k < dimension[2]; k++){
          multiIndex[2] = k; 
          for(int l = 0; l < dimension[3]; l++){
            multiIndex[3] = l;
            for(int m = 0; m < dimension[4]; m++){
              multiIndex[4] = m;
              
              // Method that takes array of multi-indices
              int enumeration = rank5Container.getEnumeration(multiIndex);
              
              // Compare with method that takes a list of multi-indices
              if( enumeration != rank5Container.getEnumeration(i,j,k,l,m) ) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << " Enumeration from multi-index array  = " << enumeration << "\n";
                *outStream << " Enumeration from multi-index list   = " << rank5Container.getEnumeration(i,j,k,l,m) << "\n";
              }
              
              // Check if access by multi-index array matches values in original dataArray 
              if(dataArray[enumeration] != rank5Container.getValue(multiIndex)) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << " Value by multi-index array = " << rank5Container.getValue(multiIndex) << "\n";
                *outStream << " Value from data array      = " << dataArray[enumeration] << "\n";
              }
              
              // Check if access by multi-index list matches values in original dataArray
              if(dataArray[enumeration] != rank5Container(i,j,k,l,m)) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << " Value by multi-index list = " << rank5Container(i,j,k,l,m) << "\n";
                *outStream << " Value from data array     = " << dataArray[enumeration] << "\n";
              }
              
              // Check if access by multi-index array matches access by multi-index list
              if( rank5Container(i,j,k,l,m) != rank5Container.getValue(multiIndex)) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << " Value by multi-index array = " << rank5Container.getValue(multiIndex) << "\n";
                *outStream << " Value by multi-index list  = " << rank5Container(i,j,k,l,m) << "\n";
              }
            }
          }
        }
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 5-d: Store zeroes and empty container                                  |\n"\
      << "===============================================================================\n";
    
    rank5Container.storeZero();
    double sum = 0.0;
    for (int i=0; i<rank5Container.getSize(); i++) {
      sum += rank5Container[i];
    }
    if( (sum != 0.0) ) {
      errorFlag++;
      *outStream << " Container size = " << rank5Container.getSize() << "\n";
      *outStream << " Container rank = " << rank5Container.getRank() << "\n";      
    }
    
    rank5Container.clear();
    if( !(rank5Container.getSize() == 0 && rank5Container.getRank() == 0)) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Container size = " << rank5Container.getSize() << "\n";
      *outStream << " Container rank = " << rank5Container.getRank() << "\n";
    }
  } //try 
  
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  };
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 6: Resize container based on another container                         |\n"\
    << "===============================================================================\n";
  
  try{
    FieldContainer<double> myContainer(1,2,3);
    FieldContainer<double> hisContainer(5,4,3,2);
    
    hisContainer.resize(myContainer);
    if ( hisContainer.getRank() != myContainer.getRank() ) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " Rank of target container   = " << hisContainer.getRank() << "\n";
      *outStream << " Rank of argument container = " << myContainer.getRank() << "\n";
    }
    
    for(int dim = 0; dim < myContainer.getRank(); dim++){
      if ( hisContainer.getDimension(dim) != myContainer.getDimension(dim) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Dimension " << dim << " of target container   = " << hisContainer.getDimension(dim) << "\n";
        *outStream << " Dimension " << dim << " of argument container = " << myContainer.getDimension(dim) << "\n";
      }
    }
  }// try
    
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
