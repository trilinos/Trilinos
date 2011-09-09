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
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  std::cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                  Example use of the FieldContainer class                    |\n" \
  << "|                                                                             |\n" \
  << "|    1) Creating and filling FieldContainer objects                           |\n" \
  << "|    2) Accessing elements in FieldContainer objects                          |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n\n";
  
  // Define variables to create and use FieldContainers
  Teuchos::Array<int> dimension;
  Teuchos::Array<int> multiIndex;
  
  std::cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 1: rank 2 multi-index: {u(p,i) | 0 <= p < 5; 0 <= i < 3 }           |\n"\
    << "===============================================================================\n\n";
  // This rank 2 multi-indexed value can be used to store values of 3D vector field
  // evaluated at 5 points, or values of gradient of scalar field evaluated at the points
  
  // Resize dimension and multiIndex for rank-2 multi-indexed value 
  dimension.resize(2);
  multiIndex.resize(2);
  
  // Load upper ranges for the two indices in the multi-indexed value
  dimension[0] = 5;
  dimension[1] = 3;
    
  // Create FieldContainer that can hold the rank-2 multi-indexed value
  FieldContainer<double> myContainer(dimension);
  
  // Fill with some data: leftmost index changes last, rightmost index changes first!  
  for(int p = 0; p < dimension[0]; p++){
    multiIndex[0] = p;
    
    for(int i = 0; i < dimension[1]; i++){
      multiIndex[1] = i;
      
      // Load value with multi-index {p,i} to container
      myContainer.setValue((double)(i+p), multiIndex);
    }
  }
  
  // Show container contents
  std::cout << myContainer;
  
  // Access by overloaded (), multiindex and []:
  multiIndex[0] = 3; 
  multiIndex[1] = 1;
  int enumeration = myContainer.getEnumeration(multiIndex);

  std::cout << "Access by ():          myContainer(" << 3 <<"," << 1 << ") = " << myContainer(3,1) << "\n";  
  std::cout << "Access by multi-index: myContainer{" << multiIndex[0] << multiIndex[1] << "} = "<< myContainer.getValue(multiIndex) <<"\n";
  std::cout << "Access by enumeration: myContainer[" << enumeration << "] = " << myContainer[enumeration] <<"\n";
  
  std::cout << "Assigning value by (): \n old value at (3,1) = " << myContainer(3,1) <<"\n";
  myContainer(3,1) = 999.99;
  std::cout << " new value at (3,1) = " << myContainer(3,1) <<"\n";
      
  
  std::cout << "\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 2: rank 3 multi-index: {u(p,i,j) | 0 <=p< 5; 0 <= i<2, 0<=j<3       |\n"\
  << "===============================================================================\n\n";
  // This rank-3 value can be used to store subset of second partial derivatives  values
  // of a scalar function at p points.
  
  // Resize dimension and multiIndex for rank-3 multi-indexed value 
  dimension.resize(3);
  multiIndex.resize(3);
  
  // Define upper ranges for the three indices in the multi-indexed value
  dimension[0] = 5;
  dimension[1] = 2;
  dimension[2] = 3;
  
  // Reset the existing container to accept rank-3 value with the specified index ranges
  myContainer.resize(dimension);
  
  // Fill with some data
  for(int p = 0; p < dimension[0]; p++){
    multiIndex[0] = p;
    for(int i = 0; i < dimension[1]; i++){
      multiIndex[1] = i;
      for(int j = 0; j < dimension[2]; j++){
        multiIndex[2] = j;
        
        // Load value with multi-index {p,i} to container
        myContainer.setValue((double)(p+i+j), multiIndex);
      }
    }
  }
  
  // Display contents
  std::cout << myContainer;
  
  // Access by overloaded (), multiindex and []:
  multiIndex[0] = 3; 
  multiIndex[1] = 1;
  multiIndex[2] = 2;
  enumeration = myContainer.getEnumeration(multiIndex);
    
  std::cout << "Access by ():          myContainer(" << 3 <<"," << 1 << "," << 2 << ") = " << myContainer(3,1,2) << "\n";  
  std::cout << "Access by multi-index: myContainer{" << multiIndex[0] << multiIndex[1] << multiIndex[2] << "} = "<< myContainer.getValue(multiIndex) <<"\n";
  std::cout << "Access by enumeration: myContainer[" << enumeration << "] = " << myContainer[enumeration] <<"\n";

  std::cout << "Assigning value by (): \n old value at (3,1,2) = " << myContainer(3,1,2) <<"\n";
  myContainer(3,1,2) = -999.999;
  std::cout << " new value at (3,1,2) = " << myContainer(3,1,2) <<"\n";
  
  
  std::cout << "\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 4: making rank-5 FieldContainer from data array and index range array |\n"\
  << "===============================================================================\n\n";
  // Initialize dimension for rank-5 multi-index value
  dimension.resize(5);
  dimension[0] = 5;
  dimension[1] = 2;
  dimension[2] = 3;
  dimension[3] = 4;
  dimension[4] = 6;
  
  // Define Teuchos::Array to store values with dimension equal to the number of multi-indexed values
  Teuchos::Array<double> dataTeuchosArray(5*2*3*4*6);
  
  // Fill with data
  int counter = 0;
  for(int i=0; i < dimension[0]; i++){
    for(int j=0; j < dimension[1]; j++){
      for(int k=0; k < dimension[2]; k++){
        for(int l = 0; l < dimension[3]; l++){
          for(int m = 0; m < dimension[4]; m++){
            dataTeuchosArray[counter] = (double)counter;
            counter++;
          }
        }
      }
    }
  }
  
  // Create FieldContainer from data array and index array and show it
  FieldContainer<double> myNewContainer(dimension, dataTeuchosArray);
  std::cout << myNewContainer;

  // Access by overloaded (), multiindex and []:
  multiIndex.resize(myNewContainer.rank());
  multiIndex[0] = 3; 
  multiIndex[1] = 1;
  multiIndex[2] = 2;
  multiIndex[3] = 2;
  multiIndex[4] = 5;
  enumeration = myNewContainer.getEnumeration(multiIndex);
    
  std::cout << "Access by ():          myNewContainer(" << 3 <<"," << 1 << "," << 2 << "," << 2 << "," << 5 << ") = " << myNewContainer(3,1,2,2,5) << "\n";  
  std::cout << "Access by multi-index: myNewContainer{" << multiIndex[0] << multiIndex[1] << multiIndex[2] << multiIndex[3] << multiIndex[4] << "} = "<< myNewContainer.getValue(multiIndex) <<"\n";
  std::cout << "Access by enumeration: myNewContainer[" << enumeration << "] = " << myNewContainer[enumeration] <<"\n";

  std::cout << "Assigning value by (): \n old value at (3,1,2,2,5) = " << myNewContainer(3,1,2,2,5) <<"\n";
  myNewContainer(3,1,2,2,5) = -888.888;
  std::cout << " new value at (3,1,2,2,5) = " << myNewContainer(3,1,2,2,5) <<"\n";
  
  
  std::cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 5: making trivial FieldContainers and storing a single zero         |\n"\
    << "===============================================================================\n\n";
  
  // Make trivial container by resetting the index range to zero rank (no indices) and then
  // using resize method
  dimension.resize(0);
  myContainer.resize(dimension);
  std::cout << myContainer;
  
  // Make trivial container by using clear method:
  myNewContainer.clear();
  std::cout << myNewContainer;
  
  // Now use initialize() to reset the container to hold a single zero
  myNewContainer.initialize();
  std::cout << myNewContainer;
  
  
  std::cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 6: Timing read and write operations using () and getValue           |\n"\
    << "===============================================================================\n\n";
  
  // Initialize dimensions for rank-5 multi-index value
  int dim0 = 10;     // number of cells
  int dim1 = 50;     // number of points
  int dim2 = 27;     // number of functions
  int dim3 = 3;      // 1st space dim
  int dim4 = 3;      // 2nd space dim
  
  FieldContainer<double> myTensorContainer(dim0, dim1, dim2, dim3, dim4);
  multiIndex.resize(myTensorContainer.rank());
  double aValue;
  
  Teuchos::Time timerGetValue("Reading and writing from rank-5 container using getValue");
  timerGetValue.start();
  for(int i0 = 0; i0 < dim0; i0++){
    multiIndex[0] = i0;
    for(int i1 = 0; i1 < dim1; i1++){ 
      multiIndex[1] = i1;
      for(int i2 = 0; i2 < dim2; i2++) {
        multiIndex[2] = i2;
        for(int i3 = 0; i3 < dim3; i3++) {
          multiIndex[3] = i3;
          for(int i4 =0; i4 < dim4; i4++) { 
            multiIndex[4] = i4;
            
            aValue = myTensorContainer.getValue(multiIndex);
            myTensorContainer.setValue(999.999,multiIndex);
            
          }
        }
      }
    }
  }
  timerGetValue.stop();
  std::cout << " Time to read and write from container using getValue: " << timerGetValue.totalElapsedTime() <<"\n";
  
  Teuchos::Time timerRound("Reading and writing from rank-5 container using ()");
  timerRound.start();
  for(int i0 = 0; i0 < dim0; i0++){
    for(int i1 = 0; i1 < dim1; i1++) { 
      for(int i2 = 0; i2 < dim2; i2++) {
        for(int i3 = 0; i3 < dim3; i3++) {
          for(int i4 =0; i4 < dim4; i4++) { 
            
            aValue = myTensorContainer(i0,i1,i2,i3,i4);
            myTensorContainer(i0,i1,i2,i3,i4) = 999.999;
            
          }
        }
      }
    }
  }
  timerRound.stop();
  std::cout << " Time to read and write from container using (): " << timerRound.totalElapsedTime() <<"\n";
    
  std::cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 6: Specialized methods of FieldContainer                            |\n"\
    << "===============================================================================\n\n";
  
   return 0;
}





