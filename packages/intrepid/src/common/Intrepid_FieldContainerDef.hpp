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

/** \file   Intrepid_FieldContainer.hpp
    \brief  Definition file for utility class to provide lexicographical containers.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

  
  //--------------------------------------------------------------------------------------------//
  //                                                                                            //
  //                 Member function definitions of the class FieldContainer                    //
  //                                                                                            //
  //--------------------------------------------------------------------------------------------//
  

template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const FieldContainer<Scalar>& right) {
  
  // Copy dimensions and data values from right
  dimensions_.assign(right.dimensions_.begin(),right.dimensions_.end());  
  data_.assign(right.data_.begin(),right.data_.end());
}

//--------------------------------------------------------------------------------------------//
//                                                                                            //
//                              Constructors of FieldContainer class                          //
//                                                                                            //
//--------------------------------------------------------------------------------------------//


template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const int dim0) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative dimension.");

#endif
  dimensions_.resize(1); 
  dimensions_[0] = dim0;  dim0_ = dim0;
  data_.resize(dim0);
}



template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const int dim0,
                                       const int dim1) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 1st dimension.");
  TEST_FOR_EXCEPTION( (0 > dim1), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 2nd dimension.");
  
#endif
  dimensions_.resize(2); 
  dimensions_[0] = dim0;  dim0_ = dim0; 
  dimensions_[1] = dim1;  dim1_ = dim1;
  data_.resize(dim0*dim1);
}



template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const int dim0,
                                       const int dim1,
                                       const int dim2) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 1st dimension.");
  TEST_FOR_EXCEPTION( (0 > dim1), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 2nd dimension.");
  TEST_FOR_EXCEPTION( (0 > dim2), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 3rd dimension.");
#endif
  dimensions_.resize(3); 
  dimensions_[0] = dim0;  dim0_ = dim0; 
  dimensions_[1] = dim1;  dim1_ = dim1;
  dimensions_[2] = dim2;  dim2_ = dim2;
  data_.resize(dim0*dim1*dim2);
}



template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const int dim0,
                                       const int dim1,
                                       const int dim2,
                                       const int dim3) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 1st dimension.");
  TEST_FOR_EXCEPTION( (0 > dim1), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 2nd dimension.");
  TEST_FOR_EXCEPTION( (0 > dim2), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 3rd dimension.");
  TEST_FOR_EXCEPTION( (0 > dim3), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 4th dimension.");  
#endif
  dimensions_.resize(4); 
  dimensions_[0] = dim0;  dim0_ = dim0; 
  dimensions_[1] = dim1;  dim1_ = dim1;
  dimensions_[2] = dim2;  dim2_ = dim2;
  dimensions_[3] = dim3;  dim3_ = dim3;
  data_.resize(dim0*dim1*dim2*dim3);
}



template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const int dim0,
                                       const int dim1,
                                       const int dim2,
                                       const int dim3,
                                       const int dim4) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 1st dimension.");
  TEST_FOR_EXCEPTION( (0 > dim1), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 2nd dimension.");
  TEST_FOR_EXCEPTION( (0 > dim2), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 3rd dimension.");
  TEST_FOR_EXCEPTION( (0 > dim3), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 4th dimension.");  
  TEST_FOR_EXCEPTION( (0 > dim4), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 5th dimension.");  
#endif
  dimensions_.resize(5); 
  dimensions_[0] = dim0;  dim0_ = dim0; 
  dimensions_[1] = dim1;  dim1_ = dim1;
  dimensions_[2] = dim2;  dim2_ = dim2;
  dimensions_[3] = dim3;  dim3_ = dim3;
  dimensions_[4] = dim4;  dim4_ = dim4;
  data_.resize(dim0*dim1*dim2*dim3*dim4);
}



template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const Teuchos::Array<int>& dimensions) {  
  
#ifdef HAVE_INTREPID_DEBUG
  for(unsigned int dim = 0; dim < dimensions.size(); dim++) {
    TEST_FOR_EXCEPTION( (0 > dimensions[dim] ), std::invalid_argument,  
                        ">>> ERROR (FieldContainer): One or more negative dimensions");  
  }
#endif
  
  // Copy dimensions and resize container storage to match them
  dimensions_.assign(dimensions.begin(),dimensions.end());  
  
  // Copy first 5 dimensions to optimize class for low rank containers
  unsigned int rank = dimensions_.size();
  switch(rank) {
    case 1:
      dim0_ = dimensions_[0]; 
      break;
      
    case 2:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      break;
      
    case 3:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      break;
      
    case 4:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      dim3_ = dimensions_[3]; 
      break;
      
    case 5:
    default:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      dim3_ = dimensions_[3]; 
      dim4_ = dimensions_[4]; 
  }

  // resize data array according to specified dimensions
  data_.resize( this -> getSize());
  
}



template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const Teuchos::Array<int>&    dimensions,
                                       const Teuchos::Array<Scalar>& data) {
 
  // Copy all dimensions
  dimensions_.assign(dimensions.begin(),dimensions.end());
  
  // Copy first 5 dimensions to optimize class for low rank containers
  unsigned int rank = dimensions_.size();
  switch(rank) {
    case 1:
      dim0_ = dimensions_[0]; 
      break;
      
    case 2:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      break;
      
    case 3:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      break;
      
    case 4:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      dim3_ = dimensions_[3]; 
      break;
      
    case 5:
    default:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      dim3_ = dimensions_[3]; 
      dim4_ = dimensions_[4]; 
  }
  
    // Validate input: size of data array must match container size specified by its dimensions
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (int)data.size() != this->getSize() ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Size of input data does not match size of this container.");
#endif
  
  // Assign data
  data_.assign(data.begin(),data.end());
}


//--------------------------------------------------------------------------------------------//
//                                                                                            //
//                            Access methods of FieldContainer class                          //
//                                                                                            //
//--------------------------------------------------------------------------------------------//


template<class Scalar>
inline int FieldContainer<Scalar>::getRank() const {
  return dimensions_.size();
}
  


template<class Scalar>
int FieldContainer<Scalar>::getSize() const {
  
  // Size equals product of all dimensions stored in dimensions_
  int rank = dimensions_.size();
  if(rank == 0) {
    return 0;
  }
  else {
    int size = dim0_;
    
    // Compute size directly to optimize method for low rank (<=5) containers
    switch(rank) {
      case 5:
        size *= dim1_*dim2_*dim3_*dim4_;
        break;
        
      case 4:
        size *= dim1_*dim2_*dim3_;
        break;
        
      case 3:
        size *= dim1_*dim2_;
        break;
        
      case 2:
        size *= dim1_;
        break;
        
      case 1:
        break;
        
      // Compute size for containers with ranks hihger than 5
      default:
        for(int r = 1; r < rank ; r++){
          size *= dimensions_[r];
        }
    }
    return size;
  }
}



template<class Scalar>
inline void FieldContainer<Scalar>::getAllDimensions(Teuchos::Array<int>& dimensions) const {
  dimensions.assign(dimensions_.begin(),dimensions_.end());
}



template<class Scalar>
inline int FieldContainer<Scalar>::getDimension(const int whichDim) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (0 > whichDim),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): dimension order cannot be negative");
  TEST_FOR_EXCEPTION( (whichDim >= this -> getRank() ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): dimension order cannot exceed rank of the container");
#endif
  return dimensions_[whichDim];
}



template<class Scalar>
inline int FieldContainer<Scalar>::getEnumeration(const int i0,
                                                  const int i1) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( ( i0 < 0) || ( i0 >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
#endif
  return i0*dim1_ + i1;
}



template<class Scalar>
inline int FieldContainer<Scalar>::getEnumeration(const int i0,
                                                  const int i1,
                                                  const int i2) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( ( i0 < 0) || ( i0 >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
#endif
  return (i0*dim1_ + i1)*dim2_ + i2;
}



template<class Scalar>
inline int FieldContainer<Scalar>::getEnumeration(const int i0,
                                                  const int i1,
                                                  const int i2,
                                                  const int i3) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( ( i0 < 0) || ( i0 >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");    
#endif
  return ( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3;
}



template<class Scalar>
inline int FieldContainer<Scalar>::getEnumeration(const int i0,
                                                  const int i1,
                                                  const int i2,
                                                  const int i3,
                                                  const int i4) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( ( i0 < 0) || ( i0 >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i4 < 0) || (i4 >= dim4_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 5th index is out of range.");    
#endif
  return ( ( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3 )*dim4_ + i4;
}




template<class Scalar>
int FieldContainer<Scalar>::getEnumeration(const Teuchos::Array<int>& multiIndex) const {

  // Check if empty multi-index.
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( multiIndex.size() == 0 ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Empty multiIndex!");

  // Check if number of multi-indices matches rank of the FieldContainer object
  TEST_FOR_EXCEPTION( ( multiIndex.size() != dimensions_.size() ),
                      std::length_error,
                      ">>> ERROR (FieldContainer): Number of multi-indices does not match rank of container.");
  
  TEST_FOR_EXCEPTION( ( ( multiIndex[0] < 0) || ( multiIndex[0] >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
#endif  

  int rank = dimensions_.size();
  int address = 0;
  switch (rank) {
    
    // Optimize enumeration computation for low rank (<= 5) containers
    case 5:
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( (multiIndex[4] < 0) || (multiIndex[4] >= dim4_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 5th index is out of range.");    
      TEST_FOR_EXCEPTION( ( (multiIndex[3] < 0) || (multiIndex[3] >= dim3_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 4th index is out of range.");    
      TEST_FOR_EXCEPTION( ( (multiIndex[2] < 0) || (multiIndex[2] >= dim2_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 3rd index is out of range.");    
      TEST_FOR_EXCEPTION( ( (multiIndex[1] < 0) || (multiIndex[1] >= dim1_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 2nd index is out of range.");    
#endif
      address = (((multiIndex[0]*dim1_ + multiIndex[1])*dim2_ + multiIndex[2])*dim3_ + multiIndex[3])*dim4_ + multiIndex[4];
      break;
      
    case 4:
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( (multiIndex[3] < 0) || (multiIndex[3] >= dim3_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 4th index is out of range.");    
      TEST_FOR_EXCEPTION( ( (multiIndex[2] < 0) || (multiIndex[2] >= dim2_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 3rd index is out of range.");    
      TEST_FOR_EXCEPTION( ( (multiIndex[1] < 0) || (multiIndex[1] >= dim1_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 2nd index is out of range.");    
#endif
      address = ((multiIndex[0]*dim1_ + multiIndex[1])*dim2_ + multiIndex[2])*dim3_ + multiIndex[3];
      break;

    case 3:
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( (multiIndex[2] < 0) || (multiIndex[2] >= dim2_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 3rd index is out of range.");    
      TEST_FOR_EXCEPTION( ( (multiIndex[1] < 0) || (multiIndex[1] >= dim1_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 2nd index is out of range.");    
#endif
      address = (multiIndex[0]*dim1_ + multiIndex[1])*dim2_ + multiIndex[2];
      break;

    case 2:
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( (multiIndex[1] < 0) || (multiIndex[1] >= dim1_) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): 2nd index is out of range.");    
#endif
      address = multiIndex[0]*dim1_ + multiIndex[1];
      break;
      
    case 1:
      address = multiIndex[0];
      break;

    default:
      
      // Arbitrary rank: compute address using Horner's nested scheme: intialize address to 0th index value
      address = multiIndex[0];
      for (int r = 0; r < rank - 1; r++){
#ifdef HAVE_INTREPID_DEBUG
        TEST_FOR_EXCEPTION( ( (multiIndex[r+1] < 0) || (multiIndex[r+1] >= dimensions_[r+1]) ),
                            std::invalid_argument,
                            ">>> ERROR (FieldContainer): Multi-index component out of range.");    
#endif
        // Add increment
        address = address*dimensions_[r+1] + multiIndex[r+1];
      }
  } // end switch(rank)

  return address;
}



template<class Scalar>
void FieldContainer<Scalar>::getMultiIndex(Teuchos::Array<int>& multiIndex,
                                           const int            valueEnum) const 
{
  
  // Verify address is in the admissible range for this FieldContainer
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Specified address is out of range.");    
#endif
  
  // make sure multiIndex has the right size to hold all multi-indices
  int rank = dimensions_.size();
  multiIndex.resize(rank);
  
  // Initializations
  int temp_enum = valueEnum;
  int temp_range = 1;
  
  // Compute product of all but the first upper bound
  for(int r = 1; r < rank ; r++){
    temp_range *=dimensions_[r];
  }
  
  // Index 0 is computed first using integer division
  multiIndex[0] = temp_enum/temp_range;
  
  // Indices 1 to (rank - 2) are computed next; will be skipped if rank <=2
  for(int r = 1; r < rank - 1; r++){
    temp_enum  -= multiIndex[r-1]*temp_range;
    temp_range /= dimensions_[r];
    multiIndex[r] = temp_enum/temp_range;
  }
  
  // Index (rank - 1) is computed last, skip if rank = 1 and keep if rank = 2
  if(rank > 1) {
    multiIndex[rank - 1] = temp_enum - multiIndex[rank - 2]*temp_range;
  }
}


//--------------------------------------------------------------------------------------------//
//                                                                                            //
//                          Methods to shape (resize) a field container                       //
//                                                                                            //
//--------------------------------------------------------------------------------------------//

template<class Scalar>
inline void FieldContainer<Scalar>::empty() {
  dimensions_.resize(0);
  
  // Reset first five dimensions:
  dim0_ = 0;
  dim1_ = 0;
  dim2_ = 0;
  dim3_ = 0;
  dim4_ = 0;
  
  // Resize data array to zero length
  data_.resize(0);
}



template<class Scalar>
void FieldContainer<Scalar>::resize(const Teuchos::Array<int>& newDimensions) {
  
  // First handle the trivial case of zero dimensions
  if( newDimensions.size() == 0) {
    dimensions_.resize(0);
    dim0_ = 0;
    dim1_ = 0;
    dim2_ = 0;
    dim3_ = 0;
    dim4_ = 0;
    data_.resize(0);
  }
  else {
    
    // Copy upper index bounds and resize container storage to match new upper bounds.
    dimensions_.assign(newDimensions.begin(),newDimensions.end());  
    
    // Copy first 5 dimensions for faster access
    unsigned int rank = dimensions_.size();
    switch(rank) {
      case 1:
        dim0_ = dimensions_[0]; 
        break;
        
      case 2:
        dim0_ = dimensions_[0]; 
        dim1_ = dimensions_[1]; 
        break;
        
      case 3:
        dim0_ = dimensions_[0]; 
        dim1_ = dimensions_[1]; 
        dim2_ = dimensions_[2]; 
        break;
        
      case 4:
        dim0_ = dimensions_[0]; 
        dim1_ = dimensions_[1]; 
        dim2_ = dimensions_[2]; 
        dim3_ = dimensions_[3]; 
        break;
        
      case 5:
      default:
        dim0_ = dimensions_[0]; 
        dim1_ = dimensions_[1]; 
        dim2_ = dimensions_[2]; 
        dim3_ = dimensions_[3]; 
        dim4_ = dimensions_[4]; 
    }
    
    // Resize data array
    data_.resize(this->getSize());
  }
}



template<class Scalar>
inline void FieldContainer<Scalar>::resize(const int dim0) {
  dimensions_.resize(1);  
  dimensions_[0] = dim0;  dim0_ = dim0;
  data_.resize(dim0); 
}



template<class Scalar>
inline void FieldContainer<Scalar>::resize(const int dim0,
                                           const int dim1) {
  dimensions_.resize(2);  
  dimensions_[0] = dim0;  dim0_ = dim0;
  dimensions_[1] = dim1;  dim1_ = dim1;
  data_.resize(dim0*dim1); 
}



template<class Scalar>
inline void FieldContainer<Scalar>::resize(const int dim0,
                                           const int dim1,
                                           const int dim2) {
  dimensions_.resize(3);  
  dimensions_[0] = dim0;  dim0_ = dim0;
  dimensions_[1] = dim1;  dim1_ = dim1;
  dimensions_[2] = dim2;  dim2_ = dim2;
  data_.resize(dim0*dim1*dim2); 
}



template<class Scalar>
inline void FieldContainer<Scalar>::resize(const int dim0,
                                           const int dim1,
                                           const int dim2,
                                           const int dim3) {
  dimensions_.resize(4);  
  dimensions_[0] = dim0;  dim0_ = dim0;
  dimensions_[1] = dim1;  dim1_ = dim1;
  dimensions_[2] = dim2;  dim2_ = dim2;
  dimensions_[3] = dim3;  dim3_ = dim3;
  data_.resize(dim0*dim1*dim2*dim3); 
}



template<class Scalar>
inline void FieldContainer<Scalar>::resize(const int dim0,
                                           const int dim1,
                                           const int dim2,
                                           const int dim3,
                                           const int dim4) {
  dimensions_.resize(5);  
  dimensions_[0] = dim0;  dim0_ = dim0;
  dimensions_[1] = dim1;  dim1_ = dim1;
  dimensions_[2] = dim2;  dim2_ = dim2;
  dimensions_[3] = dim3;  dim3_ = dim3;
  dimensions_[4] = dim4;  dim4_ = dim4;
  data_.resize(dim0*dim1*dim2*dim3*dim4); 
}



template<class Scalar>
inline void FieldContainer<Scalar>::resize(const FieldContainer<Scalar>& anotherContainer) {
  
  int newRank = anotherContainer.getRank();
  dimensions_.resize(newRank);
  for(int i = 0; i < newRank; i++){
    dimensions_[i] = anotherContainer.getDimension(i);
  }
  
  // Copy first 5 dimensions for faster access
  switch(newRank) {
    case 1:
      dim0_ = dimensions_[0]; 
      break;
      
    case 2:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      break;
      
    case 3:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      break;
      
    case 4:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      dim3_ = dimensions_[3]; 
      break;
      
    case 5:
    default:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      dim3_ = dimensions_[3]; 
      dim4_ = dimensions_[4]; 
  }
  
  // Resize data array
  data_.resize(this->getSize());  
}


//--------------------------------------------------------------------------------------------//
//                                                                                            //
//                     Methods to read and write values to FieldContainer                     //
//                                                                                            //
//--------------------------------------------------------------------------------------------//


template<class Scalar>
inline void FieldContainer<Scalar>::storeZero() {
  for (int i=0; i < this->getSize(); i++) {
    data_[i] = (Scalar)0;
  } 
}



template<class Scalar>
inline Scalar FieldContainer<Scalar>::getValue(const Teuchos::Array<int>& multiIndex) const {
  return data_[this -> getEnumeration(multiIndex)];
}



template<class Scalar>
inline void FieldContainer<Scalar>::setValue(const Scalar dataValue, 
                                             const Teuchos::Array<int>& multiIndex) {
  data_[this -> getEnumeration(multiIndex)] = dataValue; 
}



template<class Scalar>
inline void FieldContainer<Scalar>::setValue(const Scalar dataValue, 
                                             const int    order) {
  data_[order] = dataValue; 
}



template<class Scalar>
void FieldContainer<Scalar>::setValues(const Teuchos::Array<Scalar>& dataArray) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (dataArray.size() != (data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Size of argument does not match size of container.");  
#endif  
  data_.assign(dataArray.begin(),dataArray.end());
}



template<class Scalar>
inline const Scalar& FieldContainer<Scalar>::operator () (const int i0,
                                                          const int i1) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 2), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( ( i0 < 0) || ( i0 >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
#endif
  return data_[i0*dim1_ + i1]; 
}


template<class Scalar>
inline Scalar& FieldContainer<Scalar>::operator () (const int i0,
                                                    const int i1)  
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 2), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( ( i0 < 0) || ( i0 >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
#endif
  return data_[i0*dim1_ + i1]; 
}



template<class Scalar>
inline const Scalar& FieldContainer<Scalar>::operator () (const int i0,
                                                          const int i1,
                                                          const int i2) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 3), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( ( i0 < 0) || ( i0 >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
#endif
  return data_[(i0*dim1_ + i1)*dim2_ + i2]; 
}

template<class Scalar>
inline Scalar& FieldContainer<Scalar>::operator () (const int i0,
                                                    const int i1,
                                                    const int i2) 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 3), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( ( i0 < 0) || ( i0 >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
#endif
  return data_[(i0*dim1_ + i1)*dim2_ + i2]; 
}



template<class Scalar>
inline const Scalar& FieldContainer<Scalar>::operator ()  (const int i0,
                                                           const int i1,
                                                           const int i2,
                                                           const int i3) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( ( i0 < 0) || ( i0 >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");    
#endif
  return data_[( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3];
}


template<class Scalar>
inline Scalar& FieldContainer<Scalar>::operator ()  (const int i0,
                                                     const int i1,
                                                     const int i2,
                                                     const int i3) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( ( i0 < 0) || ( i0 >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");    
#endif
  return data_[( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3];
}



template<class Scalar>
inline const Scalar& FieldContainer<Scalar>::operator ()  (const int i0,
                                                           const int i1,
                                                           const int i2,
                                                           const int i3,
                                                           const int i4) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( ( i0 < 0) || ( i0 >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i4 < 0) || (i4 >= dim4_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 5th index is out of range.");    
#endif
  return data_[( ( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3 )*dim4_ + i4];
}

template<class Scalar>
inline Scalar& FieldContainer<Scalar>::operator ()  (const int i0,
                                                     const int i1,
                                                     const int i2,
                                                     const int i3,
                                                     const int i4) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( ( i0 < 0) || ( i0 >= dim0_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i4 < 0) || (i4 >= dim4_) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): 5th index is out of range.");    
#endif
  return data_[( ( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3 )*dim4_ + i4];
}



template<class Scalar>
const Scalar& FieldContainer<Scalar>::operator [] (const int address) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (address < 0) || (address >= (int)data_.size() ) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Specified address is out of range.");
#endif
  return data_[address];
}
  


template<class Scalar>
inline FieldContainer<Scalar>& FieldContainer<Scalar>::operator = (const FieldContainer<Scalar>& right)
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this == &right ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Invalid right-hand side to '='. Self-assignment prohibited.");
#endif
  data_ = right.data_;
  dimensions_ = right.dimensions_; 
  return *this;
}


//===========================================================================//
//                                                                           //
//           END of member definitions; START friends and related            //
//                                                                           //
//===========================================================================//


template<class Scalar>
std::ostream& operator << (std::ostream& os, const FieldContainer<Scalar>& container) {
  
  // Save the format state of the original ostream os.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(os);  

  os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  os.setf(std::ios_base::right);
  int myprec = os.precision();
  
  int size = container.getSize();
  int rank = container.getRank();
  Teuchos::Array<int> multiIndex(rank);
  Teuchos::Array<int> dimensions;
  container.getAllDimensions(dimensions);
  
  os<< "===============================================================================\n"\
    << "\t Container size = " << size << "   rank = " << rank << "\n" ;
  
  if( (rank == 0 ) && (size == 0) ) {
    os<< "\t Index Range    = (0) \n"    
      << "===============================================================================\n"\
      << "|                     *** This is an empty container ****                     |\n";
  }
  else {
    os<< "\t Index Range    = ";
    
    for(int r = 0; r < rank; r++){
      os << " (" << dimensions[r] <<") ";
    }
    os << "\n";
    
    os<< "===============================================================================\n"\
      << "| \t Multi-index        Address              Value                            |\n"\
      << "===============================================================================\n";
  }
  
  for(int address = 0; address < size; address++){
    container.getMultiIndex(multiIndex,address);
    os << "\t\t" ;
    for(int r = 0; r < rank; r++){
      os <<  multiIndex[r]; 
    }
    os<< std::setiosflags(std::ios::left) << std::setw(16) << address << "\t" \
      << std::setw(myprec+8) << "\t" << container[address] << "\n";
  }
  
  os<< "===============================================================================\n\n";

  // reset format state of os
  os.copyfmt(oldFormatState);

  return os;
}

// End member, friend, and related function definitions of class FieldContainer.

} // end namespace Intrepid
