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
  dim0_ = right.dim0_;
  dim1_ = right.dim1_;
  dim2_ = right.dim2_;
  dim3_ = right.dim3_;
  dim4_ = right.dim4_;
}

//--------------------------------------------------------------------------------------------//
//                                                                                            //
//                              Constructors of FieldContainer class                          //
//                                                                                            //
//--------------------------------------------------------------------------------------------//


template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const int dim0) : dim0_(dim0), dim1_(0), dim2_(0), dim3_(0), dim4_(0) 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative dimension.");

#endif
  dimensions_.resize(1); 
  dimensions_[0] = dim0_;  
  data_.resize(dim0_);
}



template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const int dim0,
                                       const int dim1) : dim0_(dim0), dim1_(dim1), dim2_(0), dim3_(0), dim4_(0)
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 1st dimension.");
  TEST_FOR_EXCEPTION( (0 > dim1), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 2nd dimension.");
  
#endif
  dimensions_.resize(2); 
  dimensions_[0] = dim0_;   
  dimensions_[1] = dim1_;  
  data_.resize(dim0_*dim1_);
}



template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const int dim0,
                                       const int dim1,
                                       const int dim2) : dim0_(dim0), dim1_(dim1), dim2_(dim2), dim3_(0), dim4_(0)
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (0 > dim0), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 1st dimension.");
  TEST_FOR_EXCEPTION( (0 > dim1), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 2nd dimension.");
  TEST_FOR_EXCEPTION( (0 > dim2), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): FieldContainer cannot have a negative 3rd dimension.");
#endif
  dimensions_.resize(3); 
  dimensions_[0] = dim0_;  
  dimensions_[1] = dim1_; 
  dimensions_[2] = dim2_;  
  data_.resize(dim0_*dim1_*dim2_);
}



template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const int dim0,
                                       const int dim1,
                                       const int dim2,
                                       const int dim3) : dim0_(dim0), dim1_(dim1), dim2_(dim2), dim3_(dim3), dim4_(0)
{
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
  dimensions_[0] = dim0_;  
  dimensions_[1] = dim1_; 
  dimensions_[2] = dim2_;  
  dimensions_[3] = dim3_; 
  data_.resize(dim0_*dim1_*dim2_*dim3_);
}



template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const int dim0,
                                       const int dim1,
                                       const int dim2,
                                       const int dim3,
                                       const int dim4) : dim0_(dim0), dim1_(dim1), dim2_(dim2), dim3_(dim3), dim4_(dim4)
{
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
  dimensions_[0] = dim0_;
  dimensions_[1] = dim1_;
  dimensions_[2] = dim2_;  
  dimensions_[3] = dim3_;  
  dimensions_[4] = dim4_;  
  data_.resize(dim0_*dim1_*dim2_*dim3_*dim4_);
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
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;
      
    case 2:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;
      
    case 3:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      dim3_ = 0;
      dim4_ = 0;
      break;
      
    case 4:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      dim3_ = dimensions_[3]; 
      dim4_ = 0;
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
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;
      
    case 2:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;
      
    case 3:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      dim3_ = 0;
      dim4_ = 0;
      break;
      
    case 4:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      dim3_ = dimensions_[3]; 
      dim4_ = 0;
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
  TEST_FOR_EXCEPTION( ( (int)data.size() != this -> getSize() ),
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
  // Important! This method is used by constructors to find out what is the needed size of data_
  // based on the specified dimensions. Therefore, it cannot be implmented by returning data_.size
  // and must be able to compute the size of the container based only on its specified dimensions
  
  // Size equals product of all dimensions stored in dimensions_
  int rank = dimensions_.size();
  
  // If container has no dimensions its size is zero
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
  TEST_FOR_EXCEPTION( (0 > whichDim), std::invalid_argument,
                      ">>> ERROR (FieldContainer): dimension order cannot be negative");
  TEST_FOR_EXCEPTION( (whichDim >= this -> getRank() ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): dimension order cannot exceed rank of the container");
#endif
  return dimensions_[whichDim];
}



template<class Scalar>
inline int FieldContainer<Scalar>::getEnumeration(const int i0) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 1), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (i0 < 0) || ( i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): index is out of range.");
#endif
  return i0;
}



template<class Scalar>
inline int FieldContainer<Scalar>::getEnumeration(const int i0,
                                                  const int i1) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 2), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (i0 < 0) || ( i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
#endif
  return i0*dim1_ + i1;
}



template<class Scalar>
inline int FieldContainer<Scalar>::getEnumeration(const int i0,
                                                  const int i1,
                                                  const int i2) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 3), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (i0 < 0) || ( i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
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
  TEST_FOR_EXCEPTION( ( this -> getRank() != 4), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (i0 < 0) || ( i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ), std::invalid_argument,
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
  TEST_FOR_EXCEPTION( ( this -> getRank() != 5), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (i0 < 0) || ( i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i4 < 0) || (i4 >= dim4_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 5th index is out of range.");    
#endif
  return ( ( (i0*dim1_ + i1 )*dim2_ + i2 )*dim3_ + i3 )*dim4_ + i4;
}




template<class Scalar>
int FieldContainer<Scalar>::getEnumeration(const Teuchos::Array<int>& multiIndex) const {

#ifdef HAVE_INTREPID_DEBUG
  // Check if number of multi-indices matches rank of the FieldContainer object
  TEST_FOR_EXCEPTION( ( multiIndex.size() != dimensions_.size() ),
                      std::invalid_argument,
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
void FieldContainer<Scalar>::getMultiIndex(int & i0,
                                           const int valueEnum) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 1), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Value enumeration is out of range.");    
#endif
  i0 = valueEnum;
}



template<class Scalar>
void FieldContainer<Scalar>::getMultiIndex(int & i0,
                                           int & i1,
                                           const int valueEnum) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 2), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Value enumeration is out of range.");    
#endif
  
  i0 = valueEnum/dim1_;
  i1 = valueEnum - i0*dim1_;
}



template<class Scalar>
void FieldContainer<Scalar>::getMultiIndex(int & i0,
                                           int & i1,
                                           int & i2,
                                           const int valueEnum) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 3), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Value enumeration is out of range.");    
#endif
  int tempDim = dim1_*dim2_;
  int tempEnu = valueEnum;
  i0 = tempEnu/tempDim;
  
  tempEnu -= i0*tempDim;
  tempDim /= dim1_;
  i1 = tempEnu/tempDim;
  
  tempEnu -= i1*tempDim;
  tempDim /= dim2_;
  i2 = tempEnu/tempDim;
}



template<class Scalar>
void FieldContainer<Scalar>::getMultiIndex(int & i0,
                                           int & i1,
                                           int & i2,
                                           int & i3,
                                           const int valueEnum) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 4), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Value enumeration is out of range.");    
#endif
  int tempDim = dim1_*dim2_*dim3_;
  int tempEnu = valueEnum;
  i0 = tempEnu/tempDim;
  
  tempEnu -= i0*tempDim;
  tempDim /= dim1_;
  i1 = tempEnu/tempDim;
  
  tempEnu -= i1*tempDim;
  tempDim /= dim2_;
  i2 = tempEnu/tempDim;
  
  tempEnu -= i2*tempDim;
  tempDim /= dim3_;
  i3 = tempEnu/tempDim;
}




template<class Scalar>
void FieldContainer<Scalar>::getMultiIndex(int & i0,
                                           int & i1,
                                           int & i2,
                                           int & i3,
                                           int & i4,
                                           const int valueEnum) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 5), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Value enumeration is out of range.");    
#endif
  int tempDim = dim1_*dim2_*dim3_*dim4_;
  int tempEnu = valueEnum;
  i0 = tempEnu/tempDim;
  
  tempEnu -= i0*tempDim;
  tempDim /= dim1_;
  i1 = tempEnu/tempDim;
  
  tempEnu -= i1*tempDim;
  tempDim /= dim2_;
  i2 = tempEnu/tempDim;
  
  tempEnu -= i2*tempDim;
  tempDim /= dim3_;
  i3 = tempEnu/tempDim;

  tempEnu -= i3*tempDim;
  tempDim /= dim4_;
  i4 = tempEnu/tempDim;
}



template<class Scalar>
void FieldContainer<Scalar>::getMultiIndex(Teuchos::Array<int>& multiIndex,
                                           const int            valueEnum) const 
{
  
  // Verify address is in the admissible range for this FieldContainer
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Value enumeration is out of range.");    
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
        dim1_ = 0;
        dim2_ = 0;
        dim3_ = 0;
        dim4_ = 0;
        break;
        
      case 2:
        dim0_ = dimensions_[0]; 
        dim1_ = dimensions_[1]; 
        dim2_ = 0;
        dim3_ = 0;
        dim4_ = 0;
        break;
        
      case 3:
        dim0_ = dimensions_[0]; 
        dim1_ = dimensions_[1]; 
        dim2_ = dimensions_[2]; 
        dim3_ = 0;
        dim4_ = 0;
        break;
        
      case 4:
        dim0_ = dimensions_[0]; 
        dim1_ = dimensions_[1]; 
        dim2_ = dimensions_[2]; 
        dim3_ = dimensions_[3]; 
        dim4_ = 0;
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
    data_.resize(this -> getSize());
  }
}



template<class Scalar>
inline void FieldContainer<Scalar>::resize(const int dim0) {
  dim0_ = dim0;  
  dim1_ = 0;  
  dim2_ = 0;  
  dim3_ = 0;  
  dim4_ = 0;
  dimensions_.resize(1);  
  dimensions_[0] = dim0_; 
  data_.resize(dim0_); 
}



template<class Scalar>
inline void FieldContainer<Scalar>::resize(const int dim0,
                                           const int dim1) {
  dim0_ = dim0;  
  dim1_ = dim1;  
  dim2_ = 0;  
  dim3_ = 0;  
  dim4_ = 0;
  dimensions_.resize(2);  
  dimensions_[0] = dim0_;  
  dimensions_[1] = dim1_;  
  data_.resize(dim0_*dim1_); 
}



template<class Scalar>
inline void FieldContainer<Scalar>::resize(const int dim0,
                                           const int dim1,
                                           const int dim2) {
  dim0_ = dim0;
  dim1_ = dim1;
  dim2_ = dim2;
  dim3_ = 0;  
  dim4_ = 0;
  dimensions_.resize(3);  
  dimensions_[0] = dim0_; 
  dimensions_[1] = dim1_;  
  dimensions_[2] = dim2_;  
  data_.resize(dim0_*dim1_*dim2_); 
}



template<class Scalar>
inline void FieldContainer<Scalar>::resize(const int dim0,
                                           const int dim1,
                                           const int dim2,
                                           const int dim3) {
  dim0_ = dim0;
  dim1_ = dim1;
  dim2_ = dim2;
  dim3_ = dim3;
  dim4_ = 0;
  dimensions_.resize(4);  
  dimensions_[0] = dim0_;  
  dimensions_[1] = dim1_;  
  dimensions_[2] = dim2_;  
  dimensions_[3] = dim3_;  
  data_.resize(dim0_*dim1_*dim2_*dim3_); 
}



template<class Scalar>
inline void FieldContainer<Scalar>::resize(const int dim0,
                                           const int dim1,
                                           const int dim2,
                                           const int dim3,
                                           const int dim4) {
  dim0_ = dim0;
  dim1_ = dim1;
  dim2_ = dim2;
  dim3_ = dim3;
  dim4_ = dim4;
  dimensions_.resize(5);  
  dimensions_[0] = dim0_;  
  dimensions_[1] = dim1_;  
  dimensions_[2] = dim2_;  
  dimensions_[3] = dim3_;  
  dimensions_[4] = dim4_;  
  data_.resize(dim0_*dim1_*dim2_*dim3_*dim4_); 
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
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;
      
    case 2:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;
      
    case 3:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      dim3_ = 0;
      dim4_ = 0;
      break;
      
    case 4:
      dim0_ = dimensions_[0]; 
      dim1_ = dimensions_[1]; 
      dim2_ = dimensions_[2]; 
      dim3_ = dimensions_[3]; 
      dim4_ = 0;
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


template<class Scalar>
void FieldContainer<Scalar>::resize(const int       numPoints,
                                    const int       numFields,
                                    const EField    fieldType,
                                    const EOperator operatorType,
                                    const int       spaceDim) {  
  // Validate input
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( numPoints < 0),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of points cannot be negative!");  
  TEST_FOR_EXCEPTION( ( numFields < 0),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of fields cannot be negative!");  
  TEST_FOR_EXCEPTION( !( (1 <=  spaceDim ) && ( spaceDim <= 3  ) ),
                      std::invalid_argument,
                      ">>> ERROR (FieldContainer): Invalid space dimension.");  
#endif  
  
  // Find out field and operator ranks
  int fieldRank    = getFieldRank(fieldType);
  int operatorRank = getOperatorRank(operatorType,fieldRank,spaceDim);
  
  // Compute rank of the container = 1(numPoints) + 1(numFields) + fieldRank + operatorRank
  int rank = 1 + 1 + fieldRank + operatorRank;
  
  // Define temp array for the dimensions
  Teuchos::Array<int> newDimensions(rank);
  
  // Dimensions 0 and 1 are number of points and number of fields, resp.
  newDimensions[0] = numPoints;
  newDimensions[1] = numFields;
  
  // The rest of the dimensions depend on whether we had VALUE, GRAD (D1), CURL, DIV or Dk, k>1
  switch(operatorType) {
    
    case OPERATOR_VALUE:
    case OPERATOR_GRAD:
    case OPERATOR_D1:
    case OPERATOR_CURL:
    case OPERATOR_DIV:
      
      // For these operators all dimensions from 2 to 2 + fieldRank + OperatorRank are bounded by spaceDim
      for(int i = 0; i < fieldRank + operatorRank; i++){
        newDimensions[2 + i] = spaceDim; 
      }
      break;
      
    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      
      // All dimensions from 2 to 2 + fieldRank, if any, are bounded by spaceDim
      for(int i = 0; i < fieldRank; i++){
        newDimensions[2 + i] = spaceDim; 
      }
      
      // We know that for Dk operatorRank = 1 and so there's just one more dimension left
      // given by the cardinality of the set of all derivatives of order k
      newDimensions[2 + fieldRank] = getDkCardinality(getOperatorOrder(operatorType),spaceDim);
      break;
      
    default:
      TEST_FOR_EXCEPTION( ( (operatorType != OPERATOR_VALUE) &&
                            (operatorType != OPERATOR_GRAD)  && 
                            (operatorType != OPERATOR_CURL)  && 
                            (operatorType != OPERATOR_DIV)   &&
                            (operatorType != OPERATOR_D1)    && 
                            (operatorType != OPERATOR_D2)    && 
                            (operatorType != OPERATOR_D3)    && 
                            (operatorType != OPERATOR_D4)    && 
                            (operatorType != OPERATOR_D5)    && 
                            (operatorType != OPERATOR_D6)    && 
                            (operatorType != OPERATOR_D7)    && 
                            (operatorType != OPERATOR_D8)    && 
                            (operatorType != OPERATOR_D9)    && 
                            (operatorType != OPERATOR_D10) ),
                          std::invalid_argument,
                          ">>> ERROR (FieldContainer): Invalid operator type");    
  }
  
  // Resize FieldContainer using the newDimensions in the array
  this -> resize(newDimensions);
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
inline const Scalar& FieldContainer<Scalar>::operator () (const int i0) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 1), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): index is out of range.");    
#endif
  return data_[i0]; 
}


template<class Scalar>
inline Scalar& FieldContainer<Scalar>::operator () (const int i0)  
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 1), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): index is out of range.");    
#endif
  return data_[i0]; 
}



template<class Scalar>
inline const Scalar& FieldContainer<Scalar>::operator () (const int i0,
                                                          const int i1) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> getRank() != 2), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
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
  TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
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
  TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
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
  TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
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
  TEST_FOR_EXCEPTION( ( this -> getRank() != 4), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ), std::invalid_argument,
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
  TEST_FOR_EXCEPTION( ( this -> getRank() != 4), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ), std::invalid_argument,
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
  TEST_FOR_EXCEPTION( ( this -> getRank() != 5), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i4 < 0) || (i4 >= dim4_) ), std::invalid_argument,
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
  TEST_FOR_EXCEPTION( ( this -> getRank() != 5), std::invalid_argument, 
                      ">>> ERROR (FieldContainer): Number of indices does not match rank of the container.");  
  TEST_FOR_EXCEPTION( ( (i0 < 0) || (i0 >= dim0_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 1st index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i1 < 0) || (i1 >= dim1_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 2nd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i2 < 0) || (i2 >= dim2_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 3rd index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i3 < 0) || (i3 >= dim3_) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): 4th index is out of range.");    
  TEST_FOR_EXCEPTION( ( (i4 < 0) || (i4 >= dim4_) ), std::invalid_argument,
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
  dim0_ = right.dim0_;
  dim1_ = right.dim1_;
  dim2_ = right.dim2_;
  dim3_ = right.dim3_;
  dim4_ = right.dim4_;
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
    << "\t Container size = " << size << "\n"
    << "\t Container rank = " << rank << "\n" ;
  
  if( (rank == 0 ) && (size == 0) ) {
    os<< "===============================================================================\n"\
      << "|                     *** This is an empty container ****                     |\n";
  }
  else {
    os<< "\t Dimensions     = ";
    
    for(int r = 0; r < rank; r++){
      os << " (" << dimensions[r] <<") ";
    }
    os << "\n";
    
    os<< "===============================================================================\n"\
      << "| \t Multi-index        Enumeration               Value                       |\n"\
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


//===========================================================================//
//                                                                           //
//         Global functions for field and operator rank, Dk enumeration,     //
//         cardinality, and conversion of enumertion to multiplicities       //
//                                                                           //
//===========================================================================//


int getFieldRank(const EField fieldType) {
  int fieldRank = -1;
  
  switch(fieldType){
    
    case FIELD_FORM_0:
    case FIELD_FORM_3:
      fieldRank = 0;
      break;
      
    case FIELD_FORM_1:
    case FIELD_FORM_2:
    case FIELD_VECTOR:
      fieldRank = 1;
      break;
      
    case FIELD_TENSOR:
      fieldRank = 2;
      break;
      
    default:
      TEST_FOR_EXCEPTION( ( (fieldType != FIELD_FORM_0) &&
                            (fieldType != FIELD_FORM_1) && 
                            (fieldType != FIELD_FORM_2) && 
                            (fieldType != FIELD_FORM_3) &&
                            (fieldType != FIELD_VECTOR) && 
                            (fieldType != FIELD_TENSOR) ),
                          std::invalid_argument,
                          ">>> ERROR (getFieldRank): Invalid field type");
  }
  return fieldRank;
}



int getOperatorOrder(const EOperator operatorType) {
  int opOrder = -1;
  
  switch(operatorType){
    
    case OPERATOR_VALUE:
      opOrder = 0;
      break;
      
    case OPERATOR_GRAD:
    case OPERATOR_CURL:
    case OPERATOR_DIV:
    case OPERATOR_D1:
      opOrder = 1;
      break;
      
    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      opOrder = (int)operatorType - (int)OPERATOR_D1 + 1;
      break;
      
    default:
      TEST_FOR_EXCEPTION( ( (operatorType != OPERATOR_VALUE) &&
                            (operatorType != OPERATOR_GRAD)  && 
                            (operatorType != OPERATOR_CURL)  && 
                            (operatorType != OPERATOR_DIV)   &&
                            (operatorType != OPERATOR_D1)    && 
                            (operatorType != OPERATOR_D2)    && 
                            (operatorType != OPERATOR_D3)    && 
                            (operatorType != OPERATOR_D4)    && 
                            (operatorType != OPERATOR_D5)    && 
                            (operatorType != OPERATOR_D6)    && 
                            (operatorType != OPERATOR_D7)    && 
                            (operatorType != OPERATOR_D8)    && 
                            (operatorType != OPERATOR_D9)    && 
                            (operatorType != OPERATOR_D10) ),
                          std::invalid_argument,
                          ">>> ERROR (getOperatorOrder): Invalid operator type");
  }
  return opOrder;
}    



int getOperatorRank(const EOperator operatorType,
                    const int       fieldRank,
                    const int       spaceDim) {
  
  // Verify arguments: field rank can be 0,1, or 2, spaceDim can be 1,2, or 3.
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( !( (0 <= fieldRank) && (fieldRank <= 2) ),
                      std::invalid_argument,
                      ">>> ERROR (getOperatorRank): Invalid field rank");
  TEST_FOR_EXCEPTION( !( (1 <= spaceDim ) && (spaceDim  <= 3) ),
                      std::invalid_argument,
                      ">>> ERROR (getOperatorRank): Invalid space dimension");
#endif
  int operatorRank = -999;
  
  // In 1D GRAD, CURL, and DIV default to d/dx; Dk defaults to d^k/dx^k, no casing needed.
  if(spaceDim == 1) {
    if(fieldRank == 0) {
      
      // By default, in 1D any operator other than VALUE has rank 1
      if(operatorType == OPERATOR_VALUE) {
        operatorRank = 0;
      }
      else {
        operatorRank = 1;
      }
    }
    
    // Only scalar fields are allowed in 1D
    else {
      TEST_FOR_EXCEPTION( ( fieldRank > 0 ),
                          std::invalid_argument,
                          ">>> ERROR (getOperatorRank): Only scalar fields are allowed in 1D");  
    } // fieldRank == 0
  } // spaceDim == 1
  
  // We are either in 2D or 3D
  else {  
    switch(operatorType) {
      case OPERATOR_VALUE:
        operatorRank = 0;
        break;
        
      case OPERATOR_GRAD:
      case OPERATOR_D1:
      case OPERATOR_D2:
      case OPERATOR_D3:
      case OPERATOR_D4:
      case OPERATOR_D5:
      case OPERATOR_D6:
      case OPERATOR_D7:
      case OPERATOR_D8:
      case OPERATOR_D9:
      case OPERATOR_D10:
        operatorRank = 1;
        break;
        
      case OPERATOR_CURL:
        
        // operator rank for vector and tensor fields equals spaceDim - 3 (-1 in 2D and 0 in 3D)   
        if(fieldRank > 0) {
          operatorRank = spaceDim - 3;
        }
        else {
          
          // CURL can be applied to scalar functions (rank = 0) in 2D and gives a vector (rank = 1)
          if(spaceDim == 2) {
            operatorRank = 3 - spaceDim;
          }
          
          // If we are here, fieldRank=0, spaceDim=3: CURL is undefined for 3D scalar functions
          else {
            TEST_FOR_EXCEPTION( ( (spaceDim == 3) && (fieldRank == 0) ),
                                std::invalid_argument,
                                ">>> ERROR (getOperatorRank): CURL cannot be applied to scalar fields in 3D");  
          }
        }
        break;
        
      case OPERATOR_DIV:
        
        // DIV can be applied to vectors and tensors and has rank -1 in 2D and 3D
        if(fieldRank > 0) {
          operatorRank = -1; 
        }
        
        // DIV cannot be applied to scalar fields except in 1D where it defaults to d/dx
        else {
          TEST_FOR_EXCEPTION( ( (spaceDim > 1) && (fieldRank == 0) ),
                              std::invalid_argument,
                              ">>> ERROR (getOperatorRank): DIV cannot be applied to scalar fields in 2D and 3D");  
        }
        break;
        
      default:
        TEST_FOR_EXCEPTION( ( (operatorType != OPERATOR_VALUE) &&
                              (operatorType != OPERATOR_GRAD)  && 
                              (operatorType != OPERATOR_CURL)  && 
                              (operatorType != OPERATOR_DIV)   &&
                              (operatorType != OPERATOR_D1)    && 
                              (operatorType != OPERATOR_D2)    && 
                              (operatorType != OPERATOR_D3)    && 
                              (operatorType != OPERATOR_D4)    && 
                              (operatorType != OPERATOR_D5)    && 
                              (operatorType != OPERATOR_D6)    && 
                              (operatorType != OPERATOR_D7)    && 
                              (operatorType != OPERATOR_D8)    && 
                              (operatorType != OPERATOR_D9)    && 
                              (operatorType != OPERATOR_D10) ),
                            std::invalid_argument,
                            ">>> ERROR (getOperatorRank): Invalid operator type");
    } // switch
  }// 2D and 3D
  
  return operatorRank; 
}



int getDkEnumeration(const int xMult,
                     const int yMult,
                     const int zMult) {
  
  if( (yMult < 0) && (zMult < 0)) {
    
#ifdef HAVE_INTREPID_DEBUG
    // We are in 1D: verify input - xMult is non-negative  and total order <= 10:
    TEST_FOR_EXCEPTION( !( (0 <= xMult) && (xMult <= INTREPID_MAX_DERIVATIVE) ),
                        std::out_of_range,
                        ">>> ERROR (getDkEnumeration): Derivative order out of range");
#endif
    
    // there's only one derivative of order xMult
    return 0;
  }
  else {
    if( zMult < 0 ) {
      
#ifdef HAVE_INTREPID_DEBUG
      // We are in 2D: verify input - xMult and yMult are non-negative and total order <= 10:
      TEST_FOR_EXCEPTION( !( (0 <= xMult) && (0 <= yMult) && 
                             ( (xMult + yMult) <= INTREPID_MAX_DERIVATIVE) ),
                          std::out_of_range,
                          ">>> ERROR (getDkEnumeration): Derivative order out of range");
#endif
      
      // enumeration is the value of yMult
      return yMult;
    }
    
    // We are in 3D: verify input - xMult, yMult and zMult are non-negative and total order <= 10:
    else {
      int order = xMult + yMult + zMult;
      
#ifdef HAVE_INTREPID_DEBUG
      // Verify input:  total order cannot exceed 10:
      TEST_FOR_EXCEPTION(  !( (0 <= xMult) && (0 <= yMult) && (0 <= zMult) && 
                              (order <= INTREPID_MAX_DERIVATIVE) ),
                           std::out_of_range,
                           ">>> ERROR (getDkEnumeration): Derivative order out of range");
#endif
      int enumeration = zMult;
      for(int i = 0; i < order - xMult + 1; i++){
        enumeration += i; 
      }
      return enumeration;
    }
  }
}



void getDkMultiplicities(Teuchos::Array<int>& partialMult,
                         const int derivativeEnum,
                         const int derivativeOrder,
                         const int spaceDim) {
  
  /* Hash table to convert enumeration of partial derivative to multiplicities of dx,dy,dz in 3D.
  Multiplicities {mx,my,mz} are arranged lexicographically in bins numbered from 0 to 10. 
  The size of bins is an arithmetic progression, i.e., 1,2,3,4,5,...,11. Conversion formula is:
  \verbatim
  mx = derivativeOrder - binNumber
  mz = derivativeEnum  - binBegin
  my = derivativeOrder - mx - mz = binNumber + binBegin - derivativeEnum
  \endverbatim
  where binBegin is the enumeration of the first element in the bin. Bin numbers and binBegin 
  values are stored in hash tables for quick access by derivative enumeration value. 
  */ 
  
  // Returns the bin number for the specified derivative enumeration
  static const int binNumber[66] = { 
    0,
    1, 1,
    2, 2, 2,
    3, 3, 3, 3,
    4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5,
    6, 6, 6, 6, 6, 6, 6,
    7, 7, 7, 7, 7, 7, 7, 7,
    8, 8, 8, 8, 8, 8, 8, 8, 8,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    10,10,10,10,10,10,10,10,10,10,10
  };
  
  // Returns the binBegin value for the specified derivative enumeration 
  static const int binBegin[66] ={ 
    0,
    1, 1,
    3, 3 ,3,
    6, 6, 6, 6,
    10,10,10,10,10,
    15,15,15,15,15,15,
    21,21,21,21,21,21,21,
    28,28,28,28,28,28,28,28,
    36,36,36,36,36,36,36,36,36,
    45,45,45,45,45,45,45,45,45,45,
    55,55,55,55,55,55,55,55,55,55,55
  };
  
  // Verify derivativeOrder and derivativeEnum values, spaceDim is verified in switch statement
#ifdef HAVE_INTREPID_DEBUG
  
  // Derivative order must be between 1 and 10
  TEST_FOR_EXCEPTION( !( (1 <= derivativeOrder) && (derivativeOrder <= INTREPID_MAX_DERIVATIVE) ),
                      std::invalid_argument,
                      ">>> ERROR (getDkMultiplicities): Invalid derivative order");
  
  // Enumeration value must be between 0 and the cardinality of the derivative set
  TEST_FOR_EXCEPTION( !( (0 <= derivativeEnum) && (derivativeEnum < getDkCardinality(derivativeOrder,spaceDim) ) ),
                      std::invalid_argument,
                      ">>> ERROR (getDkMultiplicities): Invalid derivative enumeration value for this order and space dimension");
#endif
  switch(spaceDim) {
    
    case 1:
      
      // Resize return array for multiplicity of {dx}
      partialMult.resize(1);
      
      // Multiplicity of dx equals derivativeOrder
      partialMult[0] = derivativeOrder;
      break;
      
    case 2:
      
      // Resize array for multiplicities of {dx,dy}
      partialMult.resize(2);
      
      // Multiplicity of dy equals the enumeration of the derivative; of dx - the complement
      partialMult[1] = derivativeEnum;
      partialMult[0] = derivativeOrder - derivativeEnum;
      break;
      
    case 3:
      
      // Resize array for multiplicities of {dx,dy,dz}
      partialMult.resize(3);
      
      // Recover multiplicities
      partialMult[0] = derivativeOrder - binNumber[derivativeEnum];
      partialMult[1] = binNumber[derivativeEnum] + binBegin[derivativeEnum] - derivativeEnum;
      partialMult[2] = derivativeEnum  -  binBegin[derivativeEnum];
      break;
      
    default:
      TEST_FOR_EXCEPTION( !( (0 < spaceDim ) && (spaceDim < 4) ),
                          std::invalid_argument,
                          ">>> ERROR (getDkMultiplicities): Invalid space dimension");          
  }
}



int getDkCardinality(const int derivativeOrder,
                     const int spaceDim) {
  
  // Verify derivativeOrder value, spaceDim is verified in switch statement
#ifdef HAVE_INTREPID_DEBUG
  
  // Derivative order must be between 1 and 10
  TEST_FOR_EXCEPTION( !( (1 <= derivativeOrder) && (derivativeOrder <= INTREPID_MAX_DERIVATIVE) ),
                      std::invalid_argument,
                      ">>> ERROR (getDkCardinality): Invalid derivative order");
#endif
  int cardinality = -999;
  switch(spaceDim) {
    
    case 1:
      cardinality = 1;
      break;
      
    case 2:
      cardinality = derivativeOrder + 1;
      break;
      
    case 3:
      cardinality = (derivativeOrder + 1)*(derivativeOrder + 2)/2;
      break;
      
    default:
      TEST_FOR_EXCEPTION( !( (0 < spaceDim ) && (spaceDim < 4) ),
                          std::invalid_argument,
                          ">>> ERROR (getDkcardinality): Invalid space dimension");          
  }
  return cardinality;
}


// End member, friend, and related function definitions of class FieldContainer.

} // end namespace Intrepid
