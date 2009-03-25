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
  data_.assign(dim0_, 0);
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
  data_.assign(dim0_*dim1_, 0);
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
  data_.assign(dim0_*dim1_*dim2_, 0);
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
  data_.assign(dim0_*dim1_*dim2_*dim3_, 0);
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
  data_.assign(dim0_*dim1_*dim2_*dim3_*dim4_, 0);
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
  data_.assign( this -> size(), 0);
  
}



template<class Scalar>
FieldContainer<Scalar>::FieldContainer(const Teuchos::Array<int>&    dimensions,
                                       const Teuchos::Array<Scalar>& data) {
 
  // Copy all dimensions
  dimensions_.assign(dimensions.begin(),dimensions.end());
  
  // Copy first 5 dimensions to optimize class for low rank containers
  unsigned int rank = dimensions_.size();
  switch(rank) {
    case 0:
      dim0_ = 0; 
      dim1_ = 0;
      dim2_ = 0;
      dim3_ = 0;
      dim4_ = 0;
      break;
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
  TEST_FOR_EXCEPTION( ( (int)data.size() != this -> size() ),
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
inline int FieldContainer<Scalar>::rank() const {
  return dimensions_.size();
}
  


template<class Scalar>
int FieldContainer<Scalar>::size() const {
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
template<class Vector>
inline void FieldContainer<Scalar>::dimensions(Vector& dimensions) const {
  dimensions.assign(dimensions_.begin(),dimensions_.end());
}



template<class Scalar>
inline int FieldContainer<Scalar>::dimension(const int whichDim) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (0 > whichDim), std::invalid_argument,
                      ">>> ERROR (FieldContainer): dimension order cannot be negative");
  TEST_FOR_EXCEPTION( (whichDim >= this -> rank() ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): dimension order cannot exceed rank of the container");
#endif
  return dimensions_[whichDim];
}



template<class Scalar>
inline int FieldContainer<Scalar>::getEnumeration(const int i0) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> rank() != 1), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 2), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 3), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 4), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 5), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 1), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 2), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 3), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 4), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 5), std::invalid_argument, 
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
template<class Vector>
void FieldContainer<Scalar>::getMultiIndex(Vector &             multiIndex,
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
inline void FieldContainer<Scalar>::clear() {
  dimensions_.resize(0);
  
  // Reset first five dimensions:
  dim0_ = 0;
  dim1_ = 0;
  dim2_ = 0;
  dim3_ = 0;
  dim4_ = 0;
  
  // Clears data array and sets to zero length
  data_.clear();
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
    data_.resize(this -> size());
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
  
  // Copy dimensions from the specified container
  anotherContainer.dimensions(dimensions_);
  int newRank = dimensions_.size();
  
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
  data_.resize(this->size());  
}


template<class Scalar>
void FieldContainer<Scalar>::resize(const int             numPoints,
                                    const int             numFields,
                                    const EFunctionSpace  spaceType,
                                    const EOperator       operatorType,
                                    const int             spaceDim) {  
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
  int fieldRank    = getFieldRank(spaceType);
  int operatorRank = getOperatorRank(spaceType,operatorType,spaceDim);
  
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
      newDimensions[2 + fieldRank] = getDkCardinality(operatorType,spaceDim);
      break;
      
    default:
      TEST_FOR_EXCEPTION( !(Intrepid::isValidOperator(operatorType) ), std::invalid_argument,
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
inline void FieldContainer<Scalar>::initialize(const Scalar value) {
  for (int i=0; i < this->size(); i++) {
    data_[i] = value;
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
                      ">>> ERROR (FieldContainer): Size of argument does not match the size of container.");  
#endif  
  data_.assign(dataArray.begin(),dataArray.end());
}



template<class Scalar>
void FieldContainer<Scalar>::setValues(const Scalar* dataPtr, 
                                       const int numData) 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (numData != this -> size() ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Number of data does not match the size of container.");  

#endif
  data_.assign(dataPtr, dataPtr + numData);  
}


template<class Scalar>
template<class ArrayType>
void FieldContainer<Scalar>::contractScalar(ArrayType &                     outputValues,
                                            const FieldContainer<Scalar> &  rightValues,
                                            const ECompEngine               compEngine) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (this -> dimensions_.size() != 3 ), std::invalid_argument,
                       ">>> ERROR (FieldContainer): Rank of the calling object must equal 3");
  TEST_FOR_EXCEPTION( (rightValues.rank() != 3 ), std::invalid_argument,
                     ">>> ERROR (FieldContainer): Rank of input container must equal 3!");
  TEST_FOR_EXCEPTION( (dim0_ != rightValues.dimension(0) ), std::invalid_argument,
                     ">>> ERROR (FieldContainer): First dimensions (number of integration domains) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim1_ != rightValues.dimension(1) ), std::invalid_argument,
                     ">>> ERROR (FieldContainer): Second dimensions (numbers of integration points) of the calling and input containers must agree!");
#endif
  /*
   Contracts the "point" dimension P of two rank-3 containers with dimensions (C,P,L) and (C,P,R), resp.
   and returns the result in a rank-3 container with dimensions (C,L,R). The "left" container is
   the calling FieldContainer object; the "right" container is the input argument. For a fixed 
   dimension "C", (C,L,R) represents a rectangular L X R matrix where L and R may be different.
        C - num. integration domains       dim0_ in both containers
        P - num. integration points        dim1_ in both containers
        L - num. "left" fields             dim2_ in "left" (calling) container
        R - num. "right" fields            dim2_ in "right" (input) container
  */
  // resize output container
  int numRightBfs  = rightValues.dimension(2);
  outputValues.resize(dim0_, dim2_, numRightBfs);
  
  switch(compEngine) {
    case COMP_CPP: {
      for (int cl = 0; cl < dim0_; cl++) {
        for (int lbf = 0; lbf < dim2_; lbf++) {          
          for (int rbf = 0; rbf < numRightBfs; rbf++) {
            Scalar tmpVal(0);
            for (int qp = 0; qp < dim1_; qp++) {
              tmpVal += (*this)(cl, qp, lbf)*rightValues(cl, qp, rbf);
            } // P-loop
            outputValues(cl, lbf, rbf) = tmpVal;
          } // R-loop
        } // L-loop
      } // C-loop      
    }
      break;
      
    case COMP_BLAS: {
      /*
       GEMM parameters and their values:
       TRANSA   NO_TRANS
       TRANSB   TRANS
       M        #rows(A)                              = dim2_ of calling (left) container
       N        #cols(B^T)             = numRightBfs  = dim2_ of input (right) container
       K        #cols(A)                              = dim1_ of both containers
       ALPHA    1.0
       A        left data for cell cl  = &this->.getData()[cl*skipL]
       LDA      #rows(A)                              = dim2_ of calling (left) container
       B        right data for cell cl = &rightValues.getData()[cl*skipR]
       LDB      #rows(B)               = numRightBfs  = dim2_ of input (right) container
       BETA     0.0
       C        result for cell cl     = outputValues.getData()[cl*skipOp]
       LDC      #rows(C)                              = dim2_ of input (right) container
      */
      int skipL    = dim2_*dim1_;               // size of the left data chunk per cell
      int skipR    = numRightBfs*dim1_;         // size of the right data chunk per cell
      int skipOp   = dim2_*numRightBfs;         // size of the output data chunk per cell
      double alpha = 1.0;                       // these are left unchanged by GEMM 
      double beta  = 0.0;
      
      for (int cl=0; cl < dim0_; cl++) {
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::NO_TRANS,
                    Teuchos::TRANS,
                    dim2_,
                    numRightBfs,
                    dim1_,
                    alpha,
                    &this -> getData()[cl*skipL],
                    dim2_,
                    &rightValues.getData()[cl*skipR],
                    numRightBfs,
                    beta,
                    &outputValues.getData()[cl*skipOp],
                    dim2_);
      }
    }
      break;
      
    default:
      TEST_FOR_EXCEPTION( ( (compEngine != COMP_CPP) && (compEngine != COMP_BLAS) ), std::invalid_argument,
                         ">>> ERROR (FieldContainer): Computational engine not defined!");
  } // switch(compEngine)  
} // contractScalar


template<class Scalar>
template<class ArrayType>
void FieldContainer<Scalar>::contractVector(ArrayType &                     outputValues,
                                            const FieldContainer<Scalar> &  rightValues,
                                            const ECompEngine               compEngine) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (this -> dimensions_.size() != 4 ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Rank of the calling object must equal 4");
  TEST_FOR_EXCEPTION( (rightValues.rank() != 4 ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Rank of input container must equal 4!");
  TEST_FOR_EXCEPTION( (dim0_ != rightValues.dimension(0) ), std::invalid_argument,
                     ">>> ERROR (FieldContainer): First dimensions (number of integration domains) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim1_ != rightValues.dimension(1) ), std::invalid_argument,
                     ">>> ERROR (FieldContainer): Second dimensions (numbers of integration points) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim2_ != rightValues.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Third dimensions (vector length) of the calling and input containers must agree!");
#endif
  /*
   Contracts the "point" and "space" dimensions P and D1 of two rank-4 containers with 
   dimensions (C,P,D1,L) and (C,P,D1,R), resp. and returns the result in a rank-3 container 
   with dimensions (C,L,R). The "left" container is the calling FieldContainer object; the 
   "right" container is the input argument. For a fixed index "C", (C,L,R) represents 
   a rectangular L x R matrix where L and R may be different.
        C - num. integration domains       dim0_ in both containers
        P - num. integration points        dim1_ in both containers
        D1- vector dimension               dim2_ in both containers
        L - num. "left" fields             dim3_ in "left" (calling) container
        R - num. "right" fields            dim3_ in "right" (input) container
   */
  // resize output container
  int numRightBfs  = rightValues.dimension(3);
  outputValues.resize(dim0_, dim3_, numRightBfs);
  
  switch(compEngine) {
    case COMP_CPP: {
      for (int cl = 0; cl < dim0_; cl++) {
        for (int lbf = 0; lbf < dim3_; lbf++) {
          for (int rbf = 0; rbf < numRightBfs; rbf++) {
            Scalar tmpVal(0);
            for (int qp = 0; qp < dim1_; qp++) {
              for (int iVec = 0; iVec < dim2_; iVec++) {
                tmpVal += (*this)(cl, qp, iVec, lbf)*rightValues(cl, qp, iVec, rbf);
              } //D-loop
            } // P-loop
            outputValues(cl, lbf, rbf) = tmpVal;
          } // R-loop
        } // L-loop
      } // C-loop      
    }
      break;
      
    case COMP_BLAS: {
      /*
       GEMM parameters and their values:
       TRANSA   NO_TRANS
       TRANSB   TRANS
       M        #rows(A)                              = dim3_ of calling (left) container
       N        #cols(B^T)             = numRightBfs  = dim3_ of input (right) container
       K        #cols(A)               = numData      = dim1_*dim2_
       ALPHA    1.0
       A        left data for cell cl  = &this->.getData()[cl*skipL]
       LDA      #rows(A)                              = dim3_ of calling (left) container
       B        right data for cell cl = &rightValues.getData()[cl*skipR]
       LDB      #rows(B)               = numRightBfs  = dim3_ of input (right) container
       BETA     0.0
       C        result for cell cl     = outputValues.getData()[cl*skipOp]
       LDC      #rows(C)                              = dim3_ of input (right) container
       numData = num. points * vector field dimension = dim1_ * dim2_
      */
      int numData  = dim1_*dim2_;       
      int skipL    = dim3_*numData;             // size of the left data chunk per cell
      int skipR    = numRightBfs*numData;       // size of the right data chunk per cell
      int skipOp   = numRightBfs*dim3_;         // size of the output data chunk per cell
      double alpha = 1.0;                       // these are left unchanged by GEMM 
      double beta  = 0.0;
      
      for (int cl = 0; cl < dim0_; cl++) {
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::NO_TRANS,
                    Teuchos::TRANS,
                    dim3_,
                    numRightBfs,
                    numData,
                    alpha,
                    &this -> getData()[cl*skipL],
                    dim3_,
                    &rightValues.getData()[cl*skipR],
                    numRightBfs,
                    beta,
                    &outputValues.getData()[cl*skipOp],
                    dim3_);
      }
    }
      break;
      
    default:
      TEST_FOR_EXCEPTION( ( (compEngine != COMP_CPP) && (compEngine != COMP_BLAS) ), std::invalid_argument,
                         ">>> ERROR (FieldContainer): Computational engine not defined!");
  } // switch(compEngine)  
} // contractVector



template<class Scalar>
template<class ArrayType>
void FieldContainer<Scalar>::contractTensor(ArrayType &                     outputValues,
                                            const FieldContainer<Scalar> &  rightValues,
                                            ECompEngine                     compEngine) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (this -> dimensions_.size() != 5 ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Rank of the calling object must equal 5");
  TEST_FOR_EXCEPTION( (rightValues.rank() != 5 ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Rank of input container must equal 5!");
  TEST_FOR_EXCEPTION( (dim0_ != rightValues.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): First dimensions (number of integration domains) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim1_ != rightValues.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Second dimensions (numbers of integration points) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim2_ != rightValues.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Third dimensions (1st tensor dim) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim3_ != rightValues.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Fourth dimensions (2nd tensor dim) of the calling and input containers must agree!");
#endif
  /*
   Contracts the "point" and "space" dimensions P, D1 and D2 of two rank-5 containers with
   dimensions (C,P,D1,D2,L) and (C,P,D1,D2,R), resp. and returns the result in a rank-3 
   container with dimensions (C,L,R).  The "left" container is the calling FieldContainer 
   object; the "right" container is the input argument. For a fixed index "C", (C,L,R) 
   represents a rectangular L x R matrix where L and R may be different.
        C - num. integration domains       dim0_ in both containers
        P - num. integration points        dim1_ in both containers
        D1- 1st tensor dimension           dim2_ in both containers
        D2- 2nd tensor dimension           dim3_ in both containers
        L - num. "left" fields             dim4_ in "left" (calling) container
        R - num. "right" fields            dim4_ in "right" (input) container
  */
  // resize output container
  int numRightBfs  = rightValues.dimension(4);
  outputValues.resize(dim0_, dim4_, numRightBfs);
  
  switch(compEngine) {
    case COMP_CPP: {
      for (int cl = 0; cl < dim0_; cl++) {
        for (int lbf = 0; lbf < dim4_; lbf++) {
          for (int rbf = 0; rbf < numRightBfs; rbf++) {
            Scalar tmpVal(0);
            for (int qp = 0; qp < dim1_; qp++) {
              for (int iTens1 = 0; iTens1 < dim2_; iTens1++) {
                for (int iTens2 =0; iTens2 < dim3_; iTens2++) {
                  tmpVal += (*this)(cl, qp, iTens1, iTens2, lbf)*rightValues(cl, qp, iTens1, iTens2, rbf);
                } // D2-loop
              } // D1-loop
            } // P-loop
            outputValues(cl, lbf, rbf) = tmpVal;
          } // R-loop
        } // L-loop
      } // C-loop      
    }
      break;
      
    case COMP_BLAS: {
      /*
       GEMM parameters and their values:
       TRANSA   NO_TRANS
       TRANSB   TRANS
       M        #rows(A)                              = dim4_ of calling (left) container
       N        #cols(B^T)             = numRightBfs  = dim4_ of input (right) container
       K        #cols(A)               = numData      = dim1_*dim2_*dim3_
       ALPHA    1.0
       A        left data for cell cl  = &this->.getData()[cl*skipL]
       LDA      #rows(A)                              = dim4_ of calling (left) container
       B        right data for cell cl = &rightValues.getData()[cl*skipR]
       LDB      #rows(B)               = numRightBfs  = dim4_ of input (right) container
       BETA     0.0
       C        result for cell cl     = outputValues.getData()[cl*skipOp]
       LDC      #rows(C)                              = dim4_ of input (right) container
       numData = num. points * num. tensor field components = dim1_ * (dim2_ * dim3_)
      */
      int numData  = dim1_*dim2_*dim3_;       
      int skipL    = dim4_*numData;             // size of the left data chunk per cell
      int skipR    = numRightBfs*numData;       // size of the right data chunk per cell
      int skipOp   = numRightBfs*dim4_;         // size of the output data chunk per cell
      double alpha = 1.0;                       // these are left unchanged by GEMM 
      double beta  = 0.0;
      
      for (int cl=0; cl < dim0_; cl++) {
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::NO_TRANS,
                    Teuchos::TRANS,
                    dim4_,
                    numRightBfs,
                    numData,
                    alpha,
                    &this -> getData()[cl*skipL],
                    dim4_,
                    &rightValues.getData()[cl*skipR],
                    numRightBfs,
                    beta,
                    &outputValues.getData()[cl*skipOp],
                    dim4_);
      }
    }
      break;
      
    default:
      TEST_FOR_EXCEPTION( ( (compEngine != COMP_CPP) && (compEngine != COMP_BLAS) ), std::invalid_argument,
                         ">>> ERROR (FieldContainer): Computational engine not defined!");
  } // switch(compEngine)  
} // contractTensor



template<class Scalar>
template<class ArrayType>
void FieldContainer<Scalar>::contractScalarData(ArrayType &         outputValues,
                                                const ArrayType &   inputData,
                                                const ECompEngine   compEngine) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (this -> dimensions_.size() != 3 ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Rank of the calling object must equal 3");
  TEST_FOR_EXCEPTION( (inputData.rank() != 2 ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Rank of data array must equal 2!");
  TEST_FOR_EXCEPTION( (dim0_ != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): First dimensions (number of integration domains) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim1_ != inputData.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Second dimensions (numbers of integration points) of the calling and input containers must agree!");
#endif
  /*
   Contracts the "point" dimensions of a rank-3 (calling) and  rank-2 (input) containers with dimensions 
   (C,P,F) and (C,P), respectively and returns a rank-2 container with dimensions (C,F). For a fixed index "C",
   (C,F) represents a (column) vector of length F.  The "left" container is the calling FieldContainer 
   object; the "right" container is the input data array.  
        C - num. integration domains       dim0_ in both containers
        P - num. integration points        dim1_ in both containers
        F - num. "left" fields             dim2_ in "left" (calling) container
  */
  // resize output container
  outputValues.resize(dim0_, dim2_);
  
  switch(compEngine) {
    case COMP_CPP: {
      for (int cl = 0; cl < dim0_; cl++) {
        for (int lbf = 0; lbf < dim2_; lbf++) {
          Scalar tmpVal(0);
          for (int qp = 0; qp < dim1_; qp++) {
            tmpVal += (*this)(cl, qp, lbf)*inputData(cl, qp);
          } // P-loop
          outputValues(cl, lbf) = tmpVal;
        } // L-loop
      } // C-loop      
    }
      break;
      
    case COMP_BLAS: {
      /*
       GEMM parameters and their values:
       TRANSA   NO_TRANS
       TRANSB   TRANS
       M        #rows(A)               = dim2_ of calling container
       N        #cols(B^T)             = 1
       K        #cols(A)               = dim1_ of both containers (num points - contraction dimension)
       ALPHA    1.0
       A        left data for cell cl  = &this -> getData()[cl*skipL]
       LDA      #rows(A)               = dim2_
       B        right data for cell cl = &inputData.getData()[cl*skipR]
       LDB      #rows(B)               = 1
       BETA     0.0
       C        result for cell cl     = outputValues.getData()[cl*skipOp]
       LDC      #rows(C)               = dim2_ of caling container
      */
      int skipL    = dim2_*dim1_;       // size of the left data chunk per cell
      int skipR    = dim1_;             // size of the right data chunk per cell
      int skipOp   = dim2_;             // size of the output data chunk per cell
      double alpha = 1.0;               // these are left unchanged by GEMM 
      double beta  = 0.0;
      
      for (int cl=0; cl < dim0_; cl++) {
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::NO_TRANS,
                    Teuchos::TRANS,
                    dim2_,
                    1,
                    dim1_,
                    alpha,
                    &this -> getData()[cl*skipL],
                    dim2_,
                    &inputData.getData()[cl*skipR],
                    1,
                    beta,
                    &outputValues.getData()[cl*skipOp],
                    dim2_);
      }
    }
      break;
      
    default:
      TEST_FOR_EXCEPTION( ( (compEngine != COMP_CPP) && (compEngine != COMP_BLAS) ), std::invalid_argument,
                         ">>> ERROR (FieldContainer): Computational engine not defined!");
  } // switch(compEngine)  
} // contractScalarData



template<class Scalar>
template<class ArrayType>
void FieldContainer<Scalar>::contractVectorData(ArrayType &        outputValues,
                                                const ArrayType &  inputData,
                                                const ECompEngine  compEngine) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (this -> dimensions_.size() != 4 ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Rank of the calling object must equal 4");
  TEST_FOR_EXCEPTION( (inputData.rank() != 3 ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Rank of data array must equal 3!");
  TEST_FOR_EXCEPTION( (dim0_ != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): First dimensions (number of integration domains) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim1_ != inputData.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Second dimensions (numbers of integration points) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim2_ != inputData.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Third dimensions (vector length) of the calling and input containers must agree!");
#endif
  /*
   Contracts the "point" and "space" dimensions P and D1 of rank-4 (calling) and rank-3 (input) containers 
   with dimensions (C,P,D1,F) and (C,P,D1), resp. and returns the result in a rank-2 container with 
   dimensions (C,F). The "left" container is the calling FieldContainer object; the "right" container 
   is the input argument. For a fixed index "C", (C,F) represents a (row) vector of length F. 
        C - num. integration domains       dim0_ in both containers
        P - num. integration points        dim1_ in both containers
        D1- vector dimension               dim2_ in both containers
        F - num. "left" fields             dim3_ in "left" (the calling) container
  */
  // resize output container
  outputValues.resize(dim0_, dim3_);
  
  switch(compEngine) {
    case COMP_CPP: {
      for (int cl = 0; cl < dim0_; cl++) {
          for (int lbf = 0; lbf < dim3_; lbf++) {
            Scalar tmpVal(0);
            for (int qp = 0; qp < dim1_; qp++) {
              for (int iVec = 0; iVec < dim2_; iVec++) {
                tmpVal += (*this)(cl, qp, iVec, lbf)*inputData(cl, qp, iVec);
              } //D-loop
            } // P-loop
            outputValues(cl, lbf) = tmpVal;
          } // L-loop
      } // C-loop      
    }
      break;
      
    case COMP_BLAS: {
      /*
       GEMM parameters and their values:
       TRANSA   NO_TRANS
       TRANSB   TRANS
       M        #rows(A)               = dim3_ of calling container
       N        #cols(B^T)             = 1
       K        #cols(A)               = numData        = dim1_*dim2_
       ALPHA    1.0
       A        left data for cell cl  = &this -> getData()[cl*skipL]
       LDA      #rows(A)               = dim3_
       B        right data for cell cl = &inputData.getData()[cl*skipR]
       LDB      #rows(B)               = 1
       BETA     0.0
       C        result for cell cl     = outputValues.getData()[cl*skipOp]
       LDC      #rows(C)               = dim3_ of caling container
       numData = num. points * vector field dimension = dim1_ * dim2_
      */
      int numData  = dim1_*dim2_;       
      int skipL    = numData*dim3_;       // size of the left data chunk per cell
      int skipR    = numData;             // size of the right data chunk per cell
      int skipOp   = dim3_;               // size of the output data chunk per cell
      double alpha = 1.0;                 // these are left unchanged by GEMM 
      double beta  = 0.0;
      
      for (int cl = 0; cl < dim0_; cl++) {
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::NO_TRANS,
                    Teuchos::TRANS,
                    dim3_,
                    1,
                    numData,
                    alpha,
                    &this -> getData()[cl*skipL],
                    dim3_,
                    &inputData.getData()[cl*skipR],
                    1,
                    beta,
                    &outputValues.getData()[cl*skipOp],
                    dim3_);
      }
    }
      break;
      
    default:
      TEST_FOR_EXCEPTION( ( (compEngine != COMP_CPP) && (compEngine != COMP_BLAS) ), std::invalid_argument,
                          ">>> ERROR (FieldContainer): Computational engine not defined!");
  } // switch(compEngine)  
} // contractVectorData



template<class Scalar>
template<class ArrayType>
void FieldContainer<Scalar>::contractTensorData(ArrayType &                     outputValues,
                                                const ArrayType &               inputData,
                                                ECompEngine                     compEngine) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (this -> dimensions_.size() != 5 ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Rank of the calling object must equal 5");
  TEST_FOR_EXCEPTION( (inputData.rank() != 4 ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Rank of data array must equal 4!");
  TEST_FOR_EXCEPTION( (dim0_ != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): First dimensions (number of integration domains) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim1_ != inputData.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Second dimensions (numbers of integration points) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim2_ != inputData.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Third dimensions (1st tensor dim) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim3_ != inputData.dimension(3) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Fourth dimensions (2nd tensor dim) of the calling and input containers must agree!");
#endif
  /*
   Contracts the "point" and "space" dimensions P, D1 and D2 of a rank-5 (calling) and rank-5 (input)  
   containers with dimensions (C,P,D1,D2,F) and (C,P,D1,D2), resp. and returns the result in a rank-2 
   container with dimensions (C,F). The "left" container is the calling FieldContainer object; the "right" 
   container is the input argument. For a fixed index "C", (C,F) represents a (row) vector of length F. 
        C - num. integration domains       dim0_ in both containers
        P - num. integration points        dim1_ in both containers
        D1- 1st tensor dimension           dim2_ in both containers
        D2- 2nd tensor dimension           dim3_ in both containers
        F - num. "left" fields             dim4_ in "left" (calling) container
   */
  // resize output container
  outputValues.resize(dim0_, dim4_);
  
  switch(compEngine) {
    case COMP_CPP: {
      for (int cl = 0; cl < dim0_; cl++) {
          for (int lbf = 0; lbf < dim4_; lbf++) {
            Scalar tmpVal(0);
            for (int qp = 0; qp < dim1_; qp++) {
              for (int iTens1 = 0; iTens1 < dim2_; iTens1++) {
                for (int iTens2 =0; iTens2 < dim3_; iTens2++) {
                  tmpVal += (*this)(cl, qp, iTens1, iTens2, lbf)*inputData(cl, qp, iTens1, iTens2);
                } // D2-loop
              } // D1-loop
            } // P-loop
            outputValues(cl, lbf) = tmpVal;
          } // L-loop
      } // C-loop      
    }
      break;
      
    case COMP_BLAS: {
      /*
       GEMM parameters and their values:
       TRANSA   NO_TRANS
       TRANSB   TRANS
       M        #rows(A)               = dim4_ of calling container
       N        #cols(B^T)             = 1
       K        #cols(A)               = numData        = dim1_*dim2_*dim3_
       ALPHA    1.0
       A        left data for cell cl  = &this -> getData()[cl*skipL]
       LDA      #rows(A)               = dim4_
       B        right data for cell cl = &inputData.getData()[cl*skipR]
       LDB      #rows(B)               = 1
       BETA     0.0
       C        result for cell cl     = outputValues.getData()[cl*skipOp]
       LDC      #rows(C)               = dim4_ of caling container
       numData = num. points * num. tensor field components = dim1_ * (dim2_ * dim3_)
      */
      int numData  = dim1_*dim2_*dim3_;       
      int skipL    = numData*dim4_;             // size of the left data chunk per cell
      int skipR    = numData;                   // size of the right data chunk per cell
      int skipOp   = dim4_;                     // size of the output data chunk per cell
      double alpha = 1.0;                       // these are left unchanged by GEMM 
      double beta  = 0.0;
      
      for (int cl=0; cl < dim0_; cl++) {
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::NO_TRANS,
                    Teuchos::TRANS,
                    dim4_,
                    1,
                    numData,
                    alpha,
                    &this -> getData()[cl*skipL],
                    dim4_,
                    &inputData.getData()[cl*skipR],
                    1,
                    beta,
                    &outputValues.getData()[cl*skipOp],
                    dim4_);
      }
    }
      break;
      
    default:
      TEST_FOR_EXCEPTION( ( (compEngine != COMP_CPP) && (compEngine != COMP_BLAS) ), std::invalid_argument,
                          ">>> ERROR (FieldContainer): Computational engine not defined!");
  } // switch(compEngine)  
} // contractTensorData



template<class Scalar>
template<class ArrayType>
void FieldContainer<Scalar>::multiplyScalarData(const ArrayType &  inputData)
{
  /*
   Multiplies the calling rank-3,4, or 5 FieldContainer with dimensions (C,P,F), (C,P,D1,F) or (C,P,D1,D2,F),
   representing values of a scalar, vector or a tensor set of fields, by the values in a user-specified  
   rank-2 container (C,P) representing values of a single scalar field.
        C - num. integration domains       dim0_ in both containers
        P - num. integration points        dim1_ in both containers
        Di- vector/tensor dimension        none,  dim2_,  (dim2_,dim3_)  in the calling container
        F - number of fields               dim2_, dim3_,  dim4_          in the calling container
  */
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputData.rank() != 2), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Data array has an invalid rank.");  
  TEST_FOR_EXCEPTION( (dim0_ != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): First dimensions (number of integration domains) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim1_ != inputData.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Second dimensions (numbers of integration points) of the calling and input containers must agree!");
#endif
  
  int myRank = this -> rank();
  Scalar tempData(0);
  for(int cl = 0; cl < dim0_; cl++) {
    for(int pt = 0; pt < dim1_; pt++) {
      tempData = inputData(cl, pt);
      
      switch(myRank) {
        case 3: {
          for(int bf = 0; bf < dim2_; bf++) {
            (*this)(cl, pt, bf) *= tempData;
          } //F-loop
        }// case 3
          break;
          
        case 4: {
          for(int bf = 0; bf < dim2_; bf++) {
            for( int iVec = 0; iVec < dim2_; iVec++) {
              (*this)(cl, pt, iVec, bf) *= tempData;
            } // D1-loop
          }// F-loop
        }// case 4
          break;
        case 5: {
          for(int bf = 0; bf < dim2_; bf++) {
            for( int iTens1 = 0; iTens1 < dim2_; iTens1++) {
              for( int iTens2 = 0; iTens2 < dim3_; iTens2++) {
                (*this)(cl, pt, iTens1, iTens2, bf) *= tempData;
              }// D2-loop
            } // D1-loop
          }// F-loop          
        }// case 5
          break;
        default:
          TEST_FOR_EXCEPTION( !( (myRank == 3) || (myRank ==4) || (myRank == 5) ), std::invalid_argument,
                              ">>> ERROR (FieldContainer): This method is defined only for rank-3,4 or 5 containers");  
      }// myRank
    } // P-loop
  }// C-loop
}// multiplyScalarData



template<class Scalar>
template<class ArrayType>
void FieldContainer<Scalar>::multiplyVectorData(FieldContainer<Scalar> outputValues, 
                                                const ArrayType &      inputData)
{
  /* 
   Contracts the "D1" dimension of the calling rank-4 or 5 FieldContainer, with dimensions 
   (C,P,D1,F) or (C,P,D1,D2,F), with a user-specified rank-3 container (C,P,D1). This operation is 
   equivalent to a left (row) vector multiplication of the vector or tensor field set in the calling container
   by the vector field in the user-specified container. The result is a container whose rank is one
   less than the rank of the calling container, i.e., (C,P,F) or (C,P,D2,F)
        C  - num. integration domains       dim0_ in both containers
        P  - num. integration points        dim1_ in both containers
        D1 - contracting dimension          dim2_ in both containers
        D2 - 2nd tensor dimension           dim3_ in the calling container, if its rank equals 5;
        F - number of fields                dim3_ for rank-4 and dim4_ for rank-5 calling containers
  */
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (inputData.rank() != 3), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Data array has an invalid rank.");  
  TEST_FOR_EXCEPTION( (dim0_ != inputData.dimension(0) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): First dimensions (number of integration domains) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim1_ != inputData.dimension(1) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Second dimensions (numbers of integration points) of the calling and input containers must agree!");
  TEST_FOR_EXCEPTION( (dim2_ != inputData.dimension(2) ), std::invalid_argument,
                      ">>> ERROR (FieldContainer): Contracting dimensions in data array and the calling object do not agree!");
#endif
  
  int myRank = this -> rank();
  Scalar temp(0);
  
  for(int cl = 0; cl < dim0_; cl++) {
    for(int pt = 0; pt < dim1_; pt++) {
      
      switch(myRank) {
        case 4: {
          
          // All but the contracting dimension (dim2_) are the same as in the calling container
          outputValues.resize(dim0_, dim1_, dim3_); 
          for(int bf = 0; bf < dim3_; bf++) {
            temp = 0.0;
            for( int iVec = 0; iVec < dim2_; iVec++) {
              temp += inputData(cl, pt, iVec)*(*this)(cl, pt, iVec, bf);
            } // D1-loop
            outputValues(cl, pt, bf) = temp;
          }// F-loop
        }// case 4
          
          break;
        case 5: {
          
          // All but the contracting dimension (dim2_) are the same as in the calling container
          outputValues.resize(dim0_, dim1_, dim3_, dim4_);
          for(int bf = 0; bf < dim4_; bf++) {
            for(int iTens2 = 0; iTens2 < dim3_; iTens2++) {
              temp = 0;
              for( int iTens1 = 0; iTens1 < dim2_; iTens1++) {
                temp += inputData(cl, pt, iTens1)*(*this)(cl, pt, iTens1, iTens2, bf);
              } // D1-loop
              outputValues(cl, pt, iTens2, bf) =  temp;
            }// D2-loop
          }// F-loop          
        }// case 5
          break;
        default:
          TEST_FOR_EXCEPTION( !( (myRank ==4) || (myRank == 5) ), std::invalid_argument,
                              ">>> ERROR (FieldContainer): This method is defined only for rank-4 or 5 containers");  
      }// myRank
    } // P-loop
  }// C-loop
}



template<class Scalar>
inline const Scalar& FieldContainer<Scalar>::operator () (const int i0) const 
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this -> rank() != 1), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 1), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 2), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 2), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 3), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 3), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 4), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 4), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 5), std::invalid_argument, 
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
  TEST_FOR_EXCEPTION( ( this -> rank() != 5), std::invalid_argument, 
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
Scalar& FieldContainer<Scalar>::operator [] (const int address) {
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
  
  int size = container.size();
  int rank = container.rank();
  Teuchos::Array<int> multiIndex(rank);
  Teuchos::Array<int> dimensions;
  container.dimensions(dimensions);
  
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
// End member, friend, and related function definitions of class FieldContainer.

} // end namespace Intrepid
