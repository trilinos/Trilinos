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

/** \file   Intrepid_LexContainer.hpp
    \brief  Definition file for utility class to provide lexicographical containers.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

//===========================================================================//
//                                                                           //
//          Member function definitions of the class LexContainer.           //
//                                                                           //
//===========================================================================//

template<class Scalar>
LexContainer<Scalar>::LexContainer(const LexContainer<Scalar>& right) {
  
  // Copy indexRange and data values from right
  indexRange_.assign(right.indexRange_.begin(),right.indexRange_.end());  
  data_.assign(right.data_.begin(),right.data_.end());
}



template<class Scalar>
LexContainer<Scalar>::LexContainer(const Teuchos::Array<int>& indexRange) {  
  
  // Copy upper index bounds and resize container storage to match them
  indexRange_ .assign(indexRange.begin(),indexRange.end());  
  data_.resize( this->getSize());
}



template<class Scalar>
LexContainer<Scalar>::LexContainer(const Teuchos::Array<int>& indexRange,
                                   const Teuchos::Array<Scalar>& data) {  
  
  // Copy upper bounds for indices
  indexRange_.assign(indexRange.begin(),indexRange.end());
  
  // Validate input: size of data must match size computed from upper index bounds in indexRange
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (int)data.size() != this->getSize() ),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Size of data does not match specified index range.");
#endif
  
  // Assign data
  data_.assign(data.begin(),data.end());
}



template<class Scalar>
inline int LexContainer<Scalar>::getRank() const {
  return indexRange_.size();
}
  


template<class Scalar>
inline int LexContainer<Scalar>::getSize() const {
  
  // Size equals product of all upper index bounds in indexRange_
  int rank = indexRange_.size();
  if(rank == 0) {
    return 0;
  }
  else {
    int size = indexRange_[0];
    for(int r = 1; r < rank ; r++){
      size *= indexRange_[r];
    }
    return size;
  }
}



template<class Scalar>
void LexContainer<Scalar>::getIndexRange(Teuchos::Array<int>& indexRange) const {
  
  //Resize return argument to match rank of the container
  indexRange.resize(this->getRank()); 
  
  indexRange.assign(indexRange_.begin(),indexRange_.end());
}



template<class Scalar>
int LexContainer<Scalar>::getIndexBound(const int whichIndex) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (0 > whichIndex),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): which index number cannot be negative");
  TEST_FOR_EXCEPTION( (whichIndex >= this->getRank() ),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): which index number cannot exceed container rank");
#endif
  return indexRange_[whichIndex];
}


template<class Scalar>
int LexContainer<Scalar>::getEnumeration(const Teuchos::Array<int>& multiIndex) const {

  // Check if empty multi-index.
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( multiIndex.size() == 0 ),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Empty multiIndex!");
#endif

  // Check if number of multi-indices matches rank of the LexContainer object
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( multiIndex.size() != indexRange_.size() ),
                      std::length_error,
                      ">>> ERROR (LexContainer): Number of multi-indices does not match rank of container.");
#endif
    
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( ( multiIndex[0] < 0) || ( multiIndex[0] >= indexRange_[0]) ),
                      std::out_of_range,
                      ">>> ERROR (LexContainer): Multi-index component out of range.");    
#endif  

  int rank = indexRange_.size();
  int address = 0;
  switch (rank) {
    // the first four cases allow for somewhat faster address computations
    case 4:
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( (multiIndex[3] < 0) || (multiIndex[3] >= indexRange_[3]) ),
                          std::out_of_range,
                          ">>> ERROR (LexContainer): Multi-index component out of range.");    
#endif
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( (multiIndex[2] < 0) || (multiIndex[2] >= indexRange_[2]) ),
                          std::out_of_range,
                          ">>> ERROR (LexContainer): Multi-index component out of range.");    
#endif
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( (multiIndex[1] < 0) || (multiIndex[1] >= indexRange_[1]) ),
                          std::out_of_range,
                          ">>> ERROR (LexContainer): Multi-index component out of range.");    
#endif
      address = ((multiIndex[0]*indexRange_[1]+multiIndex[1])*indexRange_[2]+multiIndex[2])*indexRange_[3]+multiIndex[3];
      break;

    case 3:
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( (multiIndex[2] < 0) || (multiIndex[2] >= indexRange_[2]) ),
                          std::out_of_range,
                          ">>> ERROR (LexContainer): Multi-index component out of range.");    
#endif
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( (multiIndex[1] < 0) || (multiIndex[1] >= indexRange_[1]) ),
                          std::out_of_range,
                          ">>> ERROR (LexContainer): Multi-index component out of range.");    
#endif
      address = (multiIndex[0]*indexRange_[1]+multiIndex[1])*indexRange_[2]+multiIndex[2];
      break;

    case 2:
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( ( (multiIndex[1] < 0) || (multiIndex[1] >= indexRange_[1]) ),
                          std::out_of_range,
                          ">>> ERROR (LexContainer): Multi-index component out of range.");    
#endif
      address = multiIndex[0]*indexRange_[1]+multiIndex[1];
      break;
    case 1:
      address = multiIndex[0];
      break;

    default:
      // Compute address using Horner's nested scheme: intialize address to 0th index value
      address = multiIndex[0];
      for (int r = 0; r < rank - 1; r++){
#ifdef HAVE_INTREPID_DEBUG
        TEST_FOR_EXCEPTION( ( (multiIndex[r+1] < 0) || (multiIndex[r+1] >= indexRange_[r+1]) ),
                            std::out_of_range,
                            ">>> ERROR (LexContainer): Multi-index component out of range.");    
#endif
        // Add increment
        address = address*indexRange_[r+1] + multiIndex[r+1];
      }
  } // end switch(rank)

  return address;
}



template<class Scalar>
int LexContainer<Scalar>::getEnumeration(int* multiIndexPtr) const {
  
  // Uses getEnumeration with Teuchos::Array argument to compute the address.
  int rank = this->getRank();
  Teuchos::Array<int> multiIndexArray(rank);
  multiIndexArray.assign(multiIndexPtr,multiIndexPtr + rank);  
  return this->getEnumeration(multiIndexArray);
}



template<class Scalar>
void LexContainer<Scalar>::getMultiIndex(Teuchos::Array<int>& multiIndex,
                                         const int valueEnum) const {
  
  // Verify address is in the admissible range for this LexContainer
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)data_.size()) ),
                      std::out_of_range,
                      ">>> ERROR (LexContainer): Specified address is out of range.");    
#endif
  
  // make sure multiIndex has the right size to hold all multi-indices
  int rank = indexRange_.size();
  multiIndex.resize(rank);
  
  // Initializations
  int temp_enum = valueEnum;
  int temp_range = 1;
  
  // Compute product of all but the first upper bound
  for(int r = 1; r < rank ; r++){
    temp_range *=indexRange_[r];
  }
  
  // Index 0 is computed first using integer division
  multiIndex[0] = temp_enum/temp_range;
  
  // Indices 1 to (rank - 2) are computed next; will be skipped if rank <=2
  for(int r = 1; r < rank - 1; r++){
    temp_enum  -= multiIndex[r-1]*temp_range;
    temp_range /= indexRange_[r];
    multiIndex[r] = temp_enum/temp_range;
  }
  
  // Index (rank - 1) is computed last, skip if rank = 1 and keep if rank = 2
  if(rank > 1) {
    multiIndex[rank - 1] = temp_enum - multiIndex[rank - 2]*temp_range;
  }
                                         }



template<class Scalar>
inline Scalar LexContainer<Scalar>::getValue(const Teuchos::Array<int>& multiIndex) const {
  return data_[this->getEnumeration(multiIndex)];
}



template<class Scalar>
inline void LexContainer<Scalar>::empty() {
  indexRange_.resize(0);
  data_.resize(0);
}



template<class Scalar>
inline void LexContainer<Scalar>::storeZero() {
  for (int i=0; i < this->getSize(); i++) {
    data_[i] = (Scalar)0;
  } 
}



template<class Scalar>
inline void LexContainer<Scalar>::resize(const Teuchos::Array<int>& newIndexRange) {
  
  // Copy upper index bounds and resize container storage to match new upper bounds.
  indexRange_.assign(newIndexRange.begin(),newIndexRange.end());  
  data_.resize(this->getSize());
}



template<class Scalar>
inline void LexContainer<Scalar>::resize(const LexContainer<Scalar>& anotherContainer) {
  
  int anotherRank = anotherContainer.getRank();
  Teuchos::Array<int> newIndexRange(anotherRank);
  for(int i = 0; i < anotherRank; i++){
    newIndexRange[i] = anotherContainer.getIndexBound(i);
  }

  // Copy upper index bounds and resize container storage to match new upper bounds.
  indexRange_.assign(newIndexRange.begin(),newIndexRange.end());  
  data_.resize(this->getSize());  
}



template<class Scalar>
inline void LexContainer<Scalar>::setValue(const Scalar dataValue, const Teuchos::Array<int>& multiIndex) {
  data_[this->getEnumeration(multiIndex)] = dataValue; 
}



template<class Scalar>
inline void LexContainer<Scalar>::setValue(const Scalar dataValue, const int index) {
  data_[index] = dataValue; 
}



template<class Scalar>
void LexContainer<Scalar>::setValues(const Teuchos::Array<Scalar>& dataArray) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (dataArray.size() != (data_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Size of argument does not match size of container.");  
#endif  
  data_.assign(dataArray.begin(),dataArray.end());
}



template<class Scalar>
inline void LexContainer<Scalar>::setValues(const Scalar* dataPtr, const int dataSize) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (dataSize != (int)data_.size() ),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Size of container does not match dataSize value.");  
#endif  
  data_.assign(dataPtr, dataPtr + dataSize); 
}



template<class Scalar>
const Scalar& LexContainer<Scalar>::operator [] (const int address) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (address < 0) || (address >= (int)data_.size() ) ),
                      std::out_of_range,
                      ">>> ERROR (LexContainer): Specified address is out of range.");
#endif
  return data_[address];
}
  


template<class Scalar>
inline LexContainer<Scalar>& LexContainer<Scalar>::operator = (const LexContainer<Scalar>& right)
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this == &right ),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Invalid right-hand side to '='. Self-assignment prohibited.");
#endif
  data_ = right.data_;
  indexRange_ = right.indexRange_; 
  return *this;
}


//===========================================================================//
//                                                                           //
//           END of member definitions; START friends and related            //
//                                                                           //
//===========================================================================//


template<class Scalar>
std::ostream& operator << (std::ostream& os, const LexContainer<Scalar>& container) {
  
  // Save the format state of the original ostream os.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(os);  

  os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  os.setf(std::ios_base::right);
  int myprec = os.precision();
  
  int size = container.getSize();
  int rank = container.getRank();
  Teuchos::Array<int> multiIndex(rank);
  Teuchos::Array<int> indexRange;
  container.getIndexRange(indexRange);
  
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
      os << " (" << indexRange[r] <<") ";
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

// End member, friend, and related function definitions of class LexContainer.

} // end namespace Intrepid
