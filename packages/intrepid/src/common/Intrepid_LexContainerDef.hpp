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
  data_.resize( this -> getSize());
}



template<class Scalar>
LexContainer<Scalar>::LexContainer(const Teuchos::Array<int>& indexRange,
                                   const Teuchos::Array<Scalar>& data) {  
  
  // Copy upper bounds for indices
  indexRange_.assign(indexRange.begin(),indexRange.end());
  
  // Validate input: size of data must match size computed from upper index bounds in indexRange
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (int)data.size() != this -> getSize() ),
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
int LexContainer<Scalar>::getAddress(Teuchos::Array<int> multiIndex) const {

  // Check if number of multi-indices matches rank of the LexContainer object
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( multiIndex.size() != indexRange_.size() ),
                      std::length_error,
                      ">>> ERROR (LexContainer): Number of multi-indices does not match rank of container.");
#endif
    
  // Compute address using Horner's nested scheme: intialize address to 0th index value
  int address = multiIndex[0];
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( ( multiIndex[0] < 0) || ( multiIndex[0] >= indexRange_[0]) ),
                      std::out_of_range,
                      ">>> ERROR (LexContainer): Multi-index component out of range.");    
#endif  
  
  int rank = indexRange_.size();
  for (int r = 0; r < rank - 1; r++){
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( (multiIndex[r+1] < 0) || (multiIndex[r+1] >= indexRange_[r+1]) ),
                        std::out_of_range,
                        ">>> ERROR (LexContainer): Multi-index component out of range.");    
#endif
    
    // Add increment
    address = address*indexRange_[r+1] + multiIndex[r+1];
  }
  return address;
}



template<class Scalar>
int LexContainer<Scalar>::getAddress(int* multiIndexPtr) const {
  
  // Uses getAddress with Teuchos::Array argument to compute the address.
  int rank = this -> getRank();
  Teuchos::Array<int> multiIndexArray(rank);
  multiIndexArray.assign(multiIndexPtr,multiIndexPtr + rank);  
  return this -> getAddress(multiIndexArray);
}



template<class Scalar>
void LexContainer<Scalar>::getMultiIndex(Teuchos::Array<int>& multiIndex,
                                         const int valueAddress) const {
  
  // Verify address is in the admissible range for this LexContainer
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (valueAddress < 0) || (valueAddress >= (int)data_.size()) ),
                      std::out_of_range,
                      ">>> ERROR (LexContainer): Specified address is out of range.");    
#endif
  
  // make sure multiIndex has the right size to hold all multi-indices
  int rank = indexRange_.size();
  multiIndex.resize(rank);
  
  // Initializations
  int temp_addr = valueAddress;
  int temp_range = 1;
  
  // Compute product of all but the first upper bound
  for(int r = 1; r < rank ; r++){
    temp_range *=indexRange_[r];
  }
  
  // Index 0 is computed first using integer division
  multiIndex[0] = temp_addr/temp_range;
  
  // Indices 1 to (rank - 2) are computed next; will be skipped if rank <=2
  for(int r = 1; r < rank - 1; r++){
    temp_addr  -= multiIndex[r-1]*temp_range;
    temp_range /= indexRange_[r];
    multiIndex[r] = temp_addr/temp_range;
  }
  
  // Index (rank - 1) is computed last, skip if rank = 1 and keep if rank = 2
  if(rank > 1) {
    multiIndex[rank - 1] = temp_addr - multiIndex[rank - 2]*temp_range;
  }
                                         }



template<class Scalar>
inline Scalar LexContainer<Scalar>::getValue(const Teuchos::Array<int>& multiIndex) const {
  return data_[this -> getAddress(multiIndex)];
}



template<class Scalar>
inline void LexContainer<Scalar>::emptyContainer() {
  indexRange_.resize(0);
  data_.resize(0);
}



template<class Scalar>
inline void LexContainer<Scalar>::storeZero() {
  indexRange_.resize(1); 
  data_.resize(1);       
  indexRange_[0] = 1;
  data_[0] = (Scalar)0.0;
}



template<class Scalar>
inline void LexContainer<Scalar>::resizeContainer(const Teuchos::Array<int>& newIndexRange) {
  
  // Copy upper index bounds and resize container storage to match new upper bounds.
  indexRange_.assign(newIndexRange.begin(),newIndexRange.end());  
  data_.resize(this -> getSize());
}



template<class Scalar>
void LexContainer<Scalar>::resizeContainer(const int numPoints,
                                           const int numFields,
                                           const int fieldRank,
                                           const int operatorOrd,
                                           const int spaceDim){  
  // Validate input
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( numPoints < 0),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Number of points cannot be negative!");  
  TEST_FOR_EXCEPTION( ( numFields < 0),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Number of fields cannot be negative!");  
  TEST_FOR_EXCEPTION( !( (0 <= fieldRank ) && ( fieldRank <= 2 ) ),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Field rank can be 0,1, or 2 only.");  
  TEST_FOR_EXCEPTION( !( (0 <= operatorOrd ) && ( operatorOrd <= 10 ) ),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Total derivative order has to be between 0 and 10.");  
  TEST_FOR_EXCEPTION( !( (0 <  spaceDim ) && ( spaceDim <= 3  ) ),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Space dimension has to be between 1 and 3.");  
  TEST_FOR_EXCEPTION(  ( (1 == spaceDim ) && ( fieldRank > 0  ) ),
                       std::invalid_argument,
                       ">>> ERROR (LexContainer): Field rank cannot be greater than 0 in one-dimension");  
#endif  
  
  // Set rank of the LexContainer
  int rank = 0;
  if(spaceDim == 1) {
    
    // Indices count: numPoints (1) + numFields (1) + (numDerivatives)
    rank = 1 + 1 + operatorOrd; 
  }
  else {
    
    // Indices count: numPoints (1) + numFields (1)  + (fieldRank) + (numDerivatives)
    rank = 1 + 1 + fieldRank + operatorOrd; 
  }
  
  // Define an array for the new index ranges
  Teuchos::Array<int> newIndexRange(rank);
  
  // Load indexRange_ with upper bounds for each index
  newIndexRange[0] = numPoints;
  newIndexRange[1] = numFields;
  
  // For vector and tensor fields index range for components is bounded by the space dimension
  for(int fRank = 0; fRank < fieldRank; fRank++){
    newIndexRange[2 + fRank] = spaceDim; 
  }
  
  // For all field types index range for partial derivatives is bounded by space dimension
  for(int dOrd = 0; dOrd < operatorOrd; dOrd++) {
    newIndexRange[2 + fieldRank + dOrd] = spaceDim;
  }
  
  // Once newIndexRange has been set, reset to change container size to reflect the new index bounds
  this -> resizeContainer(newIndexRange);
}


template<class Scalar>
void LexContainer<Scalar>::shapeContainer(const int       numPoints,
                                          const int       numFields,
                                          const EField    fieldType,
                                          const EOperator operatorType,
                                          const int       spaceDim)
{
  // Validate input
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( numPoints < 0),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Number of points cannot be negative!");  
  TEST_FOR_EXCEPTION( ( numFields < 0),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Number of fields cannot be negative!");  
  TEST_FOR_EXCEPTION( !( (0 <  spaceDim ) && ( spaceDim <= 3  ) ),
                      std::invalid_argument,
                      ">>> ERROR (LexContainer): Space dimension has to be between 1 and 3.");  
#endif  
  
  // We need to shape a container that will hold the result of (operatorType) applied to (fieldType).
  // Depending on the values of these arguments and the space dimension, the resulting field may have
  // different rank and operator order. For example, operators VALUE and all total derivatives, i.e.,
  // Dk, preserve the rank and the operator order. Other, like DIV or CURL may change one or both
  // of these values. See comments below.
  
  // rank of the range field and its derivative order
  int rangeFieldRank = -1;
  int rangeOperatorOrd  = -1;
  
  // Rank of the input field and derivative order of the specified input operator
  int fieldRank = getFieldRank(fieldType);
  int operatorOrd = getOperatorOrder(operatorType);
  
  // If space dimension = 1 only scalar fields admissible. Any differential operator in 1D defaults
  // to a 1D derivative and so only the order of the operator matters. Only necessary to verify 
  // input and copy rank and derivative order data; no need to fall through separate operator cases.
  if(spaceDim == 1) {
#ifdef HAVE_INTREPID_DEBUG
    TEST_FOR_EXCEPTION( ( fieldRank > 0 ),
                        std::invalid_argument,
                        ">>> ERROR (LexContainer): Field rank cannot be greater than 0 in one-dimension");  
#endif  
    rangeFieldRank = fieldRank;
    rangeOperatorOrd = operatorOrd;
  }
  
  // We are either in 2D or 3D
  else {
    switch(operatorType){
      
      // The range container for value and total derivative operators has the same rank as the input 
      // field and order equal to the operator order. These operators can be applied to all fields.
      case OPERATOR_VALUE:
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
        
        rangeFieldRank = fieldRank;
        rangeOperatorOrd = operatorOrd;
        break;
        
      case OPERATOR_DIV:
        
        // DIV applied to VECTOR gives SCALAR => rangeFieldRank = 0 = fieldRank - 1
        // DIV applied to TENSOR gives VECTOR => rangeFieldRank = 1 = fieldRank - 1
        // and rangeOperatorOrd = 0. Verify that field rank is 1 or 2! (cannot apply DIV to scalars!)
#ifdef HAVE_INTREPID_DEBUG
        TEST_FOR_EXCEPTION( ( !( (fieldRank == 1) || ( fieldRank == 2) ) ),
                            std::invalid_argument,
                            ">>> ERROR (LexContainer): Invalid field type for OPERATOR_DIV");  
#endif  
        rangeFieldRank = fieldRank - 1;
        rangeOperatorOrd = 0;
        break;
        
      case OPERATOR_CURL:
        
        // CURL: Rank of the range field depends on space dimension and rank of input field:
        switch(spaceDim){
          case 3:
            
            // 3D CURL applied to VECTOR returns VECTOR => rangeFieldRank =  fieldRank
            // 3D CURL applied to TENSOR returns TENOSR => rangeFieldRank =  fieldRank
            // and rangeOperatorOrd = 0. Verify that field rank is 1 or 2! (cannot apply 3D CURL to scalars)
#ifdef HAVE_INTREPID_DEBUG
            TEST_FOR_EXCEPTION( ( !( (fieldRank == 1) || ( fieldRank == 2) ) ),
                                std::invalid_argument,
                                ">>> ERROR (LexContainer): Invalid field type for 3D OPERATOR_CURL");  
#endif  
            rangeFieldRank = fieldRank;
            rangeOperatorOrd = 0;
            break;
            
          case 2:
            
            // 2D CURL: VECTOR -> SCALAR => if(fieldRank == 1) { rangeFieldRank = 0}  
            // 2D CURL: SCALAR -> VECTOR => if(fieldRank == 0) { rangeFieldRank = 1}  
            // In both cases rangeOperatorOrd = 0. Verify that field rank is 0 or 1
#ifdef HAVE_INTREPID_DEBUG
            TEST_FOR_EXCEPTION( ( !( (fieldRank == 0) || ( fieldRank == 1) ) ),
                                std::invalid_argument,
                                ">>> ERROR (LexContainer): Invalid field type for 2D OPERATOR_CURL");  
#endif  
            // range rank can be computed without logic statements by the simple formula below            
            rangeFieldRank = 1 - fieldRank;
            rangeOperatorOrd = 0;
            break;
            
          default:
            TEST_FOR_EXCEPTION( !( (spaceDim == 2) || (spaceDim == 3) ),
                                std::invalid_argument,
                                ">>> ERROR (LexContainer): Invalid space dimension");
        }// spaceDim
        break; // CURL
        
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
                            ">>> ERROR (LexContainer): Invalid operator type");
    }// operatorType
  }// !(spaceDim == 1)
  
  // Call specialized resize method
  this -> resizeContainer(numPoints,
                          numFields,
                          rangeFieldRank,
                          rangeOperatorOrd,
                          spaceDim);
}



template<class Scalar>
inline void LexContainer<Scalar>::setValue(const Scalar dataValue, const Teuchos::Array<int>& multiIndex) {
  data_[this -> getAddress(multiIndex)] = dataValue; 
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
                      ">>> ERROR (Point): address out of range.");
#endif
  return data_[address];
}
  


template<class Scalar>
inline LexContainer<Scalar>& LexContainer<Scalar>::operator = (const LexContainer<Scalar>& right)
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this == &right ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid right-hand side to '='. Self-assignment prohibited.");
  TEST_FOR_EXCEPTION( ( this -> getSize() != right.getSize() ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid size of right-hand side argument to '='.");
  TEST_FOR_EXCEPTION( ( indexRange_.size() != right.indexRange_.size() ),
                      std::invalid_argument,
                      ">>> ERROR (Point): Invalid rank of right-hand side argument to '='.");
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
  
  if( (rank == 0 ) && (size == 0) ) {
    os<< "===============================================================================\n"\
    << "\t Container size = " << size << "   rank = " << rank << "\n"
    << "===============================================================================\n"\
    << "|                     *** This is an empty container ****                     |\n";
  }
  else {
    os<< "===============================================================================\n"\
    << "\t Container size = " << size << "   rank = " << rank << "\n"
    << "===============================================================================\n"\
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
  int opOrder;
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
      opOrder = -1;
  }
  return opOrder;
}    

// End member, friend, and related function definitions of class LexContainer.

} // end namespace Intrepid
