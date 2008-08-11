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

/** \file   Intrepid_VarContainer.hpp
  \brief  Definition file for templated variable container to store values of sets of vector and 
  scalar fields and their derivatives, evaluated at a set of points.
  \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

//===========================================================================//
//                                                                           //
//          Member function definitions of the class VarContainer.           //
//                                                                           //
//===========================================================================//

template<class Scalar>
VarContainer<Scalar>::VarContainer(const VarContainer<Scalar>& right) {
  this->reset(right.getNumPoints(), right.getNumFields(), right.fieldType_, right.operatorType_, right.spaceDim_);
  this->data_   = right.data_;
}



template<class Scalar>
VarContainer<Scalar>::VarContainer(const int         numPoints,
                                   const int         numFields,
                                   const EField      fieldType,
                                   const EOperator   operatorType,
                                   const int         spaceDim) {  
  // Use reset method to adjust the internal storage and set private data
  this->reset(numPoints,numFields,fieldType,operatorType,spaceDim);
}



template<class Scalar>
inline int VarContainer<Scalar>::getNumPoints() const {
  // Number of points is upper bound for 1st index: stored in dataContainer_.indexRange_[0]
  int val = 0;
  if (this->getSize() != 0) {
    val = this->getIndexBound(0);
  }
  return val;
}



template<class Scalar>
inline int VarContainer<Scalar>::getNumFields() const {
  // Number of fields is upper bound for 2nd index: stored in dataContainer_.indexRange_[1]
  int val = 0;
  if (this->getRank() != 0) {
    val = this->getIndexBound(1);
  }
  return val;
}



template<class Scalar>
inline EField VarContainer<Scalar>::getFieldType() const {
  return fieldType_; 
}



template<class Scalar>
inline EOperator VarContainer<Scalar>::getOperatorType() const {
  return operatorType_; 
}



template<class Scalar>
inline int VarContainer<Scalar>::getSpaceDim() const {
  return spaceDim_;
}



template<class Scalar>
inline int VarContainer<Scalar>::getEnumeration(const Teuchos::Array<int>& multiIndex) const {  
  return this->LexContainer<Scalar>::getEnumeration(multiIndex);
}



template<class Scalar>
void VarContainer<Scalar>::reset(const int       numPoints,
                                 const int       numFields,
                                 const EField    fieldType,
                                 const EOperator operatorType,
                                 const int       spaceDim) {  
  // Validate input
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( numPoints < 0),
                      std::invalid_argument,
                      ">>> ERROR (VarContainer): Number of points cannot be negative!");  
  TEST_FOR_EXCEPTION( ( numFields < 0),
                      std::invalid_argument,
                      ">>> ERROR (VarContainer): Number of fields cannot be negative!");  
  TEST_FOR_EXCEPTION( !( (1 <=  spaceDim ) && ( spaceDim <= 3  ) ),
                      std::invalid_argument,
                      ">>> ERROR (VarContainer): Space dimension has to be between 1 and 3.");  
#endif  
  
  // Find out field and operator ranks
  int fieldRank    = getFieldRank(fieldType);
  int operatorRank = getOperatorRank(operatorType,fieldRank,spaceDim);
  
  // Compute rank of the container = 1(numPoints) + 1(numFields) + fieldRank + operatorRank
  int rank = 1 + 1 + fieldRank + operatorRank;
  
  // Define temp array for the index ranges
  Teuchos::Array<int> newIndexRange(rank);
  
  // Number of points and number of fields are upper bounds for indices 0 and 1, resp.
  newIndexRange[0] = numPoints;
  newIndexRange[1] = numFields;
  
  // The rest of the upper bounds depend on whether we had VALUE, GRAD (D1), CURL, DIV or Dk, k>1
  switch(operatorType) {

    case OPERATOR_VALUE:
    case OPERATOR_GRAD:
    case OPERATOR_D1:
    case OPERATOR_CURL:
    case OPERATOR_DIV:
      
      // For these operators all indices from 2 to 2 + fieldRank + OperatorRank are bounded by spaceDim
      for(int i = 0; i < fieldRank + operatorRank; i++){
        newIndexRange[2 + i] = spaceDim; 
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
      
      // All indices from 2 to 2 + fieldRank, if any, are bounded by spaceDim
      for(int i = 0; i < fieldRank; i++){
        newIndexRange[2 + i] = spaceDim; 
      }
      
      // We know that for Dk operatorRank = 1 and so there's just one more upper bound left
      // given by the cardinality of the set of all derivatives of order k
      newIndexRange[2 + fieldRank] = getDkCardinality(getOperatorOrder(operatorType),spaceDim);
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
                          ">>> ERROR (VarContainer): Invalid operator type");    
  }
  
  // Once newIndexRange has been set, resize LexContainer (this also resets the index range)
  // and change field type, operator type and space dim:
  this->resize(newIndexRange);
  fieldType_    = fieldType;
  operatorType_ = operatorType;
  spaceDim_     = spaceDim;
}



template<class Scalar>
inline VarContainer<Scalar>& VarContainer<Scalar>::operator = (const VarContainer<Scalar>& right)
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this == &right ),
                      std::invalid_argument,
                      ">>> ERROR (VarContainer): Invalid right-hand side to '='. Self-assignment prohibited.");
#endif
  this->reset(right.getNumPoints(), right.getNumFields(), right.fieldType_, right.operatorType_, right.spaceDim_);
  this->data_ = right.data_;
  return *this;
}



template<class Scalar>
void VarContainer<Scalar>::print(std::ostream& os) const {
  
  // Save the format state of the original ostream os.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(os);  
  os.setf(std::ios_base::right);
  
  int size = this -> getSize();
  int rank = this -> getRank();
  Teuchos::Array<int> multiIndex(rank);
  Teuchos::Array<int> indexRange;
  this -> getIndexRange(indexRange);
  
  os<< "===============================================================================\n"\
    << "\t Container size = " << size << "   rank = " << rank << "\n"\
    << "\t Stores " << size << " values of " \
    << EOperatorToString(this -> getOperatorType()) \
    << "(" 
    << EFieldToString(this -> getFieldType()) 
    << ") for a set of " 
    << this -> getNumFields() << " fields, evaluated at " \
    << this -> getNumPoints() << " points in " << this -> getSpaceDim() << "D. \n";
  
  
  if( ((rank == 0 ) && (size == 0)) ) {
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
    if(size == 0) {
      os<< "===============================================================================\n"\
      << "|                     *** This container has no data ****                     |\n";
    }
  }
  os<< "===============================================================================\n\n";
  
  // reset format state of os
  os.copyfmt(oldFormatState);
  
  return os;
}


//===========================================================================//
//                                                                           //
//           END of member definitions; START friends and related            //
//                                                                           //
//===========================================================================//


template<class Scalar>
std::ostream& operator << (std::ostream& os, const VarContainer<Scalar>& container) {
  
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
    << "\t Container size = " << size << "   rank = " << rank << "\n"\
    << "\t Stores " << size << " values of " \
    << EOperatorToString(container.getOperatorType()) \
    << "(" 
    << EFieldToString(container.getFieldType()) 
    << ") for a set of " 
    << container.getNumFields() << " fields, evaluated at " \
    << container.getNumPoints() << " points in " << container.getSpaceDim() << "D. \n";
  
  
  if( ((rank == 0 ) && (size == 0)) ) {
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
    if(size == 0) {
      os<< "===============================================================================\n"\
      << "|                     *** This container has no data ****                     |\n";
    }
    else {      
      os<< "===============================================================================\n"\
      << "|Enumeration            Value                  Multi-index                    |\n"\
      << "===============================================================================\n";
    }
  }
  
  for(int enumeration = 0; enumeration < size; enumeration++){
    container.getMultiIndex(multiIndex,enumeration);
    os<< std::setiosflags(std::ios::left) << std::setw(5) << enumeration << "\t" \
      << std::setw(myprec+8) << "\t" << container[enumeration] << "\t\t\t\t\t" ;
    for(int r = 0; r < rank; r++){
      os <<  multiIndex[r] << " "; 
    }
    os << "\n";
  }
  
  os<< "===============================================================================\n\n";
  
  // reset format state of os
  os.copyfmt(oldFormatState);

  return os;
}

// End member, friend, and related function definitions of class VarContainer.

} // end namespace Intrepid
