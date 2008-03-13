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
  dataContainer_ = right.dataContainer_;
  fieldType_ = right.fieldType_;
  operatorType_ = right.operatorType_;
}



template<class Scalar>
VarContainer<Scalar>::VarContainer(const int         numPoints,
                                   const int         numFields,
                                   const EField      fieldType,
                                   const EOperator   operatorType,
                                   const int         spaceDim) {  
  
  // Use reset method to adjust the internal storage and set private data
  this -> reset(numPoints,numFields,fieldType,operatorType,spaceDim);
}



template<class Scalar>
inline int VarContainer<Scalar>::getSize() const {
  return dataContainer_.getSize();
}



template<class Scalar> 
inline void VarContainer<Scalar>::getIndexRange(Teuchos::Array<int>& indexRange) const {  
  dataContainer_.getIndexRange(indexRange);
}



template<class Scalar>
inline int VarContainer<Scalar>::getIndexBound(const int whichIndex) const {
  return dataContainer_.getIndexBound(whichIndex);
}


template<class Scalar>
inline int VarContainer<Scalar>::getNumPoints() const {
  
  // Number of points is upper bound for 1st index: stored in dataContainer_.indexRange_[0]
    return dataContainer_.getIndexBound(0);
}



template<class Scalar>
inline int VarContainer<Scalar>::getNumFields() const {
  
  // Number of fields is upper bound for 2nd index: stored in dataContainer_.indexRange_[1]
  return dataContainer_.getIndexBound(1);
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
inline int VarContainer<Scalar>::getRank() const {
  return dataContainer_.getRank(); 
}



template<class Scalar>
inline int VarContainer<Scalar>::getEnumeration(const Teuchos::Array<int>& multiIndex) const {  
  return dataContainer_.getEnumeration(multiIndex);
}



template<class Scalar>
inline void VarContainer<Scalar>::getMultiIndex(Teuchos::Array<int>& multiIndex,
                                                const int            valueEnum) const {
  dataContainer_.getMultiIndex(multiIndex,valueEnum);
}



template<class Scalar>
Scalar VarContainer<Scalar>::getValue(const Teuchos::Array<int>& multiIndex) const {  
  return dataContainer_.getValue(multiIndex);
}



template<class Scalar>
inline void VarContainer<Scalar>::storeZero() {
  dataContainer_.storeZero();
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
  TEST_FOR_EXCEPTION( !( (0 <  spaceDim ) && ( spaceDim <= 3  ) ),
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
  // and change field type, operator tyepe and space dim:
  dataContainer_.resize(newIndexRange);
  fieldType_    = fieldType;
  operatorType_ = operatorType;
  spaceDim_     = spaceDim;
  
}



template<class Scalar>
inline void VarContainer<Scalar>::setValue(const Scalar               dataValue,
                                           const Teuchos::Array<int>& multiIndex) {
  dataContainer_.setValue(dataValue,multiIndex);
}



template<class Scalar>
void VarContainer<Scalar>::setValues(const Teuchos::Array<Scalar>& dataArray) {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( (dataArray.size() != (dataContainer_.size()) ),
                      std::invalid_argument,
                      ">>> ERROR (VarContainer): Size of argument does not match size of container.");  
#endif  
  dataContainer_.assign(dataArray.begin(),dataArray.end());
}



template<class Scalar>
const Scalar& VarContainer<Scalar>::operator [] (const int valueEnum) const {
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( (valueEnum < 0) || (valueEnum >= (int)dataContainer_.getSize() ) ),
                      std::out_of_range,
                      ">>> ERROR (VarContainer): value enumeration out of range.");
#endif
  return dataContainer_[valueEnum];
}
  


template<class Scalar>
inline VarContainer<Scalar>& VarContainer<Scalar>::operator = (const VarContainer<Scalar>& right)
{
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( ( this == &right ),
                      std::invalid_argument,
                      ">>> ERROR (VarContainer): Invalid right-hand side to '='. Self-assignment prohibited.");
  TEST_FOR_EXCEPTION( ( this -> getSize() != right.getSize() ),
                      std::invalid_argument,
                      ">>> ERROR (VarContainer): Invalid size of right-hand side argument to '='.");
#endif
  dataContainer_ = right.dataContainer_;
  fieldType_     = right.fieldType_;
  operatorType_  = right.OperatorType_;
  spaceDim_      = right.spaceDim_;
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
    << "\t Stores " << size << " values of " 
    << OperatorNames[this -> getOperatorType()] \
    << "(" 
    << FieldNames[this -> getFieldType()] 
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
    << "\t Stores " << size << " values of " 
    << OperatorNames[container.getOperatorType()] \
    << "(" 
    << FieldNames[container.getFieldType()] 
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


// End member, friend, and related function definitions of class VarContainer.

} // end namespace Intrepid
