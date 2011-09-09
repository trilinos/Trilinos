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

/** \file   Intrepid_Utils.cpp
    \brief  Implementation file for Intrepid_Utils.hpp.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_UTILS_CPP
#define INTREPID_UTILS_CPP

#include "Intrepid_Utils.hpp"

namespace Intrepid {
  
//--------------------------------------------------------------------------------------------//
//                                                                                            //
//    Definitions: functions for orders, cardinality and etc. of differential operators       //
//                                                                                            //
//--------------------------------------------------------------------------------------------//

int getFieldRank(const EFunctionSpace spaceType) {
  int fieldRank = -1;
  
  switch(spaceType){
    
    case FUNCTION_SPACE_HGRAD:
    case FUNCTION_SPACE_HVOL:
      fieldRank = 0;
      break;
      
    case FUNCTION_SPACE_HCURL:
    case FUNCTION_SPACE_HDIV:
    case FUNCTION_SPACE_VECTOR_HGRAD:
      fieldRank = 1;
      break;
      
    case FUNCTION_SPACE_TENSOR_HGRAD:
      fieldRank = 2;
      break;
      
    default:
      TEST_FOR_EXCEPTION( !( isValidFunctionSpace(spaceType) ), std::invalid_argument,
                          ">>> ERROR (Intrepid::getFieldRank): Invalid function space type");
  }
  return fieldRank;
}



int getOperatorRank(const EFunctionSpace spaceType,
                    const EOperator      operatorType,
                    const int            spaceDim) {
  
  int fieldRank = Intrepid::getFieldRank(spaceType);
  
  // Verify arguments: field rank can be 0,1, or 2, spaceDim can be 1,2, or 3.
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION( !( (0 <= fieldRank) && (fieldRank <= 2) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid::getOperatorRank): Invalid field rank");
  TEST_FOR_EXCEPTION( !( (1 <= spaceDim ) && (spaceDim  <= 3) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid::getOperatorRank): Invalid space dimension");
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
            TEST_FOR_EXCEPTION( ( (spaceDim == 3) && (fieldRank == 0) ), std::invalid_argument,
                                ">>> ERROR (Intrepid::getOperatorRank): CURL cannot be applied to scalar fields in 3D");  
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
          TEST_FOR_EXCEPTION( ( (spaceDim > 1) && (fieldRank == 0) ), std::invalid_argument,
                              ">>> ERROR (Intrepid::getOperatorRank): DIV cannot be applied to scalar fields in 2D and 3D");  
        }
        break;
        
      default:
        TEST_FOR_EXCEPTION( !( isValidOperator(operatorType) ), std::invalid_argument,
                            ">>> ERROR (Intrepid::getOperatorRank): Invalid operator type");
    } // switch
  }// 2D and 3D
  
  return operatorRank; 
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
      TEST_FOR_EXCEPTION( !( Intrepid::isValidOperator(operatorType) ),
                          std::invalid_argument,
                          ">>> ERROR (Intrepid::getOperatorOrder): Invalid operator type");
  }
  return opOrder;
}    



int getDkEnumeration(const int xMult,
                     const int yMult,
                     const int zMult) {
  
  if( (yMult < 0) && (zMult < 0)) {
    
#ifdef HAVE_INTREPID_DEBUG
    // We are in 1D: verify input - xMult is non-negative  and total order <= 10:
    TEST_FOR_EXCEPTION( !( (0 <= xMult) && (xMult <= INTREPID_MAX_DERIVATIVE) ), std::out_of_range,
                        ">>> ERROR (Intrepid::getDkEnumeration): Derivative order out of range");
#endif
    
    // there's only one derivative of order xMult
    return 0;
  }
  else {
    if( zMult < 0 ) {
      
#ifdef HAVE_INTREPID_DEBUG
      // We are in 2D: verify input - xMult and yMult are non-negative and total order <= 10:
      TEST_FOR_EXCEPTION( !( (0 <= xMult) && (0 <= yMult) && 
                             ( (xMult + yMult) <= INTREPID_MAX_DERIVATIVE) ), std::out_of_range,
                          ">>> ERROR (Intrepid::getDkEnumeration): Derivative order out of range");
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
                              (order <= INTREPID_MAX_DERIVATIVE) ), std::out_of_range,
                           ">>> ERROR (Intrepid::getDkEnumeration): Derivative order out of range");
#endif
      int enumeration = zMult;
      for(int i = 0; i < order - xMult + 1; i++){
        enumeration += i; 
      }
      return enumeration;
    }
  }
}



void getDkMultiplicities(Teuchos::Array<int>&  partialMult,
                         const int             derivativeEnum,
                         const EOperator       operatorType,
                         const int             spaceDim) {
  
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
    
#ifdef HAVE_INTREPID_DEBUG
  // Enumeration value must be between 0 and the cardinality of the derivative set
  TEST_FOR_EXCEPTION( !( (0 <= derivativeEnum) && (derivativeEnum < getDkCardinality(operatorType,spaceDim) ) ),
                      std::invalid_argument,
                      ">>> ERROR (Intrepid::getDkMultiplicities): Invalid derivative enumeration value for this order and space dimension");
#endif
  
  // This method should only be called for Dk operators
  int derivativeOrder;
  switch(operatorType){
    
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
      derivativeOrder = Intrepid::getOperatorOrder(operatorType);
      break;
      
    default:
      TEST_FOR_EXCEPTION(true, std::invalid_argument,
                         ">>> ERROR (Intrepid::getDkMultiplicities): operator type Dk required for this method");
  }// switch
  
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
      TEST_FOR_EXCEPTION( !( (0 < spaceDim ) && (spaceDim < 4) ), std::invalid_argument,
                          ">>> ERROR (Intrepid::getDkMultiplicities): Invalid space dimension");          
  }
}



int getDkCardinality(const EOperator operatorType,
                     const int       spaceDim) {
  
  // This should only be called for Dk operators
  int derivativeOrder;
  switch(operatorType){
    
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
      derivativeOrder = Intrepid::getOperatorOrder(operatorType);
      break;
      
    default:
      TEST_FOR_EXCEPTION(true, std::invalid_argument,
                        ">>> ERROR (Intrepid::getDkCardinality): operator type Dk required for this method");
  }// switch

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
      TEST_FOR_EXCEPTION( !( (0 < spaceDim ) && (spaceDim < 4) ), std::invalid_argument,
                          ">>> ERROR (Intrepid::getDkcardinality): Invalid space dimension");          
  }
  return cardinality;
}



//--------------------------------------------------------------------------------------------//
//                                                                                            //
//            Definitions:    Helper functions of the Basis class                             //
//                                                                                            //
//--------------------------------------------------------------------------------------------//

void setOrdinalTagData(std::vector<std::vector<std::vector<int> > >    &tagToOrdinal,
                       std::vector<std::vector<int> >                  &ordinalToTag,
                       const int                                       *tags,
                       const int                                       basisCard,
                       const int                                       tagSize,
                       const int                                       posScDim,
                       const int                                       posScOrd,
                       const int                                       posDfOrd) {


  // Resize ordinalToTag to a rank-2 array with dimensions (basisCardinality_, 4) and copy tag data
  ordinalToTag.resize(basisCard);
  for (int i = 0; i < basisCard; i++) {
    ordinalToTag[i].resize(4);
    for (int j = 0; j < tagSize; j++) {
      ordinalToTag[i][j] = tags[i*tagSize + j];
    }
  }

  // Resize tagToOrdinal to a rank-3 array with dimensions (maxScDim + 1, maxScOrd + 1 , maxDfOrd +1)
  // The 1st dimension of tagToOrdinal is the max value of the 1st column (max subcell dim) in the tag + 1
  int maxScDim = 0;
  for (int i = 0; i < basisCard; i++) {
    if (maxScDim < tags[i*tagSize + posScDim]) {
      maxScDim = tags[i*tagSize + posScDim];
    }
  }
  maxScDim += 1;

  // The 2nd dimension of tagToOrdinal is the max value of the 2nd column (max subcell id) in the tag  + 1
  int maxScOrd = 0;
  for (int i = 0; i < basisCard; i++) {
    if (maxScOrd < tags[i*tagSize + posScOrd]) {
      maxScOrd = tags[i*tagSize + posScOrd];
    }
  }
  maxScOrd += 1;

  // The 3rd dimension of tagToOrdinal is the max value of the 3rd column (max subcell DofId in the tag) + 1
  int maxDfOrd = 0;
  for (int i = 0; i < basisCard; i++) {
    if (maxDfOrd < tags[i*tagSize + posDfOrd]) {
      maxDfOrd = tags[i*tagSize + posDfOrd];
    }
  }
  maxDfOrd += 1;

  // Create rank-1 array with dimension maxDfOrd (the 3rd dimension of tagToOrdinal) filled with -1
  std::vector<int> rank1Array(maxDfOrd, -1);

  // Create rank-2 array with dimensions (maxScOrd, maxDfOrd) (2nd and 3rd dimensions of tagToOrdinal)
  std::vector<std::vector<int> > rank2Array(maxScOrd, rank1Array);

  // Resize tagToOrdinal to a rank-3 array with dimensions (maxScDim, maxScOrd, maxDfOrd)
  tagToOrdinal.assign(maxScDim, rank2Array);

  // Overwrite elements of the array corresponding to tags with local DoF Id's, leave all other = -1
  for (int i = 0; i < basisCard; i++) {
    tagToOrdinal[tags[i*tagSize]][tags[i*tagSize+1]][tags[i*tagSize+2]] = i;
  }
}



} // end namespace Intrepid

#endif
