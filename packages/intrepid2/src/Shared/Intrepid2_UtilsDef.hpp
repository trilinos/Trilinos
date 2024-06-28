// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_UtilsDef.hpp
    \brief  Definition file for Intrepid2 utilities.
    \author Created by P. Bochev and D. Ridzal and Kyungjoo Kim.
*/

//
// functions here are moved to basis class
//

#ifndef __INTREPID2_UTILS_DEF_HPP__
#define __INTREPID2_UTILS_DEF_HPP__

#include "Intrepid2_Utils.hpp"

namespace Intrepid2 {

  //--------------------------------------------------------------------------------------------//
  //                                                                                            //
  //    Definitions: functions for orders, cardinality and etc. of differential operators       //
  //                                                                                            //
  //--------------------------------------------------------------------------------------------//

  KOKKOS_INLINE_FUNCTION
  ordinal_type getFieldRank(const EFunctionSpace spaceType) {
    ordinal_type fieldRank = -1;

    switch (spaceType) {

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
      INTREPID2_TEST_FOR_ABORT( !isValidFunctionSpace(spaceType),
                                ">>> ERROR (Intrepid2::getFieldRank): Invalid function space type");
    }
    return fieldRank;
  }


  KOKKOS_INLINE_FUNCTION
  ordinal_type getOperatorRank(const EFunctionSpace spaceType,
                               const EOperator      operatorType,
                               const ordinal_type   spaceDim) {

    const auto fieldRank = Intrepid2::getFieldRank(spaceType);
#ifdef HAVE_INTREPID2_DEBUG
    // Verify arguments: field rank can be 0,1, or 2, spaceDim can be 1,2, or 3.
    INTREPID2_TEST_FOR_ABORT( !(0 <= fieldRank && fieldRank <= 2),
                              ">>> ERROR (Intrepid2::getOperatorRank): Invalid field rank");
    INTREPID2_TEST_FOR_ABORT( !(1 <= spaceDim && spaceDim  <= 3),
                              ">>> ERROR (Intrepid2::getOperatorRank): Invalid space dimension");
#endif
    ordinal_type operatorRank = -999;
    
    // In 1D GRAD, CURL, and DIV default to d/dx; Dk defaults to d^k/dx^k, no casing needed.
    if (spaceDim == 1) {
      if (fieldRank == 0) {
        // By default, in 1D any operator other than VALUE has rank 1
        if (operatorType == OPERATOR_VALUE) {
          operatorRank = 0;
        } else {
          operatorRank = 1;
        }
      }
      
      // Only scalar fields are allowed in 1D
      else {
        INTREPID2_TEST_FOR_ABORT( fieldRank > 0,
                                  ">>> ERROR (getOperatorRank): Only scalar fields are allowed in 1D");
      } // fieldRank == 0
    } // spaceDim == 1
    
    // We are either in 2D or 3D
    else {
      switch (operatorType) {
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
        if (fieldRank > 0) {
          operatorRank = spaceDim - 3;
        } else {
          
          // CURL can be applied to scalar functions (rank = 0) in 2D and gives a vector (rank = 1)
          if (spaceDim == 2) {
            operatorRank = 3 - spaceDim;
          }
          
          // If we are here, fieldRank=0, spaceDim=3: CURL is undefined for 3D scalar functions
          else {
            INTREPID2_TEST_FOR_ABORT( ( (spaceDim == 3) && (fieldRank == 0) ),
                                      ">>> ERROR (Intrepid2::getOperatorRank): CURL cannot be applied to scalar fields in 3D");
          }
        }
        break;
        
      case OPERATOR_DIV:

        // DIV can be applied to vectors and tensors and has rank -1 in 2D and 3D
        if (fieldRank > 0) {
          operatorRank = -1;
        }

        // DIV cannot be applied to scalar fields except in 1D where it defaults to d/dx
        else {
          INTREPID2_TEST_FOR_ABORT( ( (spaceDim > 1) && (fieldRank == 0) ),
                                    ">>> ERROR (Intrepid2::getOperatorRank): DIV cannot be applied to scalar fields in 2D and 3D");
        }
        break;

      default:
        INTREPID2_TEST_FOR_ABORT( !( isValidOperator(operatorType) ),
                                  ">>> ERROR (Intrepid2::getOperatorRank): Invalid operator type");
      } // switch
    }// 2D and 3D

    return operatorRank;
  }


  KOKKOS_INLINE_FUNCTION
  ordinal_type getOperatorOrder(const EOperator operatorType) {
    ordinal_type opOrder = -1;

    switch (operatorType) {

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
      opOrder = (ordinal_type)operatorType - (ordinal_type)OPERATOR_D1 + 1;
      break;

    default:
      INTREPID2_TEST_FOR_ABORT( !( Intrepid2::isValidOperator(operatorType) ),
                                ">>> ERROR (Intrepid2::getOperatorOrder): Invalid operator type");
    }
    return opOrder;
  }


  KOKKOS_INLINE_FUNCTION
  ordinal_type getDkEnumeration(const ordinal_type xMult,
                                const ordinal_type yMult,
                                const ordinal_type zMult) {

    if (yMult < 0 && zMult < 0) {

#ifdef HAVE_INTREPID2_DEBUG
      // We are in 1D: verify input - xMult is non-negative  and total order <= 10:
      INTREPID2_TEST_FOR_ABORT( !( (0 <= xMult) && (xMult <= Parameters::MaxDerivative) ),
                                ">>> ERROR (Intrepid2::getDkEnumeration): Derivative order out of range");
#endif

      // there's only one derivative of order xMult
      return 0;
    } else {
      if (zMult < 0) {

#ifdef HAVE_INTREPID2_DEBUG
        // We are in 2D: verify input - xMult and yMult are non-negative and total order <= 10:
        INTREPID2_TEST_FOR_ABORT( !(0 <= xMult && 0 <= yMult && (xMult + yMult) <= Parameters::MaxDerivative),
                                  ">>> ERROR (Intrepid2::getDkEnumeration): Derivative order out of range");
#endif

        // enumeration is the value of yMult
        return yMult;
      }

      // We are in 3D: verify input - xMult, yMult and zMult are non-negative and total order <= 10:
      else {
        const auto order = xMult + yMult + zMult;
#ifdef HAVE_INTREPID2_DEBUG
        // Verify input:  total order cannot exceed 10:
        INTREPID2_TEST_FOR_ABORT(  !( (0 <= xMult) && (0 <= yMult) && (0 <= zMult) &&
                                      (order <= Parameters::MaxDerivative) ),
                                   ">>> ERROR (Intrepid2::getDkEnumeration): Derivative order out of range");
#endif
        ordinal_type enumeration = zMult;
        const ordinal_type iend = order-xMult+1;
        for(ordinal_type i=0;i<iend;++i) {
          enumeration += i;
        }
        return enumeration;
      }
    }
  }

//   template<typename OrdinalArrayType>
//   KOKKOS_INLINE_FUNCTION
//   void getDkMultiplicities(OrdinalArrayType   partialMult,
//                            const ordinal_type derivativeEnum,
//                            const EOperator    operatorType,
//                            const ordinal_type spaceDim) {

//     /* Hash table to convert enumeration of partial derivative to multiplicities of dx,dy,dz in 3D.
//        Multiplicities {mx,my,mz} are arranged lexicographically in bins numbered from 0 to 10.
//        The size of bins is an arithmetic progression, i.e., 1,2,3,4,5,...,11. Conversion formula is:
//        \verbatim
//        mx = derivativeOrder - binNumber
//        mz = derivativeEnum  - binBegin
//        my = derivativeOrder - mx - mz = binNumber + binBegin - derivativeEnum
//        \endverbatim
//        where binBegin is the enumeration of the first element in the bin. Bin numbers and binBegin
//        values are stored in hash tables for quick access by derivative enumeration value.
//     */

//     // Returns the bin number for the specified derivative enumeration
//     const ordinal_type binNumber[66] = {
//       0,
//       1, 1,
//       2, 2, 2,
//       3, 3, 3, 3,
//       4, 4, 4, 4, 4,
//       5, 5, 5, 5, 5, 5,
//       6, 6, 6, 6, 6, 6, 6,
//       7, 7, 7, 7, 7, 7, 7, 7,
//       8, 8, 8, 8, 8, 8, 8, 8, 8,
//       9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
//       10,10,10,10,10,10,10,10,10,10,10
//     };

//     // Returns the binBegin value for the specified derivative enumeration
//     const ordinal_type binBegin[66] = {
//       0,
//       1, 1,
//       3, 3 ,3,
//       6, 6, 6, 6,
//       10,10,10,10,10,
//       15,15,15,15,15,15,
//       21,21,21,21,21,21,21,
//       28,28,28,28,28,28,28,28,
//       36,36,36,36,36,36,36,36,36,
//       45,45,45,45,45,45,45,45,45,45,
//       55,55,55,55,55,55,55,55,55,55,55
//     };

// #ifdef HAVE_INTREPID2_DEBUG
//     // Enumeration value must be between 0 and the cardinality of the derivative set
//     INTREPID2_TEST_FOR_ABORT( !( (0 <= derivativeEnum) && (derivativeEnum < getDkCardinality(operatorType,spaceDim) ) ),
//                               ">>> ERROR (Intrepid2::getDkMultiplicities): Invalid derivative enumeration value for this order and space dimension");
// #endif

//     // This method should only be called for Dk operators
//     ordinal_type derivativeOrder;
//     switch(operatorType) {

//     case OPERATOR_D1:
//     case OPERATOR_D2:
//     case OPERATOR_D3:
//     case OPERATOR_D4:
//     case OPERATOR_D5:
//     case OPERATOR_D6:
//     case OPERATOR_D7:
//     case OPERATOR_D8:
//     case OPERATOR_D9:
//     case OPERATOR_D10:
//       derivativeOrder = Intrepid2::getOperatorOrder(operatorType);
//       break;

//     default:
//       INTREPID2_TEST_FOR_ABORT(true,
//                                ">>> ERROR (Intrepid2::getDkMultiplicities): operator type Dk required for this method");
//     }// switch

// #ifdef HAVE_INTREPID2_DEBUG
//     INTREPID2_TEST_FOR_ABORT( partialMult.extent(0) != spaceDim,
//                               ">>> ERROR (Intrepid2::getDkMultiplicities): partialMult must have the same dimension as spaceDim" );
// #endif

//     switch(spaceDim) {

//     case 1:
//       // Multiplicity of dx equals derivativeOrder
//       partialMult(0) = derivativeOrder;
//       break;

//     case 2:
//       // Multiplicity of dy equals the enumeration of the derivative; of dx - the complement
//       partialMult(1) = derivativeEnum;
//       partialMult(0) = derivativeOrder - derivativeEnum;
//       break;

//     case 3:
//       // Recover multiplicities
//       partialMult(0) = derivativeOrder - binNumber[derivativeEnum];
//       partialMult(1) = binNumber[derivativeEnum] + binBegin[derivativeEnum] - derivativeEnum;
//       partialMult(2) = derivativeEnum  -  binBegin[derivativeEnum];
//       break;

//     default:
//       INTREPID2_TEST_FOR_ABORT( !( (0 < spaceDim ) && (spaceDim < 4) ),
//                                 ">>> ERROR (Intrepid2::getDkMultiplicities): Invalid space dimension");
//     }
//   }


  KOKKOS_INLINE_FUNCTION
  ordinal_type getDkCardinality(const EOperator    operatorType,
                                const ordinal_type spaceDim) {

    // This should only be called for Dk operators
    ordinal_type derivativeOrder;
    switch(operatorType) {

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
      derivativeOrder = Intrepid2::getOperatorOrder(operatorType);
      break;

    default:
      INTREPID2_TEST_FOR_ABORT(true,
                               ">>> ERROR (Intrepid2::getDkCardinality): operator type Dk required for this method");
    }// switch

    ordinal_type cardinality = -999;
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
      INTREPID2_TEST_FOR_ABORT( !( (0 < spaceDim ) && (spaceDim < 4) ),
                                ">>> ERROR (Intrepid2::getDkcardinality): Invalid space dimension");
    }

    return cardinality;
  }

} // end namespace Intrepid2

#endif
