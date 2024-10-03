// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_BasisDef.hpp
    \brief  Implementation file for the abstract base class Intrepid2::Basis.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_BASIS_DEF_HPP__
#define __INTREPID2_BASIS_DEF_HPP__

#include <stdexcept>

namespace Intrepid2 {

  //--------------------------------------------------------------------------------------------//
  //                                                                                            //
  //                            Helper functions of the Basis class                             //
  //                                                                                            //
  //--------------------------------------------------------------------------------------------//

  //
  //    functions for orders, cardinality and etc. of differential operators
  //

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

  template<EOperator operatorType>
  KOKKOS_INLINE_FUNCTION
  constexpr ordinal_type getOperatorOrder() {
    return (operatorType == OPERATOR_VALUE) ? 0 :
           ((operatorType == OPERATOR_GRAD) || (operatorType == OPERATOR_CURL) || (operatorType == OPERATOR_DIV) || (operatorType == OPERATOR_D1)) ? 1 :
           (ordinal_type)operatorType - (ordinal_type)OPERATOR_D1 + 1;
  }


  template<ordinal_type spaceDim>
  KOKKOS_INLINE_FUNCTION
  ordinal_type getDkEnumeration(const ordinal_type /*xMult*/,
                                const ordinal_type yMult,
                                const ordinal_type zMult) {
    switch(spaceDim) {

      case 1: return 0;
      case 2: return yMult;
      case 3: return zMult + (yMult+zMult)*(yMult+zMult+1)/2;

      default: {
        INTREPID2_TEST_FOR_ABORT( !( (0 < spaceDim ) && (spaceDim < 4) ),
              ">>> ERROR (Intrepid2::getDkEnumeration): Invalid space dimension");
        return 0;
      }
    }
  }

  template<ordinal_type spaceDim>
  KOKKOS_INLINE_FUNCTION
  ordinal_type getPnEnumeration(const ordinal_type p,
                                const ordinal_type q /* = 0*/,
                                const ordinal_type r /* = 0*/) {
    return (spaceDim==1) ? p :
           (spaceDim==2) ? (p+q)*(p+q+1)/2+q :
                           (p+q+r)*(p+q+r+1)*(p+q+r+2)/6+(q+r)*(q+r+1)/2+r;

  }

  template<typename value_type>
  KOKKOS_INLINE_FUNCTION
  void getJacobyRecurrenceCoeffs (
            value_type  &an,
            value_type  &bn,
            value_type  &cn,
      const ordinal_type alpha,
      const ordinal_type beta ,
      const ordinal_type n) {
    an = ( (2.0 * n + 1.0 + alpha + beta) * ( 2.0 * n + 2.0 + alpha + beta ) /
        value_type(2.0 * ( n + 1 ) * ( n + 1 + alpha + beta ) ) );
    bn = ( (alpha*alpha-beta*beta)*(2.0*n+1.0+alpha+beta) /
        value_type(2.0*(n+1.0)*(2.0*n+alpha+beta)*(n+1.0+alpha+beta) ) );
    cn = ( (n+alpha)*(n+beta)*(2.0*n+2.0+alpha+beta) /
        value_type( (n+1.0)*(n+1.0+alpha+beta)*(2.0*n+alpha+beta) ) );
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

#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_ABORT( !( (0 < spaceDim ) && (spaceDim < 8) ),
                                    ">>> ERROR (Intrepid2::getDkcardinality): Invalid space dimension");
    switch (operatorType) {
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
          break;
        default:
          INTREPID2_TEST_FOR_ABORT( true, ">>> ERROR (Intrepid2::getDkCardinality): Cannot be used for this operator ");
          break;
        }
#endif
    
    ordinal_type n = Intrepid2::getOperatorOrder(operatorType);
    return (spaceDim==1) ? 1 :
           (spaceDim==2) ?  n + 1 :
           (spaceDim==3) ? (n + 1) * (n + 2) / 2 :
           (spaceDim==4) ? (n + 1) * (n + 2) * (n + 3) / 6 :
           (spaceDim==5) ? (n + 1) * (n + 2) * (n + 3) * (n + 4) / 24 :
           (spaceDim==6) ? (n + 1) * (n + 2) * (n + 3) * (n + 4) * (n + 5) / 120 :
                           (n + 1) * (n + 2) * (n + 3) * (n + 4) * (n + 5) * (n + 6) / 720;
  }

  template<EOperator operatorType, ordinal_type spaceDim>
  KOKKOS_INLINE_FUNCTION
  constexpr ordinal_type getDkCardinality() {
    return getPnCardinality<spaceDim-1,Intrepid2::getOperatorOrder<operatorType>()>();
  }

  template<ordinal_type spaceDim>
  KOKKOS_INLINE_FUNCTION
  ordinal_type getPnCardinality (ordinal_type n) {

#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_ABORT( !( (0 <= spaceDim ) && (spaceDim < 4) ),
                                        ">>> ERROR (Intrepid2::getPnCardinality): Invalid space dimension");
#endif

    return (spaceDim==0) ? 1 :
           (spaceDim==1) ? n+1 :
           (spaceDim==2) ? (n + 1) * (n + 2) / 2 :
           (n + 1) * (n + 2) * (n + 3) / 6;
  }


  template<ordinal_type spaceDim, ordinal_type n>
  KOKKOS_INLINE_FUNCTION
  constexpr ordinal_type getPnCardinality () {

    return (spaceDim==0) ? 1 :
           (spaceDim==1) ? n+1 :
           (spaceDim==2) ? (n + 1) * (n + 2) / 2 :
           (n + 1) * (n + 2) * (n + 3) / 6;

  }




  //
  // Argument check
  //


  template<typename outputValueViewType,
           typename inputPointViewType>
  void getValues_HGRAD_Args( const outputValueViewType   outputValues,
                             const inputPointViewType    inputPoints,
                             const EOperator             operatorType,
                             const shards::CellTopology  cellTopo,
                             const ordinal_type          basisCard ) {
    const auto spaceDim = cellTopo.getDimension();

    // Verify inputPoints array
    INTREPID2_TEST_FOR_EXCEPTION( !(rank(inputPoints) == 2), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) rank = 2 required for inputPoints array");

    INTREPID2_TEST_FOR_EXCEPTION(  (inputPoints.extent(0) <= 0), std::invalid_argument,
                                   ">>> ERROR (Intrepid2::getValues_HGRAD_Args): dim 0 (number of points) > 0 required for inputPoints array");

    INTREPID2_TEST_FOR_EXCEPTION( !(inputPoints.extent(1) == spaceDim), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) dim 1 (spatial dimension) of inputPoints array  does not match cell dimension");


    // Verify that all inputPoints are in the reference cell
    /*
      INTREPID2_TEST_FOR_EXCEPTION( !CellTools<Scalar>::checkPointSetInclusion(inputPoints, cellTopo), std::invalid_argument,
      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) One or more points are outside the "
      << cellTopo <<" reference cell");
    */


    // Verify that operatorType is admissible for HGRAD fields
    INTREPID2_TEST_FOR_EXCEPTION( ( (spaceDim == 2) && (operatorType == OPERATOR_DIV) ), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) DIV is invalid operator for rank-0 (scalar) fields in 2D.");

    INTREPID2_TEST_FOR_EXCEPTION( ( (spaceDim == 3) && ( (operatorType == OPERATOR_DIV) ||
                                                         (operatorType == OPERATOR_CURL) ) ), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) DIV and CURL are invalid operators for rank-0 (scalar) fields in 3D.");


    // Check rank of outputValues (all operators are admissible in 1D) and its dim 2 when operator is
    // GRAD, CURL (only in 2D), or Dk.
    
    if(spaceDim == 1) {
      switch(operatorType){
      case OPERATOR_VALUE:
        INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 2), std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) rank = 2 required for outputValues when operator = VALUE.");
        break;
      case OPERATOR_GRAD:
      case OPERATOR_CURL:
      case OPERATOR_DIV:
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
        INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 3), std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) rank = 3 required for outputValues in 1D when operator = GRAD, CURL, DIV, or Dk.");

        INTREPID2_TEST_FOR_EXCEPTION( !(outputValues.extent(2) == 1 ),
                                      std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) dim 2 of outputValues must equal 1 when operator = GRAD, CURL, DIV, or Dk.");
        break;
      default:
        INTREPID2_TEST_FOR_EXCEPTION( (true), std::invalid_argument, ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) Invalid operator");
      }
    }
    else if(spaceDim > 1) {
      switch(operatorType){
      case OPERATOR_VALUE:
        INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 2), std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) rank = 2 required for outputValues when operator = VALUE.");
        break;
      case OPERATOR_GRAD:
      case OPERATOR_CURL:
      case OPERATOR_D1:
        INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 3), std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) rank = 3 required for outputValues in 2D and 3D when operator = GRAD, CURL (in 2D), or Dk.");

        INTREPID2_TEST_FOR_EXCEPTION( !(outputValues.extent(2) == spaceDim ),
                                      std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) dim 2 of outputValues must equal cell dimension when operator = GRAD, CURL (in 2D), or D1.");
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
        INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 3), std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) rank = 3 required for outputValues in 2D and 3D when operator = GRAD, CURL (in 2D), or Dk.");
        
        INTREPID2_TEST_FOR_EXCEPTION( !(static_cast<ordinal_type>(outputValues.extent(2)) == getDkCardinality(operatorType, spaceDim)),
                                      std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) dim 2 of outputValues must equal cardinality of the Dk multiset.");
        break;
      default:
        INTREPID2_TEST_FOR_EXCEPTION( (true), std::invalid_argument, ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) Invalid operator");
      }
    }
    
    
    // Verify dim 0 and dim 1 of outputValues
    INTREPID2_TEST_FOR_EXCEPTION( !(outputValues.extent(1) == inputPoints.extent(0) ),
                                  std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) dim 1 (number of points) of outputValues must equal dim 0 of inputPoints.");
    
    INTREPID2_TEST_FOR_EXCEPTION( !(static_cast<ordinal_type>(outputValues.extent(0)) == basisCard ),
                                  std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) dim 0 (number of basis functions) of outputValues must equal basis cardinality.");
  }
  
  
  template<typename outputValueViewType,
           typename inputPointViewType>
  void getValues_HCURL_Args( const outputValueViewType   outputValues,
                             const inputPointViewType    inputPoints,
                             const EOperator             operatorType,
                             const shards::CellTopology  cellTopo,
                             const ordinal_type          basisCard ) {
    
    const auto spaceDim = cellTopo.getDimension();
    
    // Verify that cell is 2D or 3D (this is redundant for default bases where we use correct cells,
    //  but will force any user-defined basis for HCURL spaces to use 2D or 3D cells
    INTREPID2_TEST_FOR_EXCEPTION( !( (spaceDim == 2) || (spaceDim == 3) ), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HCURL_Args) cell dimension = 2 or 3 required for HCURL basis");
    
    
    // Verify inputPoints array
    INTREPID2_TEST_FOR_EXCEPTION( !(rank(inputPoints) == 2), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HCURL_Args) rank = 2 required for inputPoints array");
    INTREPID2_TEST_FOR_EXCEPTION(  (inputPoints.extent(0) <= 0), std::invalid_argument,
                                   ">>> ERROR (Intrepid2::getValues_HCURL_Args): dim 0 (number of points) > 0 required for inputPoints array");
    
    INTREPID2_TEST_FOR_EXCEPTION( !(inputPoints.extent(1) == spaceDim), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HCURL_Args) dim 1 (spatial dimension) of inputPoints array  does not match cell dimension");
    
    // Verify that all inputPoints are in the reference cell
    /*
      INTREPID2_TEST_FOR_EXCEPTION( !CellTools<Scalar>::checkPointSetInclusion(inputPoints, cellTopo), std::invalid_argument,
      ">>> ERROR: (Intrepid2::getValues_HCURL_Args) One or more points are outside the "
      << cellTopo <<" reference cell");
    */
    
    // Verify that operatorType is admissible for HCURL fields
    INTREPID2_TEST_FOR_EXCEPTION( !( (operatorType == OPERATOR_VALUE) || (operatorType == OPERATOR_CURL) ), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HCURL_Args) operator = VALUE or CURL required for HCURL fields.");
    
    
    // Check rank of outputValues: for VALUE should always be rank-3 array with (F,P,D) layout
    switch(operatorType) {
      
    case OPERATOR_VALUE:
      INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 3), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::getValues_HCURL_Args) rank = 3 required for outputValues when operator is VALUE");
      INTREPID2_TEST_FOR_EXCEPTION( !(outputValues.extent(2) == spaceDim ),
                                    std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::getValues_HCURL_Args) dim 2 of outputValues must equal cell dimension when operator is VALUE.");
      break;
      
    case OPERATOR_CURL:
      
      // in 3D we need an (F,P,D) container because CURL gives a vector field:
      if(spaceDim == 3) {
        INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 3 ) ,
                                      std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HCURL_Args) rank = 3 required for outputValues in 3D when operator is CURL");
        INTREPID2_TEST_FOR_EXCEPTION( !(outputValues.extent(2) == spaceDim),
                                      std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HCURL_Args) dim 2 of outputValues must equal cell dimension in 3D when operator is CURL.");
      }
      // In 2D we need an (F,P) container because CURL gives a scalar field
      else if(spaceDim == 2) {
        INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 2 ) ,
                                      std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HCURL_Args) rank = 2 required for outputValues in 2D when operator is CURL");
      }
      break;
      
    default:
      INTREPID2_TEST_FOR_EXCEPTION( (true), std::invalid_argument, ">>> ERROR: (Intrepid2::getValues_HCURL_Args) Invalid operator");
    }
    
    
    // Verify dim 0 and dim 1 of outputValues
    INTREPID2_TEST_FOR_EXCEPTION( !(outputValues.extent(1) == inputPoints.extent(0) ),
                                  std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HCURL_Args) dim 1 (number of points) of outputValues must equal dim 0 of inputPoints.");
    
    INTREPID2_TEST_FOR_EXCEPTION( !(static_cast<ordinal_type>(outputValues.extent(0)) == basisCard ),
                                  std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HCURL_Args) dim 0 (number of basis functions) of outputValues must equal basis cardinality.");
    
  }
  
  
  
  template<typename outputValueViewType,
           typename inputPointViewType>
  void getValues_HDIV_Args( const outputValueViewType   outputValues,
                            const inputPointViewType    inputPoints,
                            const EOperator             operatorType,
                            const shards::CellTopology  cellTopo,
                            const ordinal_type          basisCard ) {
    
    const auto spaceDim = cellTopo.getDimension();
    
    // Verify inputPoints array
    INTREPID2_TEST_FOR_EXCEPTION( !(rank(inputPoints) == 2), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HDIV_Args) rank = 2 required for inputPoints array");
    INTREPID2_TEST_FOR_EXCEPTION(  (inputPoints.extent(0) <= 0), std::invalid_argument,
                                   ">>> ERROR (Intrepid2::getValues_HDIV_Args): dim 0 (number of points) > 0 required for inputPoints array");
    
    INTREPID2_TEST_FOR_EXCEPTION( !(inputPoints.extent(1) == spaceDim), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HDIV_Args) dim 1 (spatial dimension) of inputPoints array  does not match cell dimension");
    
    // Verify that all inputPoints are in the reference cell
    /*
      INTREPID2_TEST_FOR_EXCEPTION( !CellTools<Scalar>::checkPointSetInclusion(inputPoints, cellTopo), std::invalid_argument,
      ">>> ERROR: (Intrepid2::getValues_HDIV_Args) One or more points are outside the "
      << cellTopo <<" reference cell");
    */
    
    // Verify that operatorType is admissible for HDIV fields
    INTREPID2_TEST_FOR_EXCEPTION( !( (operatorType == OPERATOR_VALUE) || (operatorType == OPERATOR_DIV) ), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HDIV_Args) operator = VALUE or DIV required for HDIV fields.");
    
    
    // Check rank of outputValues
    switch(operatorType) {
    case OPERATOR_VALUE:
      INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 3), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::getValues_HDIV_Args) rank = 3 required for outputValues when operator is VALUE.");
      
      INTREPID2_TEST_FOR_EXCEPTION( !(outputValues.extent(2) == spaceDim ),
                                    std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::getValues_HDIV_Args) dim 2 of outputValues must equal cell dimension for operator VALUE.");
      break;
    case OPERATOR_DIV:
      INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 2), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::getValues_HDIV_Args) rank = 2 required for outputValues when operator is DIV.");
      break;
      
    default:
      INTREPID2_TEST_FOR_EXCEPTION( (true), std::invalid_argument, ">>> ERROR: (Intrepid2::getValues_HDIV_Args) Invalid operator");
    }
    
    
    // Verify dim 0 and dim 1 of outputValues
    INTREPID2_TEST_FOR_EXCEPTION( !(outputValues.extent(1) == inputPoints.extent(0) ),
                                  std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HDIV_Args) dim 1 (number of points) of outputValues must equal dim 0 of inputPoints.");
    
    INTREPID2_TEST_FOR_EXCEPTION( !(outputValues.extent(0) == static_cast<size_type>(basisCard) ),
                                  std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HDIV_Args) dim 0 (number of basis functions) of outputValues must equal basis cardinality.");
  }
  
  template<typename outputValueViewType,
           typename inputPointViewType>
  void getValues_HVOL_Args( const outputValueViewType   outputValues,
                             const inputPointViewType    inputPoints,
                             const EOperator             operatorType,
                             const shards::CellTopology  cellTopo,
                             const ordinal_type          basisCard ) {
    const auto spaceDim = cellTopo.getDimension();

    // Verify inputPoints array
    INTREPID2_TEST_FOR_EXCEPTION( !(rank(inputPoints) == 2), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) rank = 2 required for inputPoints array");

    INTREPID2_TEST_FOR_EXCEPTION(  (inputPoints.extent(0) <= 0), std::invalid_argument,
                                   ">>> ERROR (Intrepid2::getValues_HGRAD_Args): dim 0 (number of points) > 0 required for inputPoints array");

    INTREPID2_TEST_FOR_EXCEPTION( !(inputPoints.extent(1) == spaceDim), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) dim 1 (spatial dimension) of inputPoints array  does not match cell dimension");


    // Verify that all inputPoints are in the reference cell
    /*
      INTREPID2_TEST_FOR_EXCEPTION( !CellTools<Scalar>::checkPointSetInclusion(inputPoints, cellTopo), std::invalid_argument,
      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) One or more points are outside the "
      << cellTopo <<" reference cell");
    */


    // Verify that operatorType is admissible for HGRAD fields
    INTREPID2_TEST_FOR_EXCEPTION( ( (spaceDim == 2) && (operatorType == OPERATOR_DIV) ), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) DIV is invalid operator for rank-0 (scalar) fields in 2D.");

    INTREPID2_TEST_FOR_EXCEPTION( ( (spaceDim == 3) && ( (operatorType == OPERATOR_DIV) ||
                                                         (operatorType == OPERATOR_CURL) ) ), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) DIV and CURL are invalid operators for rank-0 (scalar) fields in 3D.");


    // Check rank of outputValues (all operators are admissible in 1D) and its dim 2 when operator is
    // GRAD, CURL (only in 2D), or Dk.

    if(spaceDim == 1) {
      switch(operatorType){
      case OPERATOR_VALUE:
        INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 2), std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) rank = 2 required for outputValues when operator = VALUE.");
        break;
      case OPERATOR_GRAD:
      case OPERATOR_CURL:
      case OPERATOR_DIV:
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
        INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 3), std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) rank = 3 required for outputValues in 1D when operator = GRAD, CURL, DIV, or Dk.");

        INTREPID2_TEST_FOR_EXCEPTION( !(outputValues.extent(2) == 1 ),
                                      std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) dim 2 of outputValues must equal 1 when operator = GRAD, CURL, DIV, or Dk.");
        break;
      default:
        INTREPID2_TEST_FOR_EXCEPTION( (true), std::invalid_argument, ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) Invalid operator");
      }
    }
    else if(spaceDim > 1) {
      switch(operatorType){
      case OPERATOR_VALUE:
        INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 2), std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) rank = 2 required for outputValues when operator = VALUE.");
        break;
      case OPERATOR_GRAD:
      case OPERATOR_CURL:
      case OPERATOR_D1:
        INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 3), std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) rank = 3 required for outputValues in 2D and 3D when operator = GRAD, CURL (in 2D), or Dk.");

        INTREPID2_TEST_FOR_EXCEPTION( !(outputValues.extent(2) == spaceDim ),
                                      std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) dim 2 of outputValues must equal cell dimension when operator = GRAD, CURL (in 2D), or D1.");
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
        INTREPID2_TEST_FOR_EXCEPTION( !(rank(outputValues) == 3), std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) rank = 3 required for outputValues in 2D and 3D when operator = GRAD, CURL (in 2D), or Dk.");

        INTREPID2_TEST_FOR_EXCEPTION( !(static_cast<ordinal_type>(outputValues.extent(2)) == getDkCardinality(operatorType, spaceDim)),
                                      std::invalid_argument,
                                      ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) dim 2 of outputValues must equal cardinality of the Dk multiset.");
        break;
      default:
        INTREPID2_TEST_FOR_EXCEPTION( (true), std::invalid_argument, ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) Invalid operator");
      }
    }


    // Verify dim 0 and dim 1 of outputValues
    INTREPID2_TEST_FOR_EXCEPTION( !(outputValues.extent(1) == inputPoints.extent(0) ),
                                  std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) dim 1 (number of points) of outputValues must equal dim 0 of inputPoints.");

    INTREPID2_TEST_FOR_EXCEPTION( !(static_cast<ordinal_type>(outputValues.extent(0)) == basisCard ),
                                  std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::getValues_HGRAD_Args) dim 0 (number of basis functions) of outputValues must equal basis cardinality.");
  }

  template<typename Device,
            typename outputValueType,
            typename pointValueType>
  Kokkos::DynRankView<outputValueType,Device>
  Basis<Device,outputValueType,pointValueType>::allocateOutputView( const int numPoints, const EOperator operatorType) const
  {
    const bool operatorIsDk = (operatorType >= OPERATOR_D1) && (operatorType <= OPERATOR_D10);
    const bool operatorSupported = (operatorType == OPERATOR_VALUE) || (operatorType == OPERATOR_GRAD) || (operatorType == OPERATOR_CURL) || (operatorType == OPERATOR_DIV) || operatorIsDk;
    INTREPID2_TEST_FOR_EXCEPTION(!operatorSupported, std::invalid_argument, "operator is not supported by allocateOutputView()");
    
    const int numFields = this->getCardinality();
    const int spaceDim  = this->getBaseCellTopology().getDimension() + this->getNumTensorialExtrusions();
    
    using OutputViewAllocatable = Kokkos::DynRankView<outputValueType,DeviceType>;
    
    switch (functionSpace_)
    {
      case FUNCTION_SPACE_HGRAD:
        if (operatorType == OPERATOR_VALUE)
        {
          // scalar-valued container
          OutputViewAllocatable dataView("BasisValues HGRAD VALUE data", numFields, numPoints);
          return dataView;
        }
        else if (operatorType == OPERATOR_GRAD)
        {
          OutputViewAllocatable dataView("BasisValues HGRAD GRAD data", numFields, numPoints, spaceDim);
          return dataView;
        }
        else if (operatorIsDk)
        {
          ordinal_type dkCardinality = getDkCardinality(operatorType, spaceDim);
          OutputViewAllocatable dataView("BasisValues HGRAD Dk data", numFields, numPoints, dkCardinality);
          return dataView;
        }
        else
        {
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "operator/space combination not supported by allocateOutputView()");
        }
      case FUNCTION_SPACE_HDIV:
        if (operatorType == OPERATOR_VALUE)
        {
          // vector-valued container
          OutputViewAllocatable dataView("BasisValues HDIV VALUE data", numFields, numPoints, spaceDim);
          return dataView;
        }
        else if (operatorType == OPERATOR_DIV)
        {
          // scalar-valued curl
          OutputViewAllocatable dataView("BasisValues HDIV DIV data", numFields, numPoints);
          return dataView;
        }
        else
        {
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "operator/space combination not supported by allocateOutputView()");
        }
      case FUNCTION_SPACE_HCURL:
        if (operatorType == OPERATOR_VALUE)
        {
          OutputViewAllocatable dataView("BasisValues HCURL VALUE data", numFields, numPoints, spaceDim);
          return dataView;
        }
        else if (operatorType == OPERATOR_CURL)
        {
          if (spaceDim != 2)
          {
            // vector-valued curl
            OutputViewAllocatable dataView("BasisValues HCURL CURL data", numFields, numPoints, spaceDim);
            return dataView;
          }
          else
          {
            // scalar-valued curl
            OutputViewAllocatable dataView("BasisValues HCURL CURL data (scalar)", numFields, numPoints);
            return dataView;
          }
        }
        else
        {
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "operator/space combination not supported by allocateOutputView()");
        }
      case FUNCTION_SPACE_HVOL:
        if (operatorType == OPERATOR_VALUE)
        {
          // vector-valued container
          OutputViewAllocatable dataView("BasisValues HVOL VALUE data", numFields, numPoints);
          return dataView;
        }
        else if (operatorIsDk || (operatorType == OPERATOR_GRAD))
        {
          ordinal_type dkCardinality = getDkCardinality(operatorType, spaceDim);
          OutputViewAllocatable dataView("BasisValues HVOL Dk data", numFields, numPoints, dkCardinality);
          return dataView;
        }
        else
        {
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "operator/space combination not supported by allocateOutputView()");
        }
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "operator/space combination not supported by allocateOutputView()");
    }
  }
}

#endif
