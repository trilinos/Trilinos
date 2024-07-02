// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  Intrepid2_DataDimensionInfo.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 5/31/23.
//

#ifndef Intrepid2_DataDimensionInfo_hpp
#define Intrepid2_DataDimensionInfo_hpp

/** \file  Intrepid2_DataDimensionInfo.hpp
   \brief  Defines DimensionInfo struct that allows specification of a dimension within a Data object.
   \author Created by N.V. Roberts.
*/

#include "Intrepid2_DataVariationType.hpp"

namespace Intrepid2
{
/** \struct  Intrepid2::DimensionInfo
    \brief Struct expressing all variation information about a Data object in a single dimension, including its logical extent and storage extent.
*/
  struct DimensionInfo
  {
    int logicalExtent;
    DataVariationType variationType;
    int dataExtent;
    int variationModulus; // should be equal to dataExtent variationType other than MODULAR and CONSTANT
    int blockPlusDiagonalLastNonDiagonal = -1; // only relevant for variationType == BLOCK_PLUS_DIAGONAL
  };

  //! Returns DimensionInfo for a Data container that combines (through multiplication, say, or addition) the two specified DimensionInfo specifications in one of its dimensions.
  KOKKOS_INLINE_FUNCTION
  DimensionInfo combinedDimensionInfo(const DimensionInfo &myData, const DimensionInfo &otherData)
  {
    const int myNominalExtent    = myData.logicalExtent;
#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(myNominalExtent != otherData.logicalExtent, std::invalid_argument, "both Data objects must match in their logical extent in the specified dimension");
#endif
    const DataVariationType & myVariation    = myData.variationType;
    const DataVariationType & otherVariation = otherData.variationType;
    
    const int & myVariationModulus    = myData.variationModulus;
    const int & otherVariationModulus = otherData.variationModulus;
    
    int myDataExtent    = myData.dataExtent;
    int otherDataExtent = otherData.dataExtent;
    
    DimensionInfo combinedDimensionInfo;
    combinedDimensionInfo.logicalExtent = myNominalExtent;
    
    switch (myVariation)
    {
      case CONSTANT:
        switch (otherVariation)
        {
          case CONSTANT:
          case MODULAR:
          case GENERAL:
          case BLOCK_PLUS_DIAGONAL:
            combinedDimensionInfo = otherData;
        }
        break;
      case MODULAR:
        switch (otherVariation)
        {
          case CONSTANT:
            combinedDimensionInfo = myData;
            break;
          case MODULAR:
            if (myVariationModulus == otherVariationModulus)
            {
              // in this case, expect both to have the same data extent
              INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(myDataExtent != otherDataExtent, std::logic_error, "Unexpected data extent/modulus combination");
              combinedDimensionInfo.variationType    = MODULAR;
              combinedDimensionInfo.dataExtent       = myDataExtent;
              combinedDimensionInfo.variationModulus = myVariationModulus;
            }
            else
            {
              // both modular with two different moduli
              // we could do something clever with e.g. least common multiples, but for now we just use GENERAL
              // (this is not a use case we anticipate being a common one)
              combinedDimensionInfo.variationType    = GENERAL;
              combinedDimensionInfo.dataExtent       = myNominalExtent;
              combinedDimensionInfo.variationModulus = myNominalExtent;
            }
            break;
          case BLOCK_PLUS_DIAGONAL:
            combinedDimensionInfo.variationType    = GENERAL;
            combinedDimensionInfo.dataExtent       = myNominalExtent;
            combinedDimensionInfo.variationModulus = myNominalExtent;
            break;
          case GENERAL:
            // otherData is GENERAL: its info dominates
            combinedDimensionInfo = otherData;
            break;
        }
        break;
      case BLOCK_PLUS_DIAGONAL:
        switch (otherVariation)
        {
          case CONSTANT:
            combinedDimensionInfo = myData;
            break;
          case MODULAR:
            combinedDimensionInfo.variationType    = GENERAL;
            combinedDimensionInfo.dataExtent       = myNominalExtent;
            combinedDimensionInfo.variationModulus = myNominalExtent;
            break;
          case GENERAL:
            // otherData is GENERAL: its info dominates
            combinedDimensionInfo = otherData;
            break;
          case BLOCK_PLUS_DIAGONAL:
            combinedDimensionInfo.variationType    = GENERAL;
            combinedDimensionInfo.dataExtent       = max(myDataExtent,otherDataExtent);
            combinedDimensionInfo.variationModulus = combinedDimensionInfo.dataExtent;
            // for this case, we want to take the minimum of the two Data objects' blockPlusDiagonalLastNonDiagonal as the combined object's blockPlusDiagonalLastNonDiagonal
            combinedDimensionInfo.blockPlusDiagonalLastNonDiagonal = min(myData.blockPlusDiagonalLastNonDiagonal, otherData.blockPlusDiagonalLastNonDiagonal);
        }
        break;
      case GENERAL:
        switch (otherVariation)
        {
          case CONSTANT:
          case MODULAR:
          case GENERAL:
          case BLOCK_PLUS_DIAGONAL:
            combinedDimensionInfo = myData;
        }
    }
    return combinedDimensionInfo;
  }

} // end namespace Intrepid2

#endif /* Intrepid2_DataDimensionInfo_hpp */
