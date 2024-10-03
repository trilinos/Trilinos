// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  Intrepid2_DataVariationType.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 5/31/23.
//

#ifndef Intrepid2_DataVariationType_hpp
#define Intrepid2_DataVariationType_hpp

/** \file  Intrepid2_DataVariationType.hpp
   \brief  Defines DataVariationType enum that specifies the types of variation possible within a Data object.
   \author Created by N.V. Roberts.
*/

namespace Intrepid2 {
  template<class DataScalar,typename DeviceType>
  class Data;

/** \enum  Intrepid2::DataVariationType
    \brief Enumeration to indicate how data varies in a particular dimension of an Intrepid2::Data object.
 CONSTANT indicates that the data does not vary; MODULAR indicates that it varies according to some (separately specified) modulus; BLOCK_PLUS_DIAGONAL allows specification of a matrix that has a non-diagonal block followed by a diagonal block; GENERAL indicates arbitrary variation.
 
 To give some examples for a Data object containing the Jacobian for reference-to-physical space mappings:
 - CONSTANT could be used in the point dimension for an affine transformation
 - MODULAR could be used for the cell dimension for a uniform grid that has been subdivided into simplices
 - BLOCK_PLUS_DIAGONAL could be used for the coordinate dimensions for an arbitrary 2D mesh that has been orthogonally extruded in the z dimension (resulting in diagonal entries in the final row and column of the Jacobian matrix)
 - GENERAL should be used in any dimension in which the data varies in a way not captured by the other options
*/
  enum DataVariationType
  {
    CONSTANT            /// does not vary
  , MODULAR             /// varies according to modulus of the index
  , BLOCK_PLUS_DIAGONAL /// one of two dimensions in a matrix; bottom-right part of matrix is diagonal
  , GENERAL             /// arbitrary variation
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DataVariationType_hpp */
