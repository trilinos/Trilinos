// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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
// Questions? Contact Mauro Perego  (mperego@sandia.gov) or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef Intrepid2_DataTools_h
#define Intrepid2_DataTools_h

#include "Intrepid2_Data.hpp"

/** \file  Intrepid2_DataTools.hpp
   \brief  Utility methods for manipulating Intrepid2::Data objects.
   \author Created by N.V. Roberts.
*/

namespace Intrepid2 {

class DataTools
{
public:
  //! Fills Data object of logical shape (C,P,D,D) corresponding to the pointwise product of an object of shape (C,P,D,D) with one of shape (C,P).
  //! Will also work for any "matrix" data and scalar data where the "matrix" has equal or higher rank than the scalar, but otherwise share the same shape.  E.g., (C,F,P,D1,D2) matrices could be multiplied by (C,F,P) scalars, or (C,P,D) vectors by (C,P) scalars.
  //! \param resultMatrixData [out] -  the resulting (C,P,D,D) container.  Must be allocated appropriately to store the resulting data; see the implementation of the two-argument multiplyByCPWeights(), which performs this allocation for you.
  //! \param matrixDataIn [in] - the input (C,P,D,D) container.
  //! \param scalarDataIn [in] - the input (C,P) container.
  template<class Scalar, class DeviceType>
  static void multiplyByCPWeights(Data<Scalar,DeviceType> &resultMatrixData, const Data<Scalar,DeviceType> &matrixDataIn, const Data<Scalar,DeviceType> &scalarDataIn)
  {
    const ordinal_type rank       = scalarDataIn.rank();
    const ordinal_type matrixRank = matrixDataIn.rank();
    auto extents                  = scalarDataIn.getExtents();
    auto variationTypes           = scalarDataIn.getVariationTypes();
    for (int r=rank; r<matrixRank; r++)
    {
      extents[r]        = matrixDataIn.extent_int(r);
      variationTypes[r] = CONSTANT;
    }
    
    auto scalarDataInExtended = scalarDataIn.shallowCopy(matrixRank, extents, variationTypes);
    resultMatrixData.storeInPlaceProduct(matrixDataIn,scalarDataInExtended);
  }
  
  //! Allocates and fills Data object of logical shape (C,P,D,D) corresponding to the pointwise product of an object of shape (C,P,D,D) with one of shape (C,P).
  //! Will also work for any "matrix" data and scalar data where the "matrix" has equal or higher rank than the scalar, but otherwise share the same shape.  E.g., (C,F,P,D1,D2) matrices could be multiplied by (C,F,P) scalars, or (C,P,D) vectors by (C,P) scalars.
  //! \param matrixDataIn [in] - the (C,P,D,D) container.
  //! \param scalarDataIn [in] - the (C,P) container.
  //! \return
  template<class Scalar, class DeviceType>
  static Data<Scalar,DeviceType> multiplyByCPWeights(const Data<Scalar,DeviceType> &matrixDataIn, const Data<Scalar,DeviceType> &scalarDataIn)
  {
    const ordinal_type rank       = scalarDataIn.rank();
    const ordinal_type matrixRank = matrixDataIn.rank();
    auto extents                  = scalarDataIn.getExtents();
    auto variationTypes           = scalarDataIn.getVariationTypes();
    for (int r=rank; r<matrixRank; r++)
    {
      extents[r]        = matrixDataIn.extent_int(r);
      variationTypes[r] = CONSTANT;
    }
    
    auto scalarDataInExtended = scalarDataIn.shallowCopy(matrixRank, extents, variationTypes);
    
    auto result = Data<Scalar,DeviceType>::allocateInPlaceCombinationResult(scalarDataInExtended, matrixDataIn);
    
    result.storeInPlaceProduct(matrixDataIn,scalarDataInExtended);
    return result;
  }
  
  //! Allocates and fills Data object corresponding to the transpose of matrix data, represented by the last two dimensions of the input object.
  //! \param matrixDataIn [in] - the (…,D1,D2) container.
  //! \return a (…,D2,D1) container containing the transpose of the input matrix data.
  template<class Scalar, class DeviceType>
  static Data<Scalar,DeviceType> transposeMatrix(const Data<Scalar,DeviceType> &matrixDataIn)
  {
    // A direct construction of the transpose could be more efficient, but here we take advantage of existing
    // implementations within the Data class supporting matrix-matrix multiplication.  We construct an identity
    // matrix, and left-multiply this by the transpose of the input matrix.
    const ordinal_type rank = matrixDataIn.rank();
    auto extents            = matrixDataIn.getExtents();
    auto variationTypes     = matrixDataIn.getVariationTypes();
    const auto D1           = extents[rank-2];
    
    extents[rank-2]        = D1;
    extents[rank-1]        = D1;
    variationTypes[rank-2] = BLOCK_PLUS_DIAGONAL;
    variationTypes[rank-1] = BLOCK_PLUS_DIAGONAL;
    
    Kokkos::View<Scalar*,DeviceType> identityUnderlyingView("Intrepid2::DataTools::transposeMatrix() - identity view",D1);
    Kokkos::deep_copy(identityUnderlyingView, 1.0);
    Data<Scalar,DeviceType> identityData(identityUnderlyingView,extents,variationTypes);
    
    auto result = Data<Scalar,DeviceType>::allocateMatMatResult(true, matrixDataIn, false, identityData);
    result.storeMatMat(true, matrixDataIn, false, identityData);
    
    return result;
  }
};
}

#endif /* Intrepid2_DataTools_h */
