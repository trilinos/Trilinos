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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubatureTensor.hpp
    \brief  Header file for the Intrepid2::CubatureTensor class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_TENSOR_HPP__
#define __INTREPID2_CUBATURE_TENSOR_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_CubatureDirect.hpp"

namespace Intrepid2 {

  /** \class Intrepid2::CubatureTensor
      \brief Defines tensor-product cubature (integration) rules in Intrepid.
  */
  template<typename ExecSpaceType = void>
  class CubatureTensor : public Cubature<ExecSpaceType> {
  private:

    /** \brief Array of cubature rules, stored as FieldContainers.
     */
    ordinal_type numCubatures_;

    CubatureDirect<ExecSpaceType> cubatures_[Parameters::MaxDimension];
    
    /** \brief Dimension of integration domain.
     */
    ordinal_type dimension_;
  
  public:

    CubatureTensor() = delete;
    ~CubatureTensor() = default;


    /** \brief Constructor.
        
        \param cubature1        [in]     - First direct cubature rule.
        \param cubature2        [in]     - Second direct cubature rule.
    */
    template<typename CubatureType0, 
             typename CubatureType1>
    CubatureTensor( const CubatureType0 cubature0,
                    const CubatureType1 cubature1 );
    
    /** \brief Constructor.
        
        \param cubature1        [in]     - First direct cubature rule.
        \param cubature2        [in]     - Second direct cubature rule.
        \param cubature3        [in]     - Third direct cubature rule.
    */
    template<typename CubatureType0, 
             typename CubatureType1,
             typename CubatureType2>
    CubatureTensor( const CubatureType0 cubature0,
                    const CubatureType1 cubature1,
                    const CubatureType2 cubature2 );
    
    
    /** \brief Returns cubature points and weights
        (return arrays must be pre-sized/pre-allocated).
        
        \param cubPoints       [out]     - Vector containing the cubature points.
        \param cubWeights      [out]     - Vector of corresponding cubature weights.
    */

    template<typename cubPointValueType,  class ...cubPointProperties,
             typename cubWeightValueType, class ...cubWeightProperties>
    void
    getCubature( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
                 Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights ) const;

    
    /** \brief Returns cubature points and weights.
        Method for physical space cubature, throws an exception.
        
        \param cubPoints             [out]        - Array containing the cubature points.
        \param cubWeights            [out]        - Array of corresponding cubature weights.
        \param cellCoords             [in]        - Array of cell coordinates
    */
    template<typename cubPointValueType,  class ...cubPointProperties,
             typename cubWeightValueType, class ...cubWeightProperties,
             typename cellCoordValueType, class ...cellCoordProperties>
    void
    getCubature( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
                 Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights,
                 Kokkos::DynRankView<cellCoordValueType,cellCoordProperties...> cellCoords ) const;
    
    /** \brief Returns the number of cubature points.
     */
    ordinal_type getNumPoints() const;
    
    /** \brief Returns dimension of integration domain.
     */
    ordinal_type getDimension() const;

    /** \brief Returns cubature name.
     */
    const char* getName() const;

    /** \brief Return the number of cubatures.
     */
    ordinal_type getNumCubatures() const;
    
    /** \brief Returns max. degree of polynomials that are integrated exactly.
     */
    void getAccuracy( ordinal_type *accuracy ) const;
  };
  
  
} // end namespace Intrepid2


// include templated definitions
#include <Intrepid2_CubatureTensorDef.hpp>

#endif
