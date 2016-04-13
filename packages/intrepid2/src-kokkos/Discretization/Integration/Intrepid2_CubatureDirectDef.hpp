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

/** \file   Intrepid_CubatureDirectDef.hpp
    \brief  Definition file for the Intrepid2::CubatureDirect class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

namespace Intrepid2 {

  template<typename ExecSpaceType>
  template<typename cubPointValueType,  class ...cubPointProperties,
           typename cubWeightValueType, class ...cubWeightProperties>
  void 
  CubatureDirect<ExecSpaceType>::
  getCubatureFromData( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
                       Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights,
                       const CubatureData cubData ) const {
#ifdef HAVE_INTREPID_DEBUG
    // check size of cubPoints and cubWeights
    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.dimension(0)  < this->getNumPoints() ||
                                  cubPoints.dimension(1)  < this->getDimension() ||
                                  cubWeights.dimension(0) < this->getNumPoints(), std::out_of_range,
                                  ">>> ERROR (CubatureDirect): Insufficient space allocated for cubature points or weights.");
#endif
    // need subview here
    Kokkos::deep_copy(cubPoints,  cubData.points_);
    Kokkos::deep_copy(cubWeights, cubData.weights_);
  }
  

  template<typename ExecSpaceType>
  template<typename ExecSpaceType>
  template<typename cubPointValueType,  class ...cubPointProperties,
           typename cubWeightValueType, class ...cubWeightProperties>
  void 
  CubatureDirect<ExecSpaceType>::
  getCubature( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
               Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights ) const {
    getCubatureData(cubPoints, cubWeights, getCubatureData(degree_));
  } 
  

  template<typename ExecSpaceType>
  template<typename cubPointValueType,  class ...cubPointProperties,
           typename cubWeightValueType, class ...cubWeightProperties,
           typename cellCoordValueType, class ...cellCoordProperties>
  void 
  CubatureDirect<ExecSpaceType>::
  getCubature( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
               Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights,
               Kokkos::DynRankView<cellCoordValueType,cellCoordProperties...> cellCoords ) const {
    INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                  ">>> ERROR (CubatureDirect::getCubature): Cubature defined in reference space calling method for physical space cubature.");
  }
  

  template<typename ExecSpaceType>
  ordinal_type 
  CubatureDirect<ExecSpaceType>::
  getNumPoints() const {
    return getCubatureData(degree_).numPoints_;
  } 


  template<typename ExecSpaceType>
  ordinal_type 
  CubatureDirect<ExecSpaceType>::
  getDimension() const {
    return dimension_;
  } 


  template<typename ExecSpaceType>
  ordinal_type
  CubatureDirect<ExecSpaceType>::
  getAccuracy() const {
    return degree_;
  }


  template<typename ExecSpaceType>
  ordinal_type
  CubatureDirect<ExecSpaceType>::
  getMaxAccuracy() const {
    INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                  ">>> ERROR (CubatureDirect::getMaxAccuracy): this method should be over-riden by derived classes.");
    return 0;
  }


  template<typename ExecSpaceType>
  const char*
  CubatureDirect<ExecSpaceType>::
  getName() const {
    return "CubatureDirect";
  }


    
}
