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

/** \file   Intrepid_CubatureTensorDef.hpp
    \brief  Definition file for the Intrepid2::CubatureTensor class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_TENSOR_DEF_HPP__
#define __INTREPID2_CUBATURE_TENSOR_DEF_HPP__

namespace Intrepid2 {

  template<typename SpT, typename PT, typename WT>
  template<typename CT0, typename CT1>
  CubatureTensor<SpT,PT,WT>::
  CubatureTensor( const CT0 ct0,
                  const CT1 ct1 ) 
    : numCubatures_(2),
      dimension_(ct0.getDimension()+ct1.getDimension()),
      impl_(this) {
    cubatures_[0] = ct0;
    cubatures_[1] = ct1;
  }


  template<typename SpT, typename PT, typename WT>
  template<typename CT0, typename CT1, typename CT2>
  CubatureTensor<SpT,PT,WT>::
  CubatureTensor( const CT0 ct0,
                  const CT1 ct1,
                  const CT2 ct2 ) 
    : numCubatures_(3),
      dimension_(ct0.getDimension()+ct1.getDimension()+ct2.getDimension()),
      impl_(this) { 
    cubatures_[0] = ct0;
    cubatures_[1] = ct1;
    cubatures_[2] = ct2;
  }


  template<typename SpT, typename PT, typename WT>
  template<typename cubPointValueType,  class ...cubPointProperties,
           typename cubWeightValueType, class ...cubWeightProperties>
  void 
  CubatureTensor<SpT,PT,WT>::Internal::
  getCubature( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
               Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights ) const {
#ifdef HAVE_INTREPID2_DEBUG
    // check size of cubPoints and cubWeights
    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.dimension(0)  < obj_->getNumPoints() ||
                                  cubPoints.dimension(1)  < obj_->getDimension() ||
                                  cubWeights.dimension(0) < obj_->getNumPoints(), std::out_of_range,
                                  ">>> ERROR (CubatureTensor): Insufficient space allocated for cubature points or weights.");
#endif
    typedef Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPointViewType;
    typedef Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeightViewType; 

    // mirroring and where the data is problematic... when it becomes a problem, then deal with it.
    cubPointViewType  tmpPoints [Parameters::MaxDimension];
    cubWeightViewType tmpWeights[Parameters::MaxDimension];

    // this temporary allocation can be member of cubature; for now, let's do this way.
    // this is cubature setup on the reference cell and called for tensor elements.
    for (auto k=0;k<obj_->numCubatures_;++k) {
      const auto &cub = obj_->cubatures_[k];
      tmpPoints [k] = cubPointViewType ("CubatureTensor::getCubature::tmpPoints",  cub.getNumPoints(), cub.getDimension());
      tmpWeights[k] = cubWeightViewType("CubatureTensor::getCubature::tmpWeights", cub.getNumPoints());
      cub.getCubature(tmpPoints[k], tmpWeights[k]);
    }      
    
    {
      const auto npts = obj_->getNumPoints();
      const auto dim = obj_->getDimension();
      
      const Kokkos::pair<ordinal_type,ordinal_type> pointRange(0, npts);
      const Kokkos::pair<ordinal_type,ordinal_type> dimRange(0, dim);
      
      Kokkos::deep_copy(Kokkos::subdynrankview(cubPoints, pointRange, dimRange), 0.0);
      Kokkos::deep_copy(Kokkos::subdynrankview(cubWeights, pointRange), 1.0);
    }

    // when the input containers are device space, this is better computed on host and copy to devices
    // fill tensor cubature
    {
      ordinal_type offset[Parameters::MaxDimension+1] = {};
      for (auto k=0;k<obj_->numCubatures_;++k) {
        offset[k+1] = offset[k] + obj_->cubatures_[k].getDimension();
      }
      ordinal_type ii = 0, i[3] = {};
      switch (obj_->numCubatures_) {
      case 2: {
        const ordinal_type npts[] = { obj_->cubatures_[0].getNumPoints(), obj_->cubatures_[1].getNumPoints() };
        const ordinal_type dim [] = { obj_->cubatures_[0].getDimension(), obj_->cubatures_[1].getDimension() };

        for (i[1]=0;i[1]<npts[1];++i[1])
          for (i[0]=0;i[0]<npts[0];++i[0]) {
            for (auto nc=0;nc<2;++nc) {
              cubWeights(ii) *= tmpWeights[nc](i[nc]);
              for (auto j=0;j<dim[nc];++j)  
                cubPoints(ii, offset[nc]+j) = tmpPoints[nc](i[nc], j);
            }
            ++ii;
          }
        break;
      }
      case 3: {
        const ordinal_type npts[] = { obj_->cubatures_[0].getNumPoints(), obj_->cubatures_[1].getNumPoints(), obj_->cubatures_[2].getNumPoints() };
        const ordinal_type dim [] = { obj_->cubatures_[0].getDimension(), obj_->cubatures_[1].getDimension(), obj_->cubatures_[2].getDimension() };

        for (i[2]=0;i[2]<npts[2];++i[2])
          for (i[1]=0;i[1]<npts[1];++i[1])
            for (i[0]=0;i[0]<npts[0];++i[0]) {
              for (auto nc=0;nc<3;++nc) {             
                cubWeights(ii) *= tmpWeights[nc](i[nc]);
                for (auto j=0;j<dim[nc];++j) 
                  cubPoints(ii, offset[nc]+j) = tmpPoints[nc](i[nc], j);
              }
              ++ii;
            }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( obj_->numCubatures_ != 2 || obj_->numCubatures_ != 3, std::runtime_error,
                                      ">>> ERROR (CubatureTensor::getCubature): CubatureTensor supports only 2 or 3 component direct cubatures.");
      }
      }
    }
  }

} // end namespace Intrepid2

#endif
