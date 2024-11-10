// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CubatureTensorDef.hpp
    \brief  Definition file for the Intrepid2::CubatureTensor class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_TENSOR_DEF_HPP__
#define __INTREPID2_CUBATURE_TENSOR_DEF_HPP__

namespace Intrepid2 {

  template<typename DT, typename PT, typename WT>
  template<typename cubPointValueType,  class ...cubPointProperties,
           typename cubWeightValueType, class ...cubWeightProperties>
  void 
  CubatureTensor<DT,PT,WT>::
  getCubatureImpl( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
                   Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights ) const {
#ifdef HAVE_INTREPID2_DEBUG
    // check size of cubPoints and cubWeights
    INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(cubPoints.extent(0))  < this->getNumPoints() ||
                                  static_cast<ordinal_type>(cubPoints.extent(1))  < this->getDimension() ||
                                  static_cast<ordinal_type>(cubWeights.extent(0)) < this->getNumPoints(), std::out_of_range,
                                  ">>> ERROR (CubatureTensor): Insufficient space allocated for cubature points or weights.");
#endif
    using cubPointViewType = Kokkos::DynRankView<cubPointValueType, DT>;
    using cubWeightViewType = Kokkos::DynRankView<cubWeightValueType,DT>;

    // mirroring and where the data is problematic... when it becomes a problem, then deal with it.
    cubPointViewType  tmpPoints [Parameters::MaxTensorComponents];
    cubWeightViewType tmpWeights[Parameters::MaxTensorComponents];

    // this temporary allocation can be member of cubature; for now, let's do this way.
    // this is cubature setup on the reference cell and called for tensor elements.
    for (auto k=0;k<this->numCubatures_;++k) {
      const auto &cub = this->cubatures_[k];
      tmpPoints [k] = cubPointViewType ("CubatureTensor::getCubature::tmpPoints",  cub.getNumPoints(), cub.getDimension());
      tmpWeights[k] = cubWeightViewType("CubatureTensor::getCubature::tmpWeights", cub.getNumPoints());
      cub.getCubature(tmpPoints[k], tmpWeights[k]);
    }      
    
    // when the input containers are device space, this is better computed on host and copy to devices
    // fill tensor cubature
    {
      ordinal_type offset[Parameters::MaxTensorComponents+1] = {};
      for (auto k=0;k<this->numCubatures_;++k) {
        offset[k+1] = offset[k] + this->cubatures_[k].getDimension();
      }
      ordinal_type ii = 0, i[3] = {};

      cubPointViewType cubPoints_device("cubPoints_device", cubPoints.extent(0), cubPoints.extent(1));
      cubWeightViewType cubWeights_device("cubWeights_device", cubWeights.extent(0));

      auto cubPoints_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), cubPoints_device);
      auto cubWeights_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), cubWeights_device);

      Kokkos::deep_copy(cubPoints_host, 0.0);
      Kokkos::deep_copy(cubWeights_host, 1.0);

      switch (this->numCubatures_) {
      case 2: {
        const ordinal_type npts[] = { this->cubatures_[0].getNumPoints(), this->cubatures_[1].getNumPoints() };
        const ordinal_type dim [] = { this->cubatures_[0].getDimension(), this->cubatures_[1].getDimension() };

        /// here create mirror view and copy does not work        
        typename cubWeightViewType::HostMirror tmpWeights_host[2];
        tmpWeights_host[0] = Kokkos::create_mirror_view(tmpWeights[0]);
        tmpWeights_host[1] = Kokkos::create_mirror_view(tmpWeights[1]);
        Kokkos::deep_copy(tmpWeights_host[0], tmpWeights[0]);
        Kokkos::deep_copy(tmpWeights_host[1], tmpWeights[1]);

        typename cubPointViewType::HostMirror tmpPoints_host[2];
        tmpPoints_host[0] = Kokkos::create_mirror_view(tmpPoints[0]);
        tmpPoints_host[1] = Kokkos::create_mirror_view(tmpPoints[1]);

        Kokkos::deep_copy(tmpPoints_host[0], tmpPoints[0]);
        Kokkos::deep_copy(tmpPoints_host[1], tmpPoints[1]);

        for (i[1]=0;i[1]<npts[1];++i[1])
          for (i[0]=0;i[0]<npts[0];++i[0]) {
            for (auto nc=0;nc<2;++nc) {
              cubWeights_host(ii) *= tmpWeights_host[nc](i[nc]);
              for (ordinal_type j=0;j<dim[nc];++j)  
                cubPoints_host(ii, offset[nc]+j) = tmpPoints_host[nc](i[nc], j);
            }
            ++ii;
          }
        break;
      }
      case 3: {
        const ordinal_type npts[] = { this->cubatures_[0].getNumPoints(), this->cubatures_[1].getNumPoints(), this->cubatures_[2].getNumPoints() };
        const ordinal_type dim [] = { this->cubatures_[0].getDimension(), this->cubatures_[1].getDimension(), this->cubatures_[2].getDimension() };

        /// here create mirror view and copy does not work        
        typename cubWeightViewType::HostMirror tmpWeights_host[3];
        tmpWeights_host[0] = Kokkos::create_mirror_view(tmpWeights[0]);
        tmpWeights_host[1] = Kokkos::create_mirror_view(tmpWeights[1]);
        tmpWeights_host[2] = Kokkos::create_mirror_view(tmpWeights[2]);

        Kokkos::deep_copy(tmpWeights_host[0], tmpWeights[0]);
        Kokkos::deep_copy(tmpWeights_host[1], tmpWeights[1]);
        Kokkos::deep_copy(tmpWeights_host[2], tmpWeights[2]);

        typename cubPointViewType::HostMirror tmpPoints_host[3];
        tmpPoints_host[0] = Kokkos::create_mirror_view(tmpPoints[0]);
        tmpPoints_host[1] = Kokkos::create_mirror_view(tmpPoints[1]);
        tmpPoints_host[2] = Kokkos::create_mirror_view(tmpPoints[2]);

        Kokkos::deep_copy(tmpPoints_host[0], tmpPoints[0]);
        Kokkos::deep_copy(tmpPoints_host[1], tmpPoints[1]);
        Kokkos::deep_copy(tmpPoints_host[2], tmpPoints[2]);

        for (i[2]=0;i[2]<npts[2];++i[2])
          for (i[1]=0;i[1]<npts[1];++i[1])
            for (i[0]=0;i[0]<npts[0];++i[0]) {
              for (auto nc=0;nc<3;++nc) {             
                cubWeights_host(ii) *= tmpWeights_host[nc](i[nc]);
                for (ordinal_type j=0;j<dim[nc];++j) 
                  cubPoints_host(ii, offset[nc]+j) = tmpPoints_host[nc](i[nc], j);
              }
              ++ii;
            }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( !(this->numCubatures_ == 2 || this->numCubatures_ == 3), std::runtime_error,
                                      ">>> ERROR (CubatureTensor::getCubature): CubatureTensor supports only 2 or 3 component direct cubatures.");
      }
      }
      Kokkos::deep_copy(cubPoints_device, cubPoints_host);
      Kokkos::deep_copy(cubWeights_device, cubWeights_host);

      Kokkos::deep_copy(cubPoints, cubPoints_device);
      Kokkos::deep_copy(cubWeights, cubWeights_device);
    }
  }

} // end namespace Intrepid2

#endif
