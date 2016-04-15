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

namespace Intrepid2 {

  template<typename ExecSpaceType>
  CubatureTensor<ExecSpaceType>::
  CubatureTensor( const Teuchos::RCP<CubatureDirect<ExecSpaceType> > cubature1,
                  const Teuchos::RCP<CubatureDirect<ExecSpaceType> > cubature2 ) {
    numCubatures_ = 2;
    cubatures_[0] = cubature1;
    cubatures_[1] = cubature2;

    dimension_ = 0;
    for (auto i=0;i<numCubatures;++i) {
      degree_[i]  = cubatures_[i]->getAccuracy();
      dimension_ += cubatures_[i]->getDimension();
    }
  }


  template<typename ExecSpaceType>
  CubatureTensor<ExecSpaceType>::
  CubatureTensor( const Teuchos::RCP<CubatureDirect<ExecSpaceType> > cubature1,
                  const Teuchos::RCP<CubatureDirect<ExecSpaceType> > cubature2,
                  const Teuchos::RCP<CubatureDirect<ExecSpaceType> > cubature3 ) {
    numCubatures_ = 3;
    cubatures_[0] = cubature1;
    cubatures_[1] = cubature2;
    cubatures_[2] = cubature3;

    dimension_ = 0;
    for (auto i=0;i<numCubatures;++i) {
      degree_[i]  = cubatures_[i]->getAccuracy();
      dimension_ += cubatures_[i]->getDimension();
    }
  }


  template <typename ExecSpaceType>
  template<typename cubPointValueType,  class ...cubPointProperties,
           typename cubWeightValueType, class ...cubWeightProperties>
  void 
  CubatureTensor<ExecSpaceType>::
  getCubature( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
               Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights ) const {
#ifdef HAVE_INTREPID2_DEBUG
    // check size of cubPoints and cubWeights
    INTREPID2_TEST_FOR_EXCEPTION( cubPoints.dimension(0)  < this->getNumPoints() ||
                                  cubPoints.dimension(1)  < this->getDimension() ||
                                  cubWeights.dimension(0) < this->getNumPoints(), std::out_of_range,
                                ">>> ERROR (CubatureTensor): Insufficient space allocated for cubature points or weights.");
#endif
    typedef Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  PointViewType;
    typedef Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> WeightViewType; 

    PointViewType  tmpPoints [Parameters::MaxDimension];
    WeightViewType tmpWeights[Parameters::MaxDimension];

    // cannot avoid temporary allocation here
    // this is cubature setup on the reference cell and called for tensor elements.
    for (auto k=0;k<numCubatures_;++k) {
      const auto cub = cubatures_[k];
      tmpPoints [k] = PointViewType ("CubatureTensor::getCubature::tmpPoints",  cub->getNumPoints, cub->getDimension());
      tmpWeights[k] = WeightViewType("CubatureTensor::getCubature::tmpWeights", cub->getNumPoints);
      cub->getCubature(tmpPoints[k], tmpWeights[k]);
    }      

    // reset all weights to 1.0
    //const Kokkos::pair<ordinal_type,ordinal_type> pointRange(0, this->getNumPoints());
    //Kokkos::deep_copy(Kokkos::subdynrankview(cubWeights, pointRange), 1.0);

    {
      ordinal_type ii = 0;
      for (auto k=0;k<numCubatures_;++k) {
        const auto cub  = cubatures_[k];
        const auto npts = cub->getNumPoints();
        const auto dim  = cub->getDimension();
        for (auto i=0;i<npts;++i) 
          cubWeights(ii++) = 1.0;
      }
    }

    // when the input containers are device space, this is better computed on host and copy to devices
    // fill tensor cubature
    {
      ordinal_type offset[Paramters::MaxDimension+1];
      for (auto k=0;k<numCubatures_;++k) 
        offset[k+1] = offs[k] + cubatures_[k]->getDimension();
      
      ordinal_type ii = 0;
      for (auto k=0;k<numCubatures_;++k) {
        const auto cub  = cubatures_[k];
        const auto npts = cub->getNumPoints();
        const auto dim  = cub->getDimension();
        
        const auto points  = tmpPoints[k];
        const auto weights = tmpWeights[k];
        
        const auto offs = offset[k];
        for (auto i=0;i<npts;++i) {
          for (auto j=0;j<dim;++j) 
            cubPoints(ii, offs+j) = points(i, j);
          cubWeights(ii++) *= weighst(i);
        }
      }
    }
  } 

  
  template<typename ExecSpaceType>
  template<typename cubPointValueType,  class ...cubPointProperties,
           typename cubWeightValueType, class ...cubWeightProperties,
           typename cellCoordValueType, class ...cellCoordProperties>
  void
  CubatureTensor<ExecSpaceType>::
  getCubature( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
               Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights,
               Kokkos::DynRankView<cellCoordValueType,cellCoordProperties...> cellCoords ) const {
    INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                  ">>> ERROR (CubatureTensor::getCubature): Cubature defined in reference space calling method for physical space cubature.");
  }


  template<typename ExecSpaceType>
  ordinal_type
  CubatureTensor<ExecSpaceType>::
  getNumPoints() const {
    ordinal_type numCubPoints = 1;
    for (auto i=0;i<numCubatures_;++i) 
      numCubPoints *= cubatures_[i]->getNumPoints();
    return numCubPoints;
  }


  template<typename ExecSpaceType>
  ordinal_type
  CubatureTensor<ExecSpaceType>::
  getDimension() const {
    return dimension_;
  }


  template<typename ExecSpaceType>
  ordinal_type
  CubatureTensor<ExecSpaceType>::
  getNumCubature() const {
    return numCubatures_;
  }


  template<typename ExecSpaceType>
  void
  CubatureTensor<ExecSpaceType>::
  getAccuracy( ordinal_type &accuracy[Parameters::MaxDimension],
               ordinal_type &numCubatures ) const { 
    numCubatures = numCubatures_;
    for (auto i=0;i<numCubatures_;++i) 
      accuracy[i] = degree_[i];
  }

} // end namespace Intrepid2
