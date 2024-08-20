// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CubatureTensorPyrDef.hpp
    \brief  Definition file for the Intrepid2::CubatureTensorPyr class.
    \author Created by P. Bochev, D. Ridzal and M. Perego.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_TENSOR_PYR_DEF_HPP__
#define __INTREPID2_CUBATURE_TENSOR_PYR_DEF_HPP__

namespace Intrepid2 {

  template<typename DT, typename PT, typename WT>
  template<typename cubPointValueType,  class ...cubPointProperties,
           typename cubWeightValueType, class ...cubWeightProperties>
  void
  CubatureTensorPyr<DT,PT,WT>::
  getCubatureImpl( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
                   Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights ) const {
#ifdef HAVE_INTREPID2_DEBUG
    // check size of cubPoints and cubWeights
    INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(cubPoints.extent(0))  < this->getNumPoints() ||
                                  static_cast<ordinal_type>(cubPoints.extent(1))  < this->getDimension() ||
                                  static_cast<ordinal_type>(cubWeights.extent(0)) < this->getNumPoints(), std::out_of_range,
                                  ">>> ERROR (CubatureTensor): Insufficient space allocated for cubature points or weights.");
#endif
    CubatureTensor<DT,PT,WT>::getCubatureImpl( cubPoints, cubWeights );

    typedef Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPointViewType;
    typedef Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeightViewType;

    const auto loopSize = this->getNumPoints();
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    
    typedef Functor<cubPointViewType, cubWeightViewType> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(cubPoints, cubWeights) );
  }

}

#endif
