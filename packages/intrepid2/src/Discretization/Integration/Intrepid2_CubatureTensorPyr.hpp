// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CubatureTensorPyr.hpp
    \brief  Header file for the Intrepid2::CubatureTensorPyr class.
    \author Created by P. Bochev, D. Ridzal, M. Perego. 
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_TENSOR_PYR_HPP__
#define __INTREPID2_CUBATURE_TENSOR_PYR_HPP__

#include "Intrepid2_CubatureTensor.hpp"

namespace Intrepid2 {
  
  /** \class Intrepid2::CubatureTensorPyr
      \brief Defines tensor-product cubature (integration) rules in Intrepid.
  */
  template<typename DeviceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubatureTensorPyr
    : public CubatureTensor<DeviceType,pointValueType,weightValueType> {
  public:

    template<typename cubPointViewType,
             typename cubWeightViewType>
    struct Functor {
      cubPointViewType  _cubPoints;
      cubWeightViewType _cubWeights;

      KOKKOS_INLINE_FUNCTION
      Functor( cubPointViewType  cubPoints_,
               cubWeightViewType cubWeights_)
        : _cubPoints(cubPoints_), _cubWeights(cubWeights_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type i) const {
        const auto tmp = 0.5*(1.0 - _cubPoints(i,2));
        _cubPoints(i,0) *= tmp; 
        _cubPoints(i,1) *= tmp; 
        _cubPoints(i,2)  = 1.0 - tmp; 

        _cubWeights(i)  /= 8.0;  // when line jacobi20 is used (exact)
        //_cubWeights(i)  *= (tmp*tmp*.5); // when gauss integration rule is used (need over-integration) 
      }
    };

    template<typename cubPointValueType,  class ...cubPointProperties,
             typename cubWeightValueType, class ...cubWeightProperties>
    void
    getCubatureImpl( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
                     Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights ) const;

    using typename Cubature<DeviceType,pointValueType,weightValueType>::PointViewType;
    using typename Cubature<DeviceType,pointValueType,weightValueType>::weightViewType;

    using CubatureTensor<DeviceType,pointValueType,weightValueType>::getCubature;
    using typename CubatureTensor<DeviceType,pointValueType,weightValueType>::ExecSpaceType;

    virtual
    void
    getCubature( PointViewType  cubPoints,
                 weightViewType cubWeights ) const {
      getCubatureImpl( cubPoints,
                       cubWeights );
    }

    CubatureTensorPyr()
      : CubatureTensor<DeviceType,pointValueType,weightValueType>() {}
    
    CubatureTensorPyr(const CubatureTensorPyr &b)
      : CubatureTensor<DeviceType,pointValueType,weightValueType>(b) {}

    template<typename CubatureLineType>
    CubatureTensorPyr( const CubatureLineType line ) 
      : CubatureTensor<DeviceType,pointValueType,weightValueType>(line, line, line) {}

    template<typename CubatureLineType0,
             typename CubatureLineType1,
             typename CubatureLineType2>
    CubatureTensorPyr( const CubatureLineType0 line0,
                       const CubatureLineType1 line1,
                       const CubatureLineType2 line2 ) 
      : CubatureTensor<DeviceType,pointValueType,weightValueType>(line0, line1, line2) {}
    
  };
} 

#include "Intrepid2_CubatureTensorPyrDef.hpp"

#endif
