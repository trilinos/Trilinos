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
  template<typename ExecSpaceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubatureTensorPyr
    : public CubatureTensor<ExecSpaceType,pointValueType,weightValueType> {
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

    typedef typename Cubature<ExecSpaceType,pointValueType,weightValueType>::PointViewType  PointViewType;
    typedef typename Cubature<ExecSpaceType,pointValueType,weightValueType>::weightViewType weightViewType;

    using CubatureTensor<ExecSpaceType,pointValueType,weightValueType>::getCubature;

    virtual
    void
    getCubature( PointViewType  cubPoints,
                 weightViewType cubWeights ) const {
      getCubatureImpl( cubPoints,
                       cubWeights );
    }

    CubatureTensorPyr()
      : CubatureTensor<ExecSpaceType,pointValueType,weightValueType>() {}
    
    CubatureTensorPyr(const CubatureTensorPyr &b)
      : CubatureTensor<ExecSpaceType,pointValueType,weightValueType>(b) {}

    template<typename CubatureLineType>
    CubatureTensorPyr( const CubatureLineType line ) 
      : CubatureTensor<ExecSpaceType,pointValueType,weightValueType>(line, line, line) {}

    template<typename CubatureLineType0,
             typename CubatureLineType1,
             typename CubatureLineType2>
    CubatureTensorPyr( const CubatureLineType0 line0,
                       const CubatureLineType1 line1,
                       const CubatureLineType2 line2 ) 
      : CubatureTensor<ExecSpaceType,pointValueType,weightValueType>(line0, line1, line2) {}
    
  };
} 

#include "Intrepid2_CubatureTensorPyrDef.hpp"

#endif
