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

/** \file   Intrepid2_CubatureTensorPyrDef.hpp
    \brief  Definition file for the Intrepid2::CubatureTensorPyr class.
    \author Created by P. Bochev, D. Ridzal and M. Perego.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_TENSOR_PYR_DEF_HPP__
#define __INTREPID2_CUBATURE_TENSOR_PYR_DEF_HPP__

namespace Intrepid2 {

  template<typename SpT, typename PT, typename WT>
  template<typename cubPointValueType,  class ...cubPointProperties,
           typename cubWeightValueType, class ...cubWeightProperties>
  void
  CubatureTensorPyr<SpT,PT,WT>::
  getCubatureImpl( Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPoints,
                   Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeights ) const {
#ifdef HAVE_INTREPID2_DEBUG
    // check size of cubPoints and cubWeights
    INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(cubPoints.extent(0))  < this->getNumPoints() ||
                                  static_cast<ordinal_type>(cubPoints.extent(1))  < this->getDimension() ||
                                  static_cast<ordinal_type>(cubWeights.extent(0)) < this->getNumPoints(), std::out_of_range,
                                  ">>> ERROR (CubatureTensor): Insufficient space allocated for cubature points or weights.");
#endif
    CubatureTensor<SpT,PT,WT>::getCubatureImpl( cubPoints, cubWeights );

    typedef Kokkos::DynRankView<cubPointValueType, cubPointProperties...>  cubPointViewType;
    typedef Kokkos::DynRankView<cubWeightValueType,cubWeightProperties...> cubWeightViewType;
    typedef typename ExecSpace<typename cubPointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

    const auto loopSize = this->getNumPoints();
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    
    typedef Functor<cubPointViewType, cubWeightViewType> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(cubPoints, cubWeights) );
  }

}

#endif
