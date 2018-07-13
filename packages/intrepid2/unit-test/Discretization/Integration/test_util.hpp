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


/** \file
    \brief  Test utils (compute volume, monomial, integral)
    \author Kyungjoo Kim
*/

#ifndef __INTREPID2_TEST_UTILS_HPP__
#define __INTREPID2_TEST_UTILS_HPP__

#include "Intrepid2_config.h"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

namespace Intrepid2 {

  namespace Test {

    template<typename cubWeightViewType>
    typename cubWeightViewType::value_type
    computeRefVolume(const ordinal_type      numPoints,
                     const cubWeightViewType cubWeights) {
      typename cubWeightViewType::value_type r_val = 0.0;
      for (auto i=0;i<numPoints;++i)
        r_val += cubWeights(i);

      return r_val;
    }

    // Monomial evaluation.
    // in 1D, for point p(x)    : x^xDeg
    // in 2D, for point p(x,y)  : x^xDeg * y^yDeg
    // in 3D, for point p(x,y,z): x^xDeg * y^yDeg * z^zDeg
    template<typename ValueType, 
             typename PointViewType>
    ValueType computeMonomial(PointViewType p, 
                              const ordinal_type xDeg, 
                              const ordinal_type yDeg = 0, 
                              const ordinal_type zDeg = 0) {
      ValueType r_val = 1.0;
      const ordinal_type polydeg[3] = { xDeg, yDeg, zDeg };

      const auto dim = p.extent(0);
      for (size_type i=0;i<dim;++i) 
        r_val *= std::pow(p(i),polydeg[i]);
      
      return r_val;
    }

    // Computes integrals of monomials 
    template<typename ValueType,
             typename cubatureType,
             typename cubPointViewType,
             typename cubWeightViewType>
    ValueType computeIntegralOfMonomial(cubatureType cub,
                                        cubPointViewType cubPoints,
                                        cubWeightViewType cubWeights,
                                        const ordinal_type xDeg,
                                        const ordinal_type yDeg = 0,
                                        const ordinal_type zDeg = 0) {
      ValueType r_val = 0.0;

      // get cubature 
      cub.getCubature(cubPoints, cubWeights);

      const auto dim  = cub.getDimension();
      const auto npts = cub.getNumPoints();
      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

      for (auto i=0;i<npts;++i) {
        const auto pt = Kokkos::subdynrankview(cubPoints, i, range_type(0, dim));
        r_val += computeMonomial<ValueType>(pt, xDeg, yDeg, zDeg)*cubWeights(i);
      }

      return r_val;
    }
  }
}

#endif
