// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Intrepid2_OrientationTools.hpp"

#include <cmath>

namespace Intrepid2 {

  namespace Test {

    template<typename cubWeightViewType>
    typename cubWeightViewType::value_type
    computeRefVolume(const ordinal_type      numPoints,
                     const cubWeightViewType cubWeights) {
      typename cubWeightViewType::value_type r_val;
      Kokkos::parallel_reduce("computeRefVolume",
                              Kokkos::RangePolicy<typename cubWeightViewType::execution_space>(0, numPoints),
      KOKKOS_LAMBDA (const int& i, double& lsum) {
          lsum += cubWeights(i);
      }, r_val);
      Kokkos::fence();

      return r_val;
    }

    template<typename ValueType,
             typename CubPointViewType,
             typename CubWeightViewType,
             size_t spaceDim>
    ValueType computeIntegralOfMonomial(CubPointViewType cubPoints,
                                        CubWeightViewType cubWeights,
                                        Kokkos::Array<int,spaceDim> degrees) {
      ValueType result = 0.0;
      
      const auto numPoints = cubPoints.extent_int(0); // P,D
      const int  dim       = spaceDim;
      
      using ExecutionSpace = typename CubWeightViewType::execution_space;
      
      Kokkos::parallel_reduce("computeIntegralOfMonomial",
                              Kokkos::RangePolicy(ExecutionSpace(), 0, numPoints),
      KOKKOS_LAMBDA (const int& pointOrdinal, double& localSum )
      {
        ValueType value = 1.0;
        for (int d=0; d<dim; d++)
        {
          auto x_d = cubPoints(pointOrdinal,d);
          value *= pow(x_d, degrees[d]);
        }
        localSum += value * cubWeights(pointOrdinal);
      }, result);
      
      return result;
    }


    /* 
      Integral of monomial over unit triangle from initial implementation by John Burkardt
       Integral ( over unit triangle ) x^m y^n dx dy = m! * n! / ( m + n + 2 )!
    */
    template<typename ValueType>
    ValueType analyticIntegralOfMonomialOverTri(const int xDeg, const int yDeg) {
      ValueType value(1.0);
      ValueType k=xDeg;
      for (int i=1; i<=yDeg; ++i) {
        value *= i/++k;
      }
      value /= ++k;
      return value/++k;
    }


    /* 
      Integral of monomial over unit tetraahedron from initial implemetation by John Burkardt
        Integral ( over unit tetraahedron ) x^m y^n z^p dx dy dz = m! * n! * p! / ( m + n + p + 3 )!
    */
    template<typename ValueType>
    ValueType analyticIntegralOfMonomialOverTet(const int xDeg, const int yDeg, const int zDeg) {
      ValueType value(1.0);
      ValueType k=xDeg;
      for (int i=1; i<=yDeg; ++i) {
        value *= i/++k;
      }
      for (int i=1; i<=zDeg; ++i) {
        value *= i/++k;
      }
      value /= ++k;
      value /= ++k;
      return value/++k;
    }

    // beta implemented here because apparently Apple clang does not define it.
    template<typename ArithmeticType>
    double beta(const ArithmeticType &x, const ArithmeticType &y)
    {
      auto Gamma_x   = std::tgamma(x);
      auto Gamma_y   = std::tgamma(y);
      auto Gamma_xpy = std::tgamma(x+y);
      return Gamma_x * Gamma_y / Gamma_xpy;
    }
  
    /* 
      Integral of monomial over unit pyramid from initial implementation by John Burkardt
        Integral ( over unit pyramid ) x^m y^n z^p dx dy dz
    */
    template<typename ValueType>
    ValueType analyticIntegralOfMonomialOverPyr(const int xDeg, const int yDeg, const int zDeg) {      
      ValueType value = 0.0;
      if (( xDeg % 2 == 0 ) && ( yDeg % 2 == 0 )) {
        int degCoeff = 2 + xDeg + yDeg;
        ValueType sign = 1.0;
        for (int i=0; i <= degCoeff; ++i) {
          auto comb = std::round(sign / ((degCoeff+1) * beta(degCoeff-i+1,i+1)));   // sign * nCr(degCoeff,i)
          value += comb / ( i + zDeg + 1.0 );
          sign = -sign;
        }
        value *= 4.0 / ( xDeg + 1.0 ) / ( yDeg + 1.0 );
      }
      return value;
    }

    /* 
       Checks whether quadratures are invariant to orientation, that is,
       if Q={(x_i,w_i)} is the set of points and weights of the quadrature rule, and 
       x_i gets mapped into X_i by a change in orientation, then there is a pair (x_j, w_j) in Q such that
       (X_i, w_i) = (x_j, w_j). 
    */
    template<typename ValueType,
             typename CubatureType,
             typename CubPointViewType,
             typename CubWeightViewType>
    bool IsQuadratureInvariantToOrientation(CubatureType cub,
                                            CubPointViewType cubPoints,
                                            CubWeightViewType cubWeights,
                                            const unsigned cellTopoKey) {

      
      ordinal_type numOrts = -1;
      switch (cellTopoKey)
      {
      case shards::Line<2>::key:
        numOrts = 2;
        break;
      
      case shards::Triangle<3>::key:
        numOrts = 6;
        break;

      case shards::Quadrilateral<4>::key:
        numOrts = 8;
        break;
      
      default:
        break;
      } 

      bool r_val = true;
      const ordinal_type npts = cub.getNumPoints();
      const ordinal_type dim = cub.getDimension();
      using DynRankViewHost = Kokkos::DynRankView<ValueType, Kokkos::HostSpace::execution_space>;
      DynRankViewHost cubPointsOriented("cubPointsOriented", npts,dim);
      
      // these constructors do any necessary allocation and copying to host
      using WeightValueType = typename CubWeightViewType::value_type;
      using PointValueType  = typename CubPointViewType::value_type;
      using HostDevice = Kokkos::HostSpace::device_type;
      TensorPoints<PointValueType,HostDevice> cubPointsHost(cubPoints);
      TensorData<WeightValueType,HostDevice> cubWeightsHost(cubWeights);
      
      for (ordinal_type ort=0;ort<numOrts;++ort) {

        Intrepid2::Impl::OrientationTools::mapToModifiedReference(cubPointsOriented, cubPointsHost, cellTopoKey, ort);
        
        for (ordinal_type i=0;i<npts;++i) {
          bool found = false;
          bool weightsMatch = false;
          for (ordinal_type j=0;j<npts;++j) {
            auto norm2 = std::pow(cubPointsOriented(j,0)-cubPointsHost(i,0),2);
            for(ordinal_type d=1; d<dim; ++d)
              norm2 += std::pow(cubPointsOriented(j,d)-cubPointsHost(i,d),2);
            norm2 = std::sqrt(norm2);
            //std::cout << "\ni: " << i << ", j: " << j << ", ort: " << ort << ", norm2: " << norm2 << " weights diff: " << std::abs(cubWeightsHost(i) - cubWeightsHost(j));
            if (norm2 < 1e-14) {
              found = true;
              weightsMatch = std::abs(cubWeightsHost(i) - cubWeightsHost(j))<1e-14;
              //if(!weightsMatch)
              //  std::cout << "\ni: " << i << ", j: " << j << ", ort: " << ort << ", norm2: " << norm2 << " weights diff: " << std::abs(cubWeightsHost(i) - cubWeightsHost(j));
              break;
            }
          }

          if(!(found && weightsMatch)) {
           // std::cout << "\ni: " << i << ", ort: " << ort << ", found: " << found << ", weightsMatch: " << weightsMatch;
            r_val = false;
            break;
          }
        }
        if(r_val == false)
          break;
      }

      return r_val;
    }

  }
}

#endif
