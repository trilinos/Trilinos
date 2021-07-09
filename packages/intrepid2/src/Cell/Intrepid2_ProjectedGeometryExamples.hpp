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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov),
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_ProjectedGeometryExamples.hpp
    \brief  Example functors that can be used with Intrepid2::ProjectedGeometry.

    \author Nathan V. Roberts
*/

#ifndef Intrepid2_ProjectedGeometryExamples_h
#define Intrepid2_ProjectedGeometryExamples_h

namespace Intrepid2 {
/** \class Intrepid2::ProjectedGeometryIdentityMap
    \brief Identity map; simply preserves linear geometry.  Intended primarily for tests.
*/
  template<typename Scalar, const int spaceDim>
  class ProjectedGeometryIdentityMap
  {
  public:
    //! coordinate values
    KOKKOS_INLINE_FUNCTION
    Scalar operator()(const Kokkos::Array<Scalar,spaceDim> &coords, const ordinal_type &d) const
    {
      return coords[d];
    }
    
    //! gradient of the (identity) mapping
    KOKKOS_INLINE_FUNCTION
    Scalar operator()(const Kokkos::Array<Scalar,spaceDim> &coords, const ordinal_type &d1, const ordinal_type &d2) const
    {
      return (d1 == d2) ? 1.0 : 0.0;
    }
  };
    
/** \class Intrepid2::UnitSquareToCircle
    \brief Maps unit square [-1,1]x[-1,1] to circle of radius 1.
 
   See https://www.xarg.org/2017/07/how-to-map-a-square-to-a-circle/
 */
    template<typename Scalar>
    class UnitSquareToCircle
    {
    public:
      //! coordinate values
      KOKKOS_INLINE_FUNCTION
      Scalar operator()(const Kokkos::Array<Scalar,2> &coords, const ordinal_type &d) const
      {
        const Scalar &x = coords[0];
        const Scalar &y = coords[1];
        
        Scalar x_prime = x * sqrt(1. - y * y / 2.);
        Scalar y_prime = y * sqrt(1. - x * x / 2.);
        
        return (d==0) ? x_prime : y_prime;
      }
      
      //! gradient of the mapping
      KOKKOS_INLINE_FUNCTION
      Scalar operator()(const Kokkos::Array<Scalar,2> &coords, const ordinal_type &d1, const ordinal_type &d2) const
      {
        const Scalar &x = coords[0];
        const Scalar &y = coords[1];
        
        Scalar x_prime_dx =           sqrt(1. - y * y / 2.);
        Scalar x_prime_dy = - x * y / sqrt(1. - y * y / 2.);
        Scalar y_prime_dx = - x * y / sqrt(1. - x * x / 2.);
        Scalar y_prime_dy =           sqrt(1. - x * x / 2.);
        
        if      ((d1 == 0) && (d2 == 0)) return x_prime_dx;
        else if ((d1 == 0) && (d2 == 1)) return x_prime_dy;
        else if ((d1 == 1) && (d2 == 0)) return y_prime_dx;
        else if ((d1 == 1) && (d2 == 1)) return y_prime_dy;
        
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "Unsupported dim");
        return 0;
      }
    };

/** \class Intrepid2::UnitCubeToSphere
   \brief Maps unit cube [-1,1]x[-1,1]x[-1,1] to sphere of radius 1.

  Extends the map in UnitSquareToCircle to 3D.
*/
  template<typename Scalar>
  class UnitCubeToSphere
  {
    //! coordinate values
    KOKKOS_INLINE_FUNCTION
    Scalar operator()(const Kokkos::Array<Scalar,3> &coords, const ordinal_type &d) const
    {
      const int spaceDim = 3;
      // value is x_i * sqrt(1 - sum(x_j * x_j / 2.)) where j ≠ i
      const Scalar &x_d = coords[d];
      
      Scalar radical = 1.0;
      for (int d1=0; d1<spaceDim; d1++)
      {
        const Scalar valueToSubtract = (d1 == d) ? 0.0 : coords[d1] * coords[d1] / 2.0;
        radical -= valueToSubtract;
      }
      
      return x_d * sqrt(radical);
    }
    
    //! gradient of the mapping
    KOKKOS_INLINE_FUNCTION
    Scalar operator()(const Kokkos::Array<Scalar,3> &coords, const ordinal_type &d1, const ordinal_type &d2) const
    {
      const int spaceDim = 3;
      const Scalar &x_d1 = coords[d1];
      const Scalar &x_d2 = coords[d2];
      
      Scalar radical = 1.0;
      for (int d=0; d<spaceDim; d++)
      {
        const Scalar valueToSubtract = (d1 == d) ? 0.0 : coords[d] * coords[d] / 2.0;
        radical -= valueToSubtract;
      }
      
      Scalar weight = (d1 == d2) ? 1.0 : - x_d1 * x_d2;
      return weight * sqrt(radical);
    }
  };

}
#endif /* Intrepid2_ProjectedGeometryExamples_h */
