// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_BLAS1_MP_VECTOR_HPP
#define KOKKOS_BLAS1_MP_VECTOR_HPP

#include "Sacado_MP_Vector.hpp"
#include "Kokkos_View_MP_Vector.hpp"
#include "Kokkos_InnerProductSpaceTraits_MP_Vector.hpp"
#include "Kokkos_Blas1_MV.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos Vector/MultiVector math functions
//----------------------------------------------------------------------------

namespace KokkosBlas {

// Rank-1 vector add with Sacado::MP::Vector scalar type, constant a, b
template <typename RS, typename RL, typename RD, typename RM,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
void
update(const typename Sacado::MP::Vector<XS>::value_type& av,
       const Kokkos::View< Sacado::MP::Vector<XS>*, XL, XD, XM >& x,
       const typename Sacado::MP::Vector<XS>::value_type& bv,
       const Kokkos::View< Sacado::MP::Vector<YS>*, YL, YD, YM >& y,
       const typename Sacado::MP::Vector<XS>::value_type& cv,
       const Kokkos::View< Sacado::MP::Vector<RS>*, RL, RD, RM >& r)
{
  typedef Kokkos::View< Sacado::MP::Vector<RS>*, RL, RD, RM > RVector;
  typedef Kokkos::View< Sacado::MP::Vector<XS>*, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::MP::Vector<YS>*, YL, YD, YM > YVector;

  typename RVector::flat_array_type r_flat = r;
  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;

  update( av, x_flat, bv, y_flat , cv, r_flat);

}

// Rank-1 vector add with Sacado::MP::Vector scalar type, non-constant a, b
template <typename RS, typename RL, typename RD, typename RM,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
void
update(const Sacado::MP::Vector<XS>& av,
       const Kokkos::View< Sacado::MP::Vector<XS>*, XL, XD, XM >& x,
       const Sacado::MP::Vector<XS>& bv,
       const Kokkos::View< Sacado::MP::Vector<YS>*, YL, YD, YM >& y,
       const Sacado::MP::Vector<XS>& cv,
       const Kokkos::View< Sacado::MP::Vector<RS>*, RL, RD, RM >& r)
{
  if (Sacado::is_constant(av) && Sacado::is_constant(bv)) {
   return update( av.fastAccessCoeff(0), x, bv.fastAccessCoeff(0), y, cv.fastAccessCoeff(0), r);
  }
  else {
    Kokkos::Impl::raise_error("V_Add not implemented for non-constant a or b");
  }
}

// Rank-2 vector add with Sacado::MP::Vector scalar type, constant a, b
template <typename RS, typename RL, typename RD, typename RM,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
void
update( const typename Sacado::MP::Vector<XS>::value_type& av,
        const Kokkos::View< Sacado::MP::Vector<XS>**, XL, XD, XM >& x,
        const typename Sacado::MP::Vector<XS>::value_type& bv,
        const Kokkos::View< Sacado::MP::Vector<YS>**, YL, YD, YM >& y,
        const typename Sacado::MP::Vector<XS>::value_type& cv,
        const Kokkos::View< Sacado::MP::Vector<RS>**, RL, RD, RM >& r)
{
  typedef Kokkos::View< Sacado::MP::Vector<RS>**, RL, RD, RM > RVector;
  typedef Kokkos::View< Sacado::MP::Vector<XS>**, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::MP::Vector<YS>**, YL, YD, YM > YVector;

  typename RVector::flat_array_type r_flat = r;
  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;

  update( av, x_flat, bv, y_flat, cv, r_flat );
}

// Rank-2 vector add with Sacado::MP::Vector scalar type, non-constant a, b
template <typename RS, typename RL, typename RD, typename RM,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
void
update( const Sacado::MP::Vector<XS>& av,
        const Kokkos::View< Sacado::MP::Vector<XS>**, XL, XD, XM >& x,
        const Sacado::MP::Vector<XS>& bv,
        const Kokkos::View< Sacado::MP::Vector<YS>**, YL, YD, YM >& y,
        const Sacado::MP::Vector<XS>& cv,
        const Kokkos::View< Sacado::MP::Vector<RS>**, RL, RD, RM >& r)
{
  if (Sacado::is_constant(av) && Sacado::is_constant(bv)) {
    return update( av.fastAccessCoeff(0), x, bv.fastAccessCoeff(0), y, cv.fastAccessCoeff(0), r );
  }
  else {
    Kokkos::Impl::raise_error("MV_Add not implemented for non-constant a or b");
  }
}

// Rank-1 dot product
template <typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
typename Kokkos::Details::InnerProductSpaceTraits< Sacado::MP::Vector<XS> >::dot_type
dot( const Kokkos::View< Sacado::MP::Vector<XS>*, XL, XD, XM >& x,
       const Kokkos::View< Sacado::MP::Vector<YS>*, YL, YD, YM >& y)
{
  typedef Kokkos::View< Sacado::MP::Vector<XS>*, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::MP::Vector<YS>*, YL, YD, YM > YVector;

  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;

  return dot( x_flat, y_flat );
}

// Rank-2 dot product
template <typename rVector,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
void
dot( const rVector& r,
        const Kokkos::View< Sacado::MP::Vector<XS>**, XL, XD, XM >& x,
        const Kokkos::View< Sacado::MP::Vector<YS>**, YL, YD, YM >& y)
{
  typedef Kokkos::View< Sacado::MP::Vector<XS>**, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::MP::Vector<YS>**, YL, YD, YM > YVector;

  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;

  dot( r, x_flat, y_flat );
}

} // namespace Kokkos

#endif /* #ifndef KOKKOS_MV_MP_VECTOR_HPP */
