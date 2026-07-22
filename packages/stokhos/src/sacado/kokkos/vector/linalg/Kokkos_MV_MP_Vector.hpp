// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_MV_MP_VECTOR_HPP
#define KOKKOS_MV_MP_VECTOR_HPP

#include "Sacado_MP_Vector.hpp"
#include "Kokkos_View_MP_Vector.hpp"
#include "Kokkos_InnerProductSpaceTraits_MP_Vector.hpp"
#include "Kokkos_Blas1_MP_Vector.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos Vector/MultiVector math functions
//----------------------------------------------------------------------------
/*
namespace Kokkos {

// Rank-1 vector add with Sacado::MP::Vector scalar type, constant a, b
template <typename RS, typename RL, typename RD, typename RM,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
Kokkos::View< Sacado::MP::Vector<RS>*, RL, RD, RM>
V_Add( const Kokkos::View< Sacado::MP::Vector<RS>*, RL, RD, RM >& r,
       const typename Sacado::MP::Vector<XS>::value_type& av,
       const Kokkos::View< Sacado::MP::Vector<XS>*, XL, XD, XM >& x,
       const typename Sacado::MP::Vector<XS>::value_type& bv,
       const Kokkos::View< Sacado::MP::Vector<YS>*, YL, YD, YM >& y,
       int n = -1)
{
  typedef Kokkos::View< Sacado::MP::Vector<RS>*, RL, RD, RM > RVector;
  typedef Kokkos::View< Sacado::MP::Vector<XS>*, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::MP::Vector<YS>*, YL, YD, YM > YVector;

  typename RVector::flat_array_type r_flat = r;
  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;
  if (n != -1) n = n * r.sacado_size();

  V_Add( r_flat, av, x_flat, bv, y_flat, n );

  return r;
}

// Rank-1 vector add with Sacado::MP::Vector scalar type, non-constant a, b
template <typename RS, typename RL, typename RD, typename RM,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
Kokkos::View< Sacado::MP::Vector<RS>*, RL, RD, RM>
V_Add( const Kokkos::View< Sacado::MP::Vector<RS>*, RL, RD, RM >& r,
       const Sacado::MP::Vector<XS>& av,
       const Kokkos::View< Sacado::MP::Vector<XS>*, XL, XD, XM >& x,
       const Sacado::MP::Vector<XS>& bv,
       const Kokkos::View< Sacado::MP::Vector<YS>*, YL, YD, YM >& y,
       int n = -1)
{
  if (Sacado::is_constant(av) && Sacado::is_constant(bv)) {
   return V_Add( r, av.fastAccessCoeff(0), x, bv.fastAccessCoeff(0), y, n );
  }
  else {
    Impl::raise_error("V_Add not implemented for non-constant a or b");
  }
  return r;
}

// Rank-2 vector add with Sacado::MP::Vector scalar type, constant a, b
template <typename RS, typename RL, typename RD, typename RM,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
Kokkos::View< Sacado::MP::Vector<RS>**, RL, RD, RM>
MV_Add( const Kokkos::View< Sacado::MP::Vector<RS>**, RL, RD, RM >& r,
        const typename Sacado::MP::Vector<XS>::value_type& av,
        const Kokkos::View< Sacado::MP::Vector<XS>**, XL, XD, XM >& x,
        const typename Sacado::MP::Vector<XS>::value_type& bv,
        const Kokkos::View< Sacado::MP::Vector<YS>**, YL, YD, YM >& y,
        int n = -1)
{
  typedef Kokkos::View< Sacado::MP::Vector<RS>**, RL, RD, RM > RVector;
  typedef Kokkos::View< Sacado::MP::Vector<XS>**, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::MP::Vector<YS>**, YL, YD, YM > YVector;

  typename RVector::flat_array_type r_flat = r;
  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;
  if (n != -1) n = n * r.sacado_size();

  MV_Add( r_flat, av, x_flat, bv, y_flat, n );

  return r;
}

// Rank-2 vector add with Sacado::MP::Vector scalar type, non-constant a, b
template <typename RS, typename RL, typename RD, typename RM,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
Kokkos::View< Sacado::MP::Vector<RS>**, RL, RD, RM>
MV_Add( const Kokkos::View< Sacado::MP::Vector<RS>**, RL, RD, RM >& r,
        const Sacado::MP::Vector<XS>& av,
        const Kokkos::View< Sacado::MP::Vector<XS>**, XL, XD, XM >& x,
        const Sacado::MP::Vector<XS>& bv,
        const Kokkos::View< Sacado::MP::Vector<YS>**, YL, YD, YM >& y,
        int n = -1)
{
  if (Sacado::is_constant(av) && Sacado::is_constant(bv)) {
    return MV_Add( r, av.fastAccessCoeff(0), x, bv.fastAccessCoeff(0), y, n );
  }
  else {
    Impl::raise_error("MV_Add not implemented for non-constant a or b");
  }
  return r;
}

// Rank-1 dot product
template <typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
typename Details::InnerProductSpaceTraits< Sacado::MP::Vector<XS> >::dot_type
V_Dot( const Kokkos::View< Sacado::MP::Vector<XS>*, XL, XD, XM >& x,
       const Kokkos::View< Sacado::MP::Vector<YS>*, YL, YD, YM >& y,
       int n = -1 )
{
  typedef Kokkos::View< Sacado::MP::Vector<XS>*, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::MP::Vector<YS>*, YL, YD, YM > YVector;

  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;
  if (n != -1) n = n * x.sacado_size();

  return V_Dot( x_flat, y_flat, n );
}

// Rank-2 dot product
template <typename rVector,
          typename XS, typename XL, typename XD, typename XM,
          typename YS, typename YL, typename YD, typename YM>
void
MV_Dot( const rVector& r,
        const Kokkos::View< Sacado::MP::Vector<XS>**, XL, XD, XM >& x,
        const Kokkos::View< Sacado::MP::Vector<YS>**, YL, YD, YM >& y,
        int n = -1 )
{
  typedef Kokkos::View< Sacado::MP::Vector<XS>**, XL, XD, XM > XVector;
  typedef Kokkos::View< Sacado::MP::Vector<YS>**, YL, YD, YM > YVector;

  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;
  if (n != -1) n = n * x.sacado_size();

  MV_Dot( r, x_flat, y_flat, n );
}

} // namespace Kokkos
*/
#endif /* #ifndef KOKKOS_MV_MP_VECTOR_HPP */
