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

#include "Sacado_ConfigDefs.h"

#include "Sacado_MP_Vector.hpp"
#include "Kokkos_View_MP_Vector.hpp"
#include "Kokkos_InnerProductSpaceTraits_MP_Vector.hpp"
#include "KokkosBlas.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos Vector/MultiVector math functions
//----------------------------------------------------------------------------

namespace KokkosBlas {

#if defined(HAVE_STOKHOS_ENSEMBLE_REDUCT)

template <typename XD, typename ... XP,
          typename YD, typename ... YP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<YD,YP...> >::value,
  typename Kokkos::Details::InnerProductSpaceTraits<
    typename Kokkos::View<XD,XP...>::non_const_value_type >::dot_type
  >::type
dot(const Kokkos::View<XD,XP...>& x,
    const Kokkos::View<YD,YP...>& y)
{
  typedef Kokkos::View<XD,XP...> XVector;
  typedef Kokkos::View<YD,YP...> YVector;

  typename Kokkos::FlatArrayType<XVector>::type x_flat = x;
  typename Kokkos::FlatArrayType<YVector>::type y_flat = y;

  return dot( x_flat, y_flat );
}

template <typename RV,
          typename XD, typename ... XP,
          typename YD, typename ... YP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<YD,YP...> >::value >::type
dot(const RV& r,
    const Kokkos::View<XD,XP...>& x,
    const Kokkos::View<YD,YP...>& y)
{
  typedef Kokkos::View<XD,XP...> XVector;
  typedef Kokkos::View<YD,YP...> YVector;

  typename Kokkos::FlatArrayType<XVector>::type x_flat = x;
  typename Kokkos::FlatArrayType<YVector>::type y_flat = y;

  dot( r, x_flat, y_flat );
}

template <typename XD, typename ... XP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value >::type
fill(const Kokkos::View<XD,XP...>& x,
     const typename Kokkos::View<XD,XP...>::non_const_value_type& val) {
  typedef Kokkos::View<XD,XP...> XVector;

  // Use the existing fill() implementation if we can
  if (Sacado::is_constant(val)) {
     typename Kokkos::FlatArrayType<XVector>::type x_flat = x;
     fill( x_flat, val.coeff(0) );
  }
  else {
    Kokkos::deep_copy(x, val);
  }
}

template <typename RV,
          typename XD, typename ... XP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value >::type
nrm2_squared(
  const RV& r,
  const Kokkos::View<XD,XP...>& x)
{
  typedef Kokkos::View<XD,XP...> XVector;

  typename Kokkos::FlatArrayType<XVector>::type x_flat = x;

  nrm2_squared( r, x_flat );
}

template <typename RV,
          typename XD, typename ... XP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value >::type
nrm1(
  const RV& r,
  const Kokkos::View<XD,XP...>& x)
{
  typedef Kokkos::View<XD,XP...> XVector;

  typename Kokkos::FlatArrayType<XVector>::type x_flat = x;

  nrm1( r, x_flat );
}

template <typename RV,
          typename XD, typename ... XP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value >::type
nrmInf(
  const RV& r,
  const Kokkos::View<XD,XP...>& x)
{
  typedef Kokkos::View<XD,XP...> XVector;

  typename Kokkos::FlatArrayType<XVector>::type x_flat = x;

  nrmInf( r, x_flat );
}

template <typename AV,
          typename XD, typename ... XP,
          typename BV,
          typename YD, typename ... YP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<YD,YP...> >::value >::type
axpby(const AV& a,
      const Kokkos::View<XD,XP...>& x,
      const BV& b,
      const Kokkos::View<YD,YP...>& y)
{
  typedef Kokkos::View<XD,XP...> XVector;
  typedef Kokkos::View<YD,YP...> YVector;

  if (!Sacado::is_constant(a) || !Sacado::is_constant(b)) {
    Kokkos::Impl::raise_error("axpby not implemented for non-constant a or b");
  }

  typename Kokkos::FlatArrayType<XVector>::type x_flat = x;
  typename Kokkos::FlatArrayType<YVector>::type y_flat = y;
  auto aa = Sacado::Value<AV>::eval(a);
  auto bb = Sacado::Value<BV>::eval(b);
  axpby( aa, x_flat, bb, y_flat );
}

// Currently not handling scal() when AV is a view

template <typename RD, typename ... RP,
          typename XD, typename ... XP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<RD,RP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value >::type
scal(const Kokkos::View<RD,RP...>& r,
     const typename Kokkos::View<XD,XP...>::non_const_value_type& a,
     const Kokkos::View<XD,XP...>& x)
{
  typedef Kokkos::View<RD,RP...> RVector;
  typedef Kokkos::View<XD,XP...> XVector;

  if (!Sacado::is_constant(a)) {
    Kokkos::Impl::raise_error("scal not implemented for non-constant a");
  }

  typename Kokkos::FlatArrayType<XVector>::type x_flat = x;
  typename Kokkos::FlatArrayType<RVector>::type r_flat = r;
  scal( r_flat, a.coeff(0), x_flat );
}

// abs -- can't do this one by flattening.  Hold out for refactoring of scalar
// types in Kokkos

// We have a special verision of update for scalar alpha/beta/gamma since it
// is used in TrilinosCouplings CG solve (even though Tpetra doesn't).
template <typename XD, typename ... XP,
          typename YD, typename ... YP,
          typename ZD, typename ... ZP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<YD,YP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<ZD,ZP...> >::value >::type
update(
  const typename Kokkos::View<XD,XP...>::array_type::non_const_value_type& alpha,
  const Kokkos::View<XD,XP...>& x,
  const typename Kokkos::View<YD,YP...>::array_type::non_const_value_type& beta,
  const Kokkos::View<YD,YP...>& y,
  const typename Kokkos::View<ZD,ZP...>::array_type::non_const_value_type& gamma,
  const Kokkos::View<ZD,ZP...>& z)
{
  typedef Kokkos::View<XD,XP...> XVector;
  typedef Kokkos::View<YD,YP...> YVector;
  typedef Kokkos::View<ZD,ZP...> ZVector;

  typename Kokkos::FlatArrayType<XVector>::type x_flat = x;
  typename Kokkos::FlatArrayType<YVector>::type y_flat = y;
  typename Kokkos::FlatArrayType<ZVector>::type z_flat = z;

  update( alpha, x_flat, beta, y_flat, gamma, z_flat);

}

template <typename XD, typename ... XP,
          typename YD, typename ... YP,
          typename ZD, typename ... ZP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<YD,YP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<ZD,ZP...> >::value >::type
update(
  const typename Kokkos::View<XD,XP...>::non_const_value_type& alpha,
  const Kokkos::View<XD,XP...>& x,
  const typename Kokkos::View<YD,YP...>::non_const_value_type& beta,
  const Kokkos::View<YD,YP...>& y,
  const typename Kokkos::View<ZD,ZP...>::non_const_value_type& gamma,
  const Kokkos::View<ZD,ZP...>& z)
{
  if (!Sacado::is_constant(alpha) || !Sacado::is_constant(beta) ||
      !Sacado::is_constant(gamma)) {
     Kokkos::Impl::raise_error(
       "update not implemented for non-constant alpha, beta, gamma");
  }

  update( alpha.coeff(0), x, beta.coeff(0), y, gamma.coeff(0), z );
}

template <typename RD, typename ... RP,
          typename XD, typename ... XP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<RD,RP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value >::type
reciprocal(
  const Kokkos::View<RD,RP...>& r,
  const Kokkos::View<XD,XP...>& x)
{
  typedef Kokkos::View<RD,RP...> RVector;
  typedef Kokkos::View<XD,XP...> XVector;

  typename Kokkos::FlatArrayType<XVector>::type x_flat = x;
  typename Kokkos::FlatArrayType<RVector>::type r_flat = r;
  reciprocal( r_flat, x_flat );
}

template <typename RD, typename ... RP,
          typename XD, typename ... XP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<RD,RP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value >::type
sum(
  const Kokkos::View<RD,RP...>& r,
  const Kokkos::View<XD,XP...>& x)
{
  typedef Kokkos::View<RD,RP...> RVector;
  typedef Kokkos::View<XD,XP...> XVector;

  typename Kokkos::FlatArrayType<XVector>::type x_flat = x;
  typename Kokkos::FlatArrayType<RVector>::type r_flat = r;
  sum( r_flat, x_flat );
}

template <typename RD, typename ... RP,
          typename XD, typename ... XP,
          typename WD, typename ... WP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<RD,RP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<XD,XP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<WD,WP...> >::value >::type
nrm2w_squared(
  const Kokkos::View<RD,RP...>& r,
  const Kokkos::View<XD,XP...>& x,
  const Kokkos::View<WD,WP...>& w)
{
  typedef Kokkos::View<RD,RP...> RVector;
  typedef Kokkos::View<XD,XP...> XVector;
  typedef Kokkos::View<WD,WP...> WVector;

  typename Kokkos::FlatArrayType<XVector>::type x_flat = x;
  typename Kokkos::FlatArrayType<RVector>::type r_flat = r;
  typename Kokkos::FlatArrayType<WVector>::type w_flat = w;
  nrm2w_squared( r_flat, x_flat, w_flat );
}

template <typename CD, typename ... CP,
          typename AD, typename ... AP,
          typename BD, typename ... BP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<CD,CP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<AD,AP...> >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View<BD,BP...> >::value >::type
mult(
  const typename Kokkos::View<CD,CP...>::const_value_type& c,
  const Kokkos::View<CD,CP...>& C,
  const typename Kokkos::View<AD,AP...>::const_value_type& ab,
  const Kokkos::View<AD,AP...>& A,
  const Kokkos::View<BD,BP...>& B)
{
  if (!Sacado::is_constant(c) || !Sacado::is_constant(ab)) {
     Kokkos::Impl::raise_error("mult not implemented for non-constant c, ab");
  }

  typedef Kokkos::View<CD,CP...> CVector;
  typedef Kokkos::View<AD,AP...> AVector;
  typedef Kokkos::View<BD,BP...> BVector;

  typename Kokkos::FlatArrayType<CVector>::type C_flat = C;
  typename Kokkos::FlatArrayType<AVector>::type A_flat = A;
  typename Kokkos::FlatArrayType<BVector>::type B_flat = B;
  mult( c.coeff(0), C_flat, ab.coeff(0), A_flat, B_flat );
}

#endif

} // namespace KokkosBlas

#endif /* #ifndef KOKKOS_MV_MP_VECTOR_HPP */
