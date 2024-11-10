// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_BLAS1_UQ_PCE_HPP
#define KOKKOS_BLAS1_UQ_PCE_HPP

#include "Sacado_UQ_PCE.hpp"
#include "Kokkos_View_UQ_PCE.hpp"
#include "Kokkos_InnerProductSpaceTraits_UQ_PCE.hpp"
#include "KokkosBlas.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos Vector/MultiVector math functions
//----------------------------------------------------------------------------

namespace KokkosBlas {

template <typename XD, typename ... XP,
          typename YD, typename ... YP>
typename std::enable_if<
  Kokkos::is_view_uq_pce< Kokkos::View<XD,XP...> >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View<YD,YP...> >::value,
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
  Kokkos::is_view_uq_pce< Kokkos::View<XD,XP...> >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View<YD,YP...> >::value >::type
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
  Kokkos::is_view_uq_pce< Kokkos::View<XD,XP...> >::value >::type
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
  Kokkos::is_view_uq_pce< Kokkos::View<XD,XP...> >::value >::type
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
  Kokkos::is_view_uq_pce< Kokkos::View<XD,XP...> >::value >::type
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
  Kokkos::is_view_uq_pce< Kokkos::View<XD,XP...> >::value >::type
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
  Kokkos::is_view_uq_pce< Kokkos::View<XD,XP...> >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View<YD,YP...> >::value >::type
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
  Kokkos::is_view_uq_pce< Kokkos::View<RD,RP...> >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View<XD,XP...> >::value >::type
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
  Kokkos::is_view_uq_pce< Kokkos::View<XD,XP...> >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View<YD,YP...> >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View<ZD,ZP...> >::value >::type
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
  Kokkos::is_view_uq_pce< Kokkos::View<XD,XP...> >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View<YD,YP...> >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View<ZD,ZP...> >::value >::type
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

// Mean-based implementation of reciprocal()
namespace Impl {

template<class RS, class ... RP,
         class XS, class ... XP,
         class SizeType>
struct MV_Reciprocal_Functor<
  Kokkos::View<Sacado::UQ::PCE<RS>**,RP...>,
  Kokkos::View<const Sacado::UQ::PCE<XS>**,XP...>,
  SizeType>
{
  typedef Kokkos::View<Sacado::UQ::PCE<RS>**,RP...> RMV;
  typedef Kokkos::View<const Sacado::UQ::PCE<XS>**,XP...> XMV;
  typedef typename RMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::ArithTraits<typename Kokkos::IntrinsicScalarType<XMV>::type> ATS;

  const size_type numCols;
  RMV R_;
  XMV X_;

  MV_Reciprocal_Functor (const RMV& R, const XMV& X) :
    numCols (X.extent(1)), R_ (R), X_ (X)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type j = 0; j < numCols; ++j) {
      R_(i,j) = ATS::one () / X_(i,j).fastAccessCoeff(0);
    }
  }
};

template<class RS, class ... RP,
         class SizeType>
struct MV_ReciprocalSelf_Functor<
  Kokkos::View<Sacado::UQ::PCE<RS>**,RP...>,
  SizeType>
{
  typedef Kokkos::View<Sacado::UQ::PCE<RS>**,RP...> RMV;
  typedef typename RMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::ArithTraits<typename Kokkos::IntrinsicScalarType<RMV>::type> ATS;

  const size_type numCols;
  RMV R_;

  MV_ReciprocalSelf_Functor (const RMV& R) :
    numCols (R.extent(1)), R_ (R)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type j = 0; j < numCols; ++j) {
      R_(i,j) = ATS::one () / R_(i,j).fastAccessCoeff(0);
    }
  }
};

template<class RS, class ... RP,
         class XS, class ... XP,
         class SizeType>
struct V_Reciprocal_Functor<
  Kokkos::View<Sacado::UQ::PCE<RS>*,RP...>,
  Kokkos::View<const Sacado::UQ::PCE<XS>*,XP...>,
  SizeType>
{
  typedef Kokkos::View<Sacado::UQ::PCE<RS>*,RP...> RV;
  typedef Kokkos::View<const Sacado::UQ::PCE<XS>*,XP...> XV;
  typedef typename RV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::ArithTraits<typename Kokkos::IntrinsicScalarType<XV>::type> ATS;

  RV R_;
  XV X_;

  V_Reciprocal_Functor (const RV& R, const XV& X) : R_ (R), X_ (X)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    R_(i) = ATS::one () / X_(i).fastAccessCoeff(0);
  }
};

template<class RS, class ... RP,
         class SizeType>
struct V_ReciprocalSelf_Functor<
  Kokkos::View<Sacado::UQ::PCE<RS>*,RP...>,
  SizeType>
{
  typedef Kokkos::View<Sacado::UQ::PCE<RS>*,RP...> RV;
  typedef typename RV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::ArithTraits<typename Kokkos::IntrinsicScalarType<RV>::type> ATS;

  RV R_;

  V_ReciprocalSelf_Functor (const RV& R) : R_ (R)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    R_(i) = ATS::one () / R_(i).fastAccessCoeff(0);
  }
};

} // namespace Impl

template <typename RD, typename ... RP,
          typename XD, typename ... XP>
typename std::enable_if<
  Kokkos::is_view_uq_pce< Kokkos::View<RD,RP...> >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View<XD,XP...> >::value >::type
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
  Kokkos::is_view_uq_pce< Kokkos::View<RD,RP...> >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View<XD,XP...> >::value &&
  Kokkos::is_view_uq_pce< Kokkos::View<WD,WP...> >::value >::type
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

// Mean-based implementation of mult()
namespace Impl {

template<class CS, class ... CP,
         class AS, class ... AP,
         class BS, class ... BP,
         int scalar_ab, int scalar_c, class SizeType>
struct MV_MultFunctor<
  Kokkos::View<Sacado::UQ::PCE<CS>**,CP...>,
  Kokkos::View<const Sacado::UQ::PCE<AS>*,AP...>,
  Kokkos::View<const Sacado::UQ::PCE<BS>**,BP...>,
  scalar_ab, scalar_c, SizeType>
{
  typedef Kokkos::View<Sacado::UQ::PCE<CS>**,CP...> CMV;
  typedef Kokkos::View<const Sacado::UQ::PCE<AS>*,AP...> AV;
  typedef Kokkos::View<const Sacado::UQ::PCE<BS>**,BP...> BMV;
  typedef typename CMV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename Kokkos::IntrinsicScalarType<CMV>::type> ATS;

  const size_type m_n;
  const size_type m_pce;
  const typename Kokkos::IntrinsicScalarType<CMV>::type m_c;
  CMV m_C;
  const typename Kokkos::IntrinsicScalarType<AV>::type m_ab;
  AV m_A;
  BMV m_B;

  MV_MultFunctor (typename CMV::const_value_type& c,
                  const CMV& C,
                  typename AV::const_value_type& ab,
                  const AV& A,
                  const BMV& B) :
    m_n (C.extent(1)),
    m_pce (dimension_scalar(C)),
    m_c (c.coeff(0)), m_C (C), m_ab (ab.coeff(0)), m_A (A), m_B (B)
  {
    if (!Sacado::is_constant(c) || !Sacado::is_constant(ab)) {
      Kokkos::Impl::raise_error("mult not implemented for non-constant c, ab");
    }
  }

  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type& i) const
  {
    if (scalar_c == 0) {
      if (scalar_ab == 0) {
        for (size_type j = 0; j < m_n; ++j) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
          for (size_type l=0; l<m_pce; ++l)
            m_C(i,j).fastAccessCoeff(l) = ATS::zero ();
        }
      }
      else { // ab != 0, c == 0
        typename Kokkos::IntrinsicScalarType<AV>::type Ai = m_A(i).fastAccessCoeff(0);
        for (size_type j = 0; j < m_n; ++j) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
          for (size_type l=0; l<m_pce; ++l)
            m_C(i,j).fastAccessCoeff(l) =
              m_ab * Ai * m_B(i,j).fastAccessCoeff(l);
        }
      }
    } else { // c != 0
      if (scalar_ab == 0) {
        for (size_type j = 0; j < m_n; ++j) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
          for (size_type l=0; l<m_pce; ++l)
            m_C(i,j).fastAccessCoeff(l) = m_c * m_C(i,j).fastAccessCoeff(l);
        }
      }
      else { // m_ab != 0, and m_c != 0
        typename Kokkos::IntrinsicScalarType<AV>::type Ai = m_A(i).fastAccessCoeff(0);
        for (size_type j = 0; j < m_n; ++j) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
          for (size_type l=0; l<m_pce; ++l)
            m_C(i,j).fastAccessCoeff(l) =
              m_c * m_C(i,j).fastAccessCoeff(l) + m_ab * Ai * m_B(i,j).fastAccessCoeff(l);
        }
      }
    }
  }
};

template<class CS, class ... CP,
         class AS, class ... AP,
         class BS, class ... BP,
         int scalar_ab, int scalar_c, class SizeType>
struct V_MultFunctor<
  Kokkos::View<Sacado::UQ::PCE<CS>*,CP...>,
  Kokkos::View<const Sacado::UQ::PCE<AS>*,AP...>,
  Kokkos::View<const Sacado::UQ::PCE<BS>*,BP...>,
  scalar_ab, scalar_c, SizeType>
{
  typedef Kokkos::View<Sacado::UQ::PCE<CS>*,CP...> CV;
  typedef Kokkos::View<const Sacado::UQ::PCE<AS>*,AP...> AV;
  typedef Kokkos::View<const Sacado::UQ::PCE<BS>*,BP...> BV;
  typedef typename CV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::ArithTraits<typename Kokkos::IntrinsicScalarType<CV>::type> ATS;

  const size_type m_pce;
  const typename Kokkos::IntrinsicScalarType<CV>::type m_c;
  CV m_C;
  const typename Kokkos::IntrinsicScalarType<AV>::type m_ab;
  AV m_A;
  BV m_B;

  V_MultFunctor (typename CV::const_value_type& c,
                 const CV& C,
                 typename AV::const_value_type& ab,
                 const AV& A,
                 const BV& B) :
    m_pce (dimension_scalar(C)),
    m_c (c.coeff(0)), m_C (C), m_ab (ab.coeff(0)), m_A (A), m_B (B)
  {
    if (!Sacado::is_constant(c) || !Sacado::is_constant(ab)) {
      Kokkos::Impl::raise_error("mult not implemented for non-constant c, ab");
    }
  }

  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type& i) const
  {
    if (scalar_c == 0) {
      if (scalar_ab == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type l=0; l<m_pce; ++l)
          m_C(i).fastAccessCoeff(l) = ATS::zero ();
      }
      else { // ab != 0, c == 0
        typename Kokkos::IntrinsicScalarType<AV>::type Ai = m_A(i).fastAccessCoeff(0);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type l=0; l<m_pce; ++l)
          m_C(i).fastAccessCoeff(l) = m_ab * Ai * m_B(i).fastAccessCoeff(l);
      }
    } else { // c != 0
      if (scalar_ab == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type l=0; l<m_pce; ++l)
          m_C(i).fastAccessCoeff(l) = m_c * m_C(i).fastAccessCoeff(l);
      }
      else { // m_ab != 0, and m_c != 0
        typename Kokkos::IntrinsicScalarType<AV>::type Ai = m_A(i).fastAccessCoeff(0);
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type l=0; l<m_pce; ++l)
          m_C(i).fastAccessCoeff(l) =
            m_c * m_C(i).fastAccessCoeff(l) + m_ab * Ai * m_B(i).fastAccessCoeff(l);
      }
    }
  }
};

} // namespace Impl

} // namespace KokkosBlas

#endif /* #ifndef KOKKOS_MV_UQ_PCE_HPP */
