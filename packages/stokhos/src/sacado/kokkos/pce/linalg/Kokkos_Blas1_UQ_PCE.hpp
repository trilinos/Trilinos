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

#ifndef KOKKOS_BLAS1_UQ_PCE_HPP
#define KOKKOS_BLAS1_UQ_PCE_HPP

#include "Sacado_UQ_PCE.hpp"
#include "Kokkos_View_UQ_PCE.hpp"
#include "Kokkos_InnerProductSpaceTraits_UQ_PCE.hpp"
#include "Kokkos_Blas1_MV.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos Vector/MultiVector math functions
//----------------------------------------------------------------------------

namespace KokkosBlas {

template <typename XT, typename XL, typename XD, typename XM,
          typename YT, typename YL, typename YD, typename YM>
typename Kokkos::Details::InnerProductSpaceTraits< typename Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >::non_const_value_type >::dot_type
dot(const Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >& x,
    const Kokkos::View< YT,YL,YD,YM,Kokkos::Impl::ViewPCEContiguous >& y)
{
  typedef Kokkos::Impl::ViewPCEContiguous S;
  typedef Kokkos::View< XT,XL,XD,XM,S > XVector;
  typedef Kokkos::View< YT,YL,YD,YM,S > YVector;

  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;

  return dot( x_flat, y_flat );
}

template <typename RV,
          typename XT, typename XL, typename XD, typename XM,
          typename YT, typename YL, typename YD, typename YM>
void
dot(const RV& r,
    const Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >& x,
    const Kokkos::View< YT,YL,YD,YM,Kokkos::Impl::ViewPCEContiguous >& y)
{
  typedef Kokkos::Impl::ViewPCEContiguous S;
  typedef Kokkos::View< XT,XL,XD,XM,S > XVector;
  typedef Kokkos::View< YT,YL,YD,YM,S > YVector;

  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;

  dot( r, x_flat, y_flat );
}

template <typename XT, typename XL, typename XD, typename XM>
void
fill(const Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >& x,
     const typename Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >::non_const_value_type& val) {
  typedef Kokkos::Impl::ViewPCEContiguous S;
  typedef Kokkos::View< XT,XL,XD,XM,S > XVector;

  // Use the existing fill() implementation if we can
  if (Sacado::is_constant(val)) {
     typename XVector::flat_array_type x_flat = x;
     fill( x_flat, val.coeff(0) );
  }
  else {
    Kokkos::deep_copy(x, val);
  }
}

template <typename RV,
          typename XT, typename XL, typename XD, typename XM>
void
nrm2_squared(
  const RV& r,
  const Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >& x)
{
  typedef Kokkos::Impl::ViewPCEContiguous S;
  typedef Kokkos::View< XT,XL,XD,XM,S > XVector;

  typename XVector::flat_array_type x_flat = x;

  nrm2_squared( r, x_flat );
}

template <typename RV,
          typename XT, typename XL, typename XD, typename XM>
void
nrm1(
  const RV& r,
  const Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >& x)
{
  typedef Kokkos::Impl::ViewPCEContiguous S;
  typedef Kokkos::View< XT,XL,XD,XM,S > XVector;

  typename XVector::flat_array_type x_flat = x;

  nrm1( r, x_flat );
}

template <typename RV,
          typename XT, typename XL, typename XD, typename XM>
void
nrmInf(
  const RV& r,
  const Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >& x)
{
  typedef Kokkos::Impl::ViewPCEContiguous S;
  typedef Kokkos::View< XT,XL,XD,XM,S > XVector;

  typename XVector::flat_array_type x_flat = x;

  nrmInf( r, x_flat );
}

template <typename AV,
          typename XT, typename XL, typename XD, typename XM,
          typename BV,
          typename YT, typename YL, typename YD, typename YM>
void
axpby(const AV& a,
      const Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >& x,
      const BV& b,
      const Kokkos::View< YT,YL,YD,YM,Kokkos::Impl::ViewPCEContiguous >& y)
{
  typedef Kokkos::Impl::ViewPCEContiguous S;
  typedef Kokkos::View< XT,XL,XD,XM,S > XVector;
  typedef Kokkos::View< YT,YL,YD,YM,S > YVector;

  if (!Sacado::is_constant(a) || !Sacado::is_constant(b)) {
    Kokkos::Impl::raise_error("axpby not implemented for non-constant a or b");
  }

  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;
  auto aa = Sacado::Value<AV>::eval(a);
  auto bb = Sacado::Value<BV>::eval(b);
  axpby( aa, x_flat, bb, y_flat );
}

// Currently not handling scal() when AV is a view

template <typename RT, typename RL, typename RD, typename RM,
          typename XT, typename XL, typename XD, typename XM>
void
scal(const Kokkos::View< RT,RL,RD,RM,Kokkos::Impl::ViewPCEContiguous >& r,
     const typename Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >::non_const_value_type& a,
     const Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >& x)
{
  typedef Kokkos::Impl::ViewPCEContiguous S;
  typedef Kokkos::View< XT,XL,XD,XM,S > XVector;
  typedef Kokkos::View< RT,RL,RD,RM,S > RVector;

  if (!Sacado::is_constant(a)) {
    Kokkos::Impl::raise_error("scal not implemented for non-constant a");
  }

  typename XVector::flat_array_type x_flat = x;
  typename RVector::flat_array_type r_flat = r;
  scal( r_flat, a.coeff(0), x_flat );
}

// abs -- can't do this one by flattening.  Hold out for refactoring of scalar
// types in Kokkos

// We have a special verision of update for scalar alpha/beta/gamma since it
// is used in TrilinosCouplings CG solve (even though Tpetra doesn't).
template <typename XT, typename XL, typename XD, typename XM,
          typename YT, typename YL, typename YD, typename YM,
          typename ZT, typename ZL, typename ZD, typename ZM>
void
update(
  const typename Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >::intrinsic_scalar_type& alpha,
  const Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >& x,
  const typename Kokkos::View< YT,YL,YD,YM,Kokkos::Impl::ViewPCEContiguous >::intrinsic_scalar_type& beta,
  const Kokkos::View< YT,YL,YD,YM,Kokkos::Impl::ViewPCEContiguous >& y,
  const typename Kokkos::View< ZT,ZL,ZD,ZM,Kokkos::Impl::ViewPCEContiguous >::intrinsic_scalar_type& gamma,
  const Kokkos::View< ZT,ZL,ZD,ZM,Kokkos::Impl::ViewPCEContiguous >& z)
{
  typedef Kokkos::Impl::ViewPCEContiguous S;
  typedef Kokkos::View< XT,XL,XD,XM,S > XVector;
  typedef Kokkos::View< YT,YL,YD,YM,S > YVector;
  typedef Kokkos::View< ZT,ZL,ZD,ZM,S > ZVector;

  typename XVector::flat_array_type x_flat = x;
  typename YVector::flat_array_type y_flat = y;
  typename ZVector::flat_array_type z_flat = z;

  update( alpha, x_flat, beta, y_flat, gamma, z_flat);

}

template <typename XT, typename XL, typename XD, typename XM,
          typename YT, typename YL, typename YD, typename YM,
          typename ZT, typename ZL, typename ZD, typename ZM>
void
update(
  const typename Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >::non_const_value_type& alpha,
  const Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >& x,
  const typename Kokkos::View< YT,YL,YD,YM,Kokkos::Impl::ViewPCEContiguous >::non_const_value_type& beta,
  const Kokkos::View< YT,YL,YD,YM,Kokkos::Impl::ViewPCEContiguous >& y,
  const typename Kokkos::View< ZT,ZL,ZD,ZM,Kokkos::Impl::ViewPCEContiguous >::non_const_value_type& gamma,
  const Kokkos::View< ZT,ZL,ZD,ZM,Kokkos::Impl::ViewPCEContiguous >& z)
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

template<class RT, class RL, class RD, class RM,
         class XT, class XL, class XD, class XM,
         class SizeType>
struct MV_Reciprocal_Functor<
  Kokkos::View<RT,RL,RD,RM,Kokkos::Impl::ViewPCEContiguous>,
  Kokkos::View<XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous>,
  SizeType>
{
  typedef Kokkos::View<RT,RL,RD,RM,Kokkos::Impl::ViewPCEContiguous> RMV;
  typedef Kokkos::View<XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous> XMV;
  typedef typename RMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::intrinsic_scalar_type> ATS;

  const size_type numCols;
  RMV R_;
  XMV X_;

  MV_Reciprocal_Functor (const RMV& R, const XMV& X) :
    numCols (X.dimension_1 ()), R_ (R), X_ (X)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type j = 0; j < numCols; ++j) {
      R_(i,j) = ATS::one () / X_(i,j).fastAccessCoeff(0);
    }
  }
};

template<class RT, class RL, class RD, class RM,
         class SizeType>
struct MV_ReciprocalSelf_Functor<
  Kokkos::View<RT,RL,RD,RM,Kokkos::Impl::ViewPCEContiguous>,
  SizeType>
{
  typedef Kokkos::View<RT,RL,RD,RM,Kokkos::Impl::ViewPCEContiguous> RMV;
  typedef typename RMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename RMV::intrinsic_scalar_type> ATS;

  const size_type numCols;
  RMV R_;

  MV_ReciprocalSelf_Functor (const RMV& R) :
    numCols (R.dimension_1 ()), R_ (R)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type j = 0; j < numCols; ++j) {
      R_(i,j) = ATS::one () / R_(i,j).fastAccessCoeff(0);
    }
  }
};

template<class RT, class RL, class RD, class RM,
         class XT, class XL, class XD, class XM,
         class SizeType>
struct V_Reciprocal_Functor<
  Kokkos::View<RT,RL,RD,RM,Kokkos::Impl::ViewPCEContiguous>,
  Kokkos::View<XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous>,
  SizeType>
{
  typedef Kokkos::View<RT,RL,RD,RM,Kokkos::Impl::ViewPCEContiguous> RV;
  typedef Kokkos::View<XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous> XV;
  typedef typename RV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename XV::intrinsic_scalar_type> ATS;

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

template<class RT, class RL, class RD, class RM,
         class SizeType>
struct V_ReciprocalSelf_Functor<
  Kokkos::View<RT,RL,RD,RM,Kokkos::Impl::ViewPCEContiguous>,
  SizeType>
{
  typedef Kokkos::View<RT,RL,RD,RM,Kokkos::Impl::ViewPCEContiguous> RV;
  typedef typename RV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::intrinsic_scalar_type> ATS;

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

template <typename RT, typename RL, typename RD, typename RM,
          typename XT, typename XL, typename XD, typename XM>
void
sum(
  const Kokkos::View< RT,RL,RD,RM,Kokkos::Impl::ViewPCEContiguous >& r,
  const Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >& x)
{
  typedef Kokkos::Impl::ViewPCEContiguous S;
  typedef Kokkos::View< XT,XL,XD,XM,S > XVector;
  typedef Kokkos::View< RT,RL,RD,RM,S > RVector;

  typename XVector::flat_array_type x_flat = x;
  typename RVector::flat_array_type r_flat = r;
  sum( r_flat, x_flat );
}

template <typename RT, typename RL, typename RD, typename RM,
          typename XT, typename XL, typename XD, typename XM,
          typename WT, typename WL, typename WD, typename WM>
void
nrm2w_squared(
  const Kokkos::View< RT,RL,RD,RM,Kokkos::Impl::ViewPCEContiguous >& r,
  const Kokkos::View< XT,XL,XD,XM,Kokkos::Impl::ViewPCEContiguous >& x,
  const Kokkos::View< WT,WL,WD,WM,Kokkos::Impl::ViewPCEContiguous >& w)
{
  typedef Kokkos::Impl::ViewPCEContiguous S;
  typedef Kokkos::View< XT,XL,XD,XM,S > XVector;
  typedef Kokkos::View< RT,RL,RD,RM,S > RVector;
  typedef Kokkos::View< WT,WL,WD,WM,S > WVector;

  typename XVector::flat_array_type x_flat = x;
  typename RVector::flat_array_type r_flat = r;
  typename WVector::flat_array_type w_flat = w;
  nrm2w_squared( r_flat, x_flat, w_flat );
}

// Mean-based implementation of mult()
namespace Impl {

template<class CT, class CL, class CD, class CM,
         class AT, class AL, class AD, class AM,
         class BT, class BL, class BD, class BM,
         int scalar_ab, int scalar_c, class SizeType>
struct MV_MultFunctor<
  Kokkos::View<CT,CL,CD,CM,Kokkos::Impl::ViewPCEContiguous>,
  Kokkos::View<AT,AL,AD,AM,Kokkos::Impl::ViewPCEContiguous>,
  Kokkos::View<BT,BL,BD,BM,Kokkos::Impl::ViewPCEContiguous>,
  scalar_ab, scalar_c, SizeType>
{
  typedef Kokkos::View<CT,CL,CD,CM,Kokkos::Impl::ViewPCEContiguous> CMV;
  typedef Kokkos::View<AT,AL,AD,AM,Kokkos::Impl::ViewPCEContiguous> AV;
  typedef Kokkos::View<BT,BL,BD,BM,Kokkos::Impl::ViewPCEContiguous> BMV;
  typedef typename CMV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename CMV::intrinsic_scalar_type> ATS;

  const size_type m_n;
  const size_type m_pce;
  const typename CMV::intrinsic_scalar_type m_c;
  CMV m_C;
  const typename AV::intrinsic_scalar_type m_ab;
  AV m_A;
  BMV m_B;

  MV_MultFunctor (typename CMV::const_value_type& c,
                  const CMV& C,
                  typename AV::const_value_type& ab,
                  const AV& A,
                  const BMV& B) :
    m_n (C.dimension_1 ()),
    m_pce (C.sacado_size()),
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
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
          for (size_type l=0; l<m_pce; ++l)
            m_C(i,j).fastAccessCoeff(l) = ATS::zero ();
        }
      }
      else { // ab != 0, c == 0
        typename AV::intrinsic_scalar_type Ai = m_A(i).fastAccessCoeff(0);
        for (size_type j = 0; j < m_n; ++j) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
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
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
          for (size_type l=0; l<m_pce; ++l)
            m_C(i,j).fastAccessCoeff(l) = m_c * m_C(i,j).fastAccessCoeff(l);
        }
      }
      else { // m_ab != 0, and m_c != 0
        typename AV::intrinsic_scalar_type Ai = m_A(i).fastAccessCoeff(0);
        for (size_type j = 0; j < m_n; ++j) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
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

template<class CT, class CL, class CD, class CM,
         class AT, class AL, class AD, class AM,
         class BT, class BL, class BD, class BM,
         int scalar_ab, int scalar_c, class SizeType>
struct V_MultFunctor<
  Kokkos::View<CT,CL,CD,CM,Kokkos::Impl::ViewPCEContiguous>,
  Kokkos::View<AT,AL,AD,AM,Kokkos::Impl::ViewPCEContiguous>,
  Kokkos::View<BT,BL,BD,BM,Kokkos::Impl::ViewPCEContiguous>,
  scalar_ab, scalar_c, SizeType>
{
  typedef Kokkos::View<CT,CL,CD,CM,Kokkos::Impl::ViewPCEContiguous> CV;
  typedef Kokkos::View<AT,AL,AD,AM,Kokkos::Impl::ViewPCEContiguous> AV;
  typedef Kokkos::View<BT,BL,BD,BM,Kokkos::Impl::ViewPCEContiguous> BV;
  typedef typename CV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename CV::intrinsic_scalar_type> ATS;

  const size_type m_pce;
  const typename CV::intrinsic_scalar_type m_c;
  CV m_C;
  const typename AV::intrinsic_scalar_type m_ab;
  AV m_A;
  BV m_B;

  V_MultFunctor (typename CV::const_value_type& c,
                 const CV& C,
                 typename AV::const_value_type& ab,
                 const AV& A,
                 const BV& B) :
    m_pce (C.sacado_size()),
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
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type l=0; l<m_pce; ++l)
          m_C(i).fastAccessCoeff(l) = ATS::zero ();
      }
      else { // ab != 0, c == 0
        typename AV::intrinsic_scalar_type Ai = m_A(i).fastAccessCoeff(0);
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type l=0; l<m_pce; ++l)
          m_C(i).fastAccessCoeff(l) = m_ab * Ai * m_B(i).fastAccessCoeff(l);
      }
    } else { // c != 0
      if (scalar_ab == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type l=0; l<m_pce; ++l)
          m_C(i).fastAccessCoeff(l) = m_c * m_C(i).fastAccessCoeff(l);
      }
      else { // m_ab != 0, and m_c != 0
        typename AV::intrinsic_scalar_type Ai = m_A(i).fastAccessCoeff(0);
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
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
