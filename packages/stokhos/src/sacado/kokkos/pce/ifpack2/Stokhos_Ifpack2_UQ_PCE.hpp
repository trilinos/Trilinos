// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_IFPACK2_UQ_PCE_HPP
#define STOKHOS_IFPACK2_UQ_PCE_HPP

// This header file should be included whenever compiling any Ifpack2
// code with Stokhos scalar types

// MP includes and specializations
#include "Stokhos_Tpetra_UQ_PCE.hpp"

// Specialization of LocalReciprocalThreshold functor in
// Ifpack2_Details_Chebyshev_def.hpp
namespace Ifpack2 {
namespace Details {

template <typename XV, class SizeType>
struct V_ReciprocalThresholdSelfFunctor;

// Mean-based implementation of ReciprocalThreshold
template <typename S, typename ... P, class SizeType>
struct V_ReciprocalThresholdSelfFunctor<
  Kokkos::View< Sacado::UQ::PCE<S>*,P... >,
  SizeType >
{
  typedef Kokkos::View< Sacado::UQ::PCE<S>*,P... > XVector;
  typedef typename XVector::array_type array_type;

  typedef typename array_type::execution_space   execution_space;
  typedef SizeType                                     size_type;
  typedef typename array_type::non_const_value_type   value_type;
  typedef Kokkos::ArithTraits<value_type>                    KAT;
  typedef typename KAT::mag_type                        mag_type;

  const XVector    m_x;
  const value_type m_min_val;
  const value_type m_min_val_mag;
  const size_type  m_n_pce;

  V_ReciprocalThresholdSelfFunctor(
    const XVector& x,
    const typename XVector::non_const_value_type& min_val) :
    m_x(x),
    m_min_val(min_val.fastAccessCoeff(0)),
    m_min_val_mag(KAT::abs(m_min_val)),
    m_n_pce(Kokkos::dimension_scalar(x)) {}
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i) const
  {
    value_type z = m_x(i).fastAccessCoeff(0);
    if (KAT::abs(z) < m_min_val_mag)
      z = m_min_val;
    else
      z = KAT::one() / z;
    m_x(i) = z;
  }
};

template<class XV, class SizeType>
struct LocalReciprocalThreshold;

template<class S, class ...P, class SizeType>
struct LocalReciprocalThreshold<
  Kokkos::View< Sacado::UQ::PCE<S>*,P... >, SizeType > {
  typedef Kokkos::View< Sacado::UQ::PCE<S>*,P... > XV;

  static void
  compute (const XV& X,
           const typename XV::non_const_value_type& minVal)
  {
    if (!Sacado::is_constant(minVal)) {
      Kokkos::Impl::raise_error(
        "LocalReciprocalThreshold not implemented for non-constant minVal");
    }

    if (Kokkos::dimension_scalar(X) == 1) {
      typedef typename Kokkos::FlatArrayType<XV>::type Flat_XV;
      Flat_XV flat_X = X;
      LocalReciprocalThreshold< Flat_XV, SizeType>::compute( flat_X,
                                                             minVal.coeff(0) );
    }
    else {
      V_ReciprocalThresholdSelfFunctor<XV, SizeType> op (X, minVal);
      Kokkos::parallel_for( X.extent(0), op );
    }
  }
};

}
}

#endif // STOKHOS_IFPACK2_UQ_PCE_HPP
