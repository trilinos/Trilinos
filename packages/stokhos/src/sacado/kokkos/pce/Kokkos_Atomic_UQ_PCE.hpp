// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_ATOMIC_UQ_ATOMIC_HPP
#define KOKKOS_ATOMIC_UQ_ATOMIC_HPP

#include "Sacado_UQ_PCE.hpp"
#include "Kokkos_Core.hpp"

//----------------------------------------------------------------------------
// Overloads of Kokkos atomic functions for Sacado::UQ::PCE scalar type
//----------------------------------------------------------------------------

namespace Kokkos {

template <typename Storage>
KOKKOS_INLINE_FUNCTION
void
atomic_assign(volatile Sacado::UQ::PCE<Storage>* const dest,
              const Sacado::UQ::PCE<Storage>& src )
{
  typedef typename Storage::ordinal_type ordinal_type;
  typedef typename Storage::volatile_pointer pointer;
  pointer dest_c = dest->coeff();
  const ordinal_type sz = dest->size();
  if (src.hasFastAccess(sz))
    for (ordinal_type i=0; i<sz; ++i)
      atomic_exchange(dest_c+i, src.fastAccessCoeff(i));
  else
    for (ordinal_type i=0; i<sz; ++i)
      atomic_exchange(dest_c+i, src.coeff(i));
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
void
atomic_add(volatile Sacado::UQ::PCE<Storage>* const dest,
           const Sacado::UQ::PCE<Storage>& src)
{
  typedef typename Storage::ordinal_type ordinal_type;
  typedef typename Storage::volatile_pointer pointer;
  pointer dest_c = dest->coeff();
  const ordinal_type sz = dest->size();
  if (src.hasFastAccess(sz))
    for (ordinal_type i=0; i<sz; ++i)
      atomic_add(dest_c+i, src.fastAccessCoeff(i));
  else
    for (ordinal_type i=0; i<sz; ++i)
      atomic_add(dest_c+i, src.coeff(i));
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_ATOMIC_UQ_PCE_HPP */
