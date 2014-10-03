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

#ifndef KOKKOS_ATOMIC_MP_VECTOR_HPP
#define KOKKOS_ATOMIC_MP_VECTOR_HPP

#include "Sacado_MP_Vector.hpp"
#include "Kokkos_Core.hpp"

//----------------------------------------------------------------------------
// Overloads of Kokkos atomic functions for Sacado::MP::Vector scalar type
//----------------------------------------------------------------------------

namespace Kokkos {

template <typename Storage1, typename Storage2>
KOKKOS_INLINE_FUNCTION
void
atomic_assign(Sacado::MP::Vector<Storage1>* const dest,
              const Sacado::MP::Vector<Storage2>& src )
{
  typedef typename Storage1::ordinal_type ordinal_type;
  typedef typename Storage1::pointer pointer1;
  pointer1 dest_c = dest->coeff();
  const ordinal_type sz = dest->size();
  if (src.hasFastAccess(sz))
    for (ordinal_type i=0; i<sz; ++i)
      atomic_exchange(dest_c+i, src.fastAccessCoeff(i));
  else
    for (ordinal_type i=0; i<sz; ++i)
      atomic_exchange(dest_c+i, src.coeff(i));
}

template <typename Storage1, typename Storage2>
KOKKOS_INLINE_FUNCTION
void
atomic_add(Sacado::MP::Vector<Storage1>* const dest,
           const Sacado::MP::Vector<Storage2>& src)
{
  typedef typename Storage1::ordinal_type ordinal_type;
  typedef typename Storage1::pointer pointer1;
  pointer1 dest_c = dest->coeff();
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

#endif /* #ifndef KOKKOS_ATOMIC_MP_VECTOR_HPP */
