/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOSBLAS1_IMPL_DOT_IMPL_HPP_
#define KOKKOSBLAS1_IMPL_DOT_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

/// \brief Functor that implements the single-vector, two-argument
///   version of KokkosBlas::dot (dot product of two vectors).
///
/// \tparam XVector Type of the first vector x; 1-D View
/// \tparam YVector Type of the second vector y; 1-D View
/// \tparam SizeType Type of the row index used in the dot product.
///   For best performance, use int instead of size_t here.
template<class AV, class XVector, class YVector, typename SizeType>
struct DotFunctor
{
  typedef typename XVector::execution_space        execution_space;
  typedef SizeType                                 size_type;
  typedef typename AV::non_const_value_type   avalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<avalue_type>  IPT;
  typedef typename IPT::dot_type                   value_type;

  XVector  m_x;
  YVector  m_y;

  DotFunctor (const XVector& x, const YVector& y) : m_x (x), m_y (y) {}

  void run(const char* label, AV result) {
    Kokkos::RangePolicy<execution_space,size_type> policy(0,m_x.extent(0));
    Kokkos::parallel_reduce(label,policy,*this,result);
  }

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &i, value_type& sum) const
  {
    Kokkos::Details::updateDot(sum, m_x(i), m_y(i)); // sum += m_x(i) * m_y(i)
  }

  KOKKOS_INLINE_FUNCTION void
  init (volatile value_type& update) const
  {
    update = Kokkos::Details::ArithTraits<value_type>::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (value_type& update,
        const value_type& source) const
  {
    update += source ;
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const
  {
    update += source ;
  }
};


} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOSBLAS1_IMPL_DOT_IMPL_HPP_
