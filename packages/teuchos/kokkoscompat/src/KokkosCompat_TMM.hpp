/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSCOMPAT_TMM_HPP
#define KOKKOSCOMPAT_TMM_HPP

/// \file KokkosCompat_TMM.hpp
/// \brief Compatibility between Kokkos Array View and Teuchos memory management classes.
///
/// \warning This file and everything in it should be considered an
///   implementation detail of Tpetra.  We make no promises of
///   backwards compatibility for any interface in this file, nor do
///   we even promise that this header file will continue to exist.

#include <TeuchosKokkosCompat_config.h>
#include <Kokkos_Core.hpp>
#include <Teuchos_ArrayView.hpp>

#if 0

namespace Kokkos {
namespace Impl {

  // template<class DestViewType>
  // struct ViewAssignment<DestViewType, Teuchos::ArrayView<typename DestViewType::value_type> {}; // ???

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  template<class ValueType>
  struct ViewAssignment<View<ValueType, Threads>, Teuchos::ArrayView<const ValueType> > {
    ViewAssignment (View<ValueType, Threads>& dst,
                    const Teuchos::ArrayView<const ValueType>& src,
                    typename enable_if< (
                      is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout ,
                               typename ViewTraits<ST,SL,SD,SM>::array_layout >::value
                      &&
                      ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) == unsigned(ViewTraits<ST,SL,SD,SM>::rank) )
                      &&
                      ( unsigned(ViewTraits<DT,DL,DD,DM>::rank_dynamic) >= unsigned(ViewTraits<ST,SL,SD,SM>::rank_dynamic) )
                    ) >::type * = 0 )



      View<ValueType, LayoutRight, Threads, MemoryUnmanaged> srcView; // (src.getRawPtr (), src.size ());
      srcView.m_shape.N0 = src.size ();


    }

    static void deep_copy (const DestViewType& dst, const Teuchos::ArrayView<const ValueType>& src);
  };


} // namespace Impl
} // namespace Kokkos

namespace Kokkos {
  namespace Compat {

    template<class ViewType>
    inline void
    deep_copy (const ViewType& dst, const Teuchos::ArrayView<DT>& src);

    template<class ValueType>
    inline void
    deep_copy (const Kokkos::View<ValueType, Kokkos::Threads>& dst,
               const Teuchos::ArrayView<ValueType>& src)
    {
      using Kokkos::LayoutRight;
      using Kokkos::Threads;
      using Kokkos::MemoryUnmanaged;
      using Kokkos::View;

      View<ValueType, LayoutRight, Threads, MemoryUnmanaged> srcView (src.getRawPtr (), src.size ());
      Kokkos::deep_copy (dst, srcView);
    }

  } // namespace Compat
} // namespace Kokkos

#endif // 0

#endif // KOKKOSCOMPAT_TMM_HPP
