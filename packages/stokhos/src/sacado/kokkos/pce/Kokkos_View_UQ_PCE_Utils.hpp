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

#ifndef KOKKOS_VIEW_UQ_PCE_UTILS_HPP
#define KOKKOS_VIEW_UQ_PCE_UTILS_HPP

#include "Kokkos_View_Utils.hpp"
#include "Kokkos_View.hpp"

namespace Sacado {
  namespace UQ {
    template <typename Storage >
    class PCE;
  }
}

namespace Kokkos {

// Type name for a local, unmanaged view with possibly a different static size
template <typename ViewType,
          unsigned LocalSize,
          unsigned Rank = ViewType::Rank,
          bool isStatic = ViewType::is_static>
struct LocalUQPCEView {};

template <typename ViewType,
          unsigned LocalSize>
struct LocalUQPCEView< ViewType, LocalSize, 1, true > {
  typedef typename ViewType::value_type vector_type;
  typedef typename ViewType::array_layout array_layout;
  typedef typename ViewType::execution_space execution_space;
  typedef typename ViewType::device_type device_type;
  typedef typename vector_type::storage_type storage_type;
  typedef typename storage_type::template apply_N<LocalSize> StorageApply;
  typedef typename StorageApply::type local_storage_type;
  typedef Sacado::UQ::PCE< local_storage_type > local_value_type;

  typedef Kokkos::View< local_value_type*,
                        array_layout,
                        device_type,
                        Kokkos::MemoryUnmanaged > type;
};

template <typename ViewType,
          unsigned LocalSize>
struct LocalUQPCEView<ViewType, LocalSize, 1, false> {
  typedef typename ViewType::value_type vector_type;
  typedef typename ViewType::array_layout array_layout;
  typedef typename ViewType::execution_space execution_space;
  typedef typename ViewType::device_type device_type;

  typedef Kokkos::View< vector_type*,
                        array_layout,
                        device_type,
                        Kokkos::MemoryUnmanaged > type;
};

namespace Impl {

template< class OldStorageType , class Device >
struct RebindStokhosStorageDevice< Sacado::UQ::PCE< OldStorageType > , Device >
{
  typedef typename
    OldStorageType::template apply<
      typename OldStorageType::ordinal_type ,
      typename OldStorageType::value_type ,
      Device >
    NewStorageApply ;

  typedef typename NewStorageApply::type NewStorageType ;
  typedef typename Sacado::UQ::PCE< OldStorageType >::template apply< NewStorageType > NewVectorApply ;

  typedef typename NewVectorApply::type type ;
};

template< class OldStorageType , class Device >
struct RebindStokhosStorageDevice< const Sacado::UQ::PCE< OldStorageType > , Device >
{
  typedef typename
    OldStorageType::template apply<
      typename OldStorageType::ordinal_type ,
      typename OldStorageType::value_type ,
      Device >
    NewStorageApply ;

  typedef typename NewStorageApply::type NewStorageType ;
  typedef typename Sacado::UQ::PCE< OldStorageType >::template apply< NewStorageType > NewVectorApply ;

  typedef const typename NewVectorApply::type type ;
};

} // namespace Impl

// Whether a given type is a view with scalar type Sacado::UQ::PCE
template <typename view_type>
struct is_view_uq_pce { static const bool value = false; };

// Typename of the Cijk tensor in a view
template <typename view_type, typename Enabled = void>
struct CijkType {};

} // namespace Kokkos

#endif // KOKKOS_UQ_PCE_UTILS_HPP
