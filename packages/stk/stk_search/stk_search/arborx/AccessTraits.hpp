// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef LOCAL_COARSE_SEARCH_ARBORX_ACESS_TRAITS_HPP
#define LOCAL_COARSE_SEARCH_ARBORX_ACESS_TRAITS_HPP

#include "Kokkos_Core.hpp"
#include "ArborX_AccessTraits.hpp"
#include "StkToArborX.hpp"
#include <type_traits>

namespace stk::search::impl {

template <typename ViewType>
struct ViewWrapperForArborXTraits
{
  static_assert(Kokkos::is_view_v<ViewType>);
  using value_type = typename ViewType::value_type;
  using memory_space = typename ViewType::memory_space;

  ViewType view;
};

template <typename ViewType>
ViewWrapperForArborXTraits<ViewType> wrap_view_for_arborx(ViewType view)
{
  return {view};
}

namespace impl2 {
template <typename T>
struct is_pair : std::false_type
{};

template <typename T, typename U>
struct is_pair<std::pair<T, U>> : std::true_type
{};
}

template <typename T>
constexpr bool is_pair_v = impl2::is_pair<std::remove_cv_t<T>>::value;



}


namespace ArborX {

template <typename ViewType>
using ViewWrapperForArborXTraits = stk::search::impl::ViewWrapperForArborXTraits<ViewType>;

// Partial template specialization supports anything that stk can convert to an ArborX box
// Use the ViewWrapperForArborXTraits so this specialization will match only be used from within stk
// (and not match specializations provided by other users of ArborX)
template <typename ViewType>
struct AccessTraits<ViewWrapperForArborXTraits<ViewType>, ArborX::PrimitivesTag>
{
  using StkShape     = typename ViewType::value_type::first_type;
  using ArborXShape  = typename stk::search::impl::StkToArborX<StkShape>::ArborXType;
  using memory_space = typename ViewType::memory_space;

  static KOKKOS_INLINE_FUNCTION std::size_t size(ViewWrapperForArborXTraits<ViewType> const& primitives)
  {
    return primitives.view.extent(0);
  }

  static KOKKOS_INLINE_FUNCTION ArborXShape get(ViewWrapperForArborXTraits<ViewType> const& primitives, std::size_t i)
  {
    StkShape stkBox;
    if constexpr (stk::search::impl::is_pair_v<typename ViewType::value_type>)
    {
      stkBox = primitives.view(i).first;
    } else
    {
      stkBox = primitives.view(i).box;
    }

  return stk::search::impl::StkToArborX<StkShape>(stkBox);
  }
};


template <typename ViewType>
struct AccessTraits<ViewWrapperForArborXTraits<ViewType>, ArborX::PredicatesTag>
{
  using StkShape     = typename ViewType::value_type::first_type;
  using ArborXShape  = typename stk::search::impl::StkToArborX<StkShape>::ArborXType;
  using ArborXPredicateWithIndex = ArborX::PredicateWithAttachment<ArborX::Intersects<ArborXShape>, int>;
  using memory_space = typename ViewType::memory_space;

  static KOKKOS_FUNCTION std::size_t size(ViewWrapperForArborXTraits<ViewType> const& predicates)
  {
    return predicates.view.extent(0);
  }

  static KOKKOS_FUNCTION ArborXPredicateWithIndex get(ViewWrapperForArborXTraits<ViewType> const & predicates, std::size_t i)
  {

    ArborXShape arborXBox;

    if constexpr (stk::search::impl::is_pair_v<typename ViewType::value_type>)
    {
      arborXBox = stk::search::impl::StkToArborX<StkShape>(predicates.view(i).first);
    } else
    {
      arborXBox = stk::search::impl::StkToArborX<StkShape>(predicates.view(i).box);
    }

    return ArborXPredicateWithIndex(arborXBox, i);
  }
};

}


#endif
