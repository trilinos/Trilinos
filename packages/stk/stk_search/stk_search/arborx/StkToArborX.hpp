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

#ifndef STK_TO_ARBORX_HPP
#define STK_TO_ARBORX_HPP

#include <tuple>
#include <vector>
#include <type_traits>
#include <limits>
#include <math.h>

#include "ArborX.hpp"
#include "stk_search/BoxIdent.hpp"
#include "stk_search/CommonSearchUtil.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_StdAlgorithms.hpp"

namespace stk::search
{
namespace impl
{
template <typename T, typename = void>
struct StkToArborX {
};

template <typename T>
struct StkToArborX<T, std::enable_if_t<is_stk_point<T>>> {
  using ValueType = typename T::value_type;
  using ArborXType = ArborX::Point;
  using StkType = Point<ValueType>;

  KOKKOS_INLINE_FUNCTION
  StkToArborX() = delete;

  KOKKOS_INLINE_FUNCTION
  StkToArborX(StkType const& src) : data(convert(src)) {}

  KOKKOS_INLINE_FUNCTION
  operator ArborXType() { return data; }

  KOKKOS_INLINE_FUNCTION
  auto convert(StkType const& src)
  {
    // ArborX only accepts single precision floating points
    return ArborXType(static_cast<float>(src[0]), static_cast<float>(src[1]), static_cast<float>(src[2]));
  }

  ArborXType data;
};

template <typename T>
struct StkToArborX<T, std::enable_if_t<is_stk_box<T>>> {
  using ValueType = typename T::value_type;
  using ArborXType = ArborX::Box;
  using StkType = Box<ValueType>;
  using PointType = StkToArborX<Point<ValueType>>;

  KOKKOS_INLINE_FUNCTION
  StkToArborX() = delete;

  KOKKOS_INLINE_FUNCTION
  StkToArborX(StkType const& src) : data(convert(src)) {}

  KOKKOS_INLINE_FUNCTION
  operator ArborXType() { return data; }

  KOKKOS_INLINE_FUNCTION
  auto convert(StkType const& src)
  {
    if (std::is_same_v<ValueType, float>) {
      return ArborXType(ArborX::Point(src.get_x_min(), src.get_y_min(), src.get_z_min()),
                        ArborX::Point(src.get_x_max(), src.get_y_max(), src.get_z_max()));
    }
    else {
      return ArborXType(ArborX::Point(src.get_expanded_x_min(), src.get_expanded_y_min(), src.get_expanded_z_min()),
                        ArborX::Point(src.get_expanded_x_max(), src.get_expanded_y_max(), src.get_expanded_z_max()));
    }
  }

  ArborXType data;
};

template <typename T>
struct StkToArborX<T, std::enable_if_t<is_stk_sphere<T>>> {
  using ValueType = typename T::value_type;
  using ArborXType = ArborX::Box;
  using StkType = Sphere<ValueType>;
  using PointType = StkToArborX<Point<ValueType>>;

  KOKKOS_INLINE_FUNCTION
  StkToArborX() = delete;

  KOKKOS_INLINE_FUNCTION
  StkToArborX(StkType const& src) : data(convert(src)) {}

  KOKKOS_INLINE_FUNCTION
  operator ArborXType() { return data; }

  // ArborX requires primitives to decay to either a Point or a Box
  KOKKOS_INLINE_FUNCTION
  auto convert(StkType const& src)
  {

    if (std::is_same_v<ValueType, float>) {
      return ArborXType(ArborX::Point(src.get_x_min(), src.get_y_min(), src.get_z_min()),
                        ArborX::Point(src.get_x_max(), src.get_y_max(), src.get_z_max()));
    }
    else {
      return ArborXType(ArborX::Point(src.get_expanded_x_min(), src.get_expanded_y_min(), src.get_expanded_z_min()),
                        ArborX::Point(src.get_expanded_x_max(), src.get_expanded_y_max(), src.get_expanded_z_max()));
    }
  }

  ArborXType data;
};

}  // namespace impl
}  // namespace stk::search

#endif
