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

#ifndef STK_MESH_DO_OP_HPP
#define STK_MESH_DO_OP_HPP

#include "stk_mesh/base/Types.hpp"
#include <Kokkos_Core.hpp>

namespace stk {
namespace mesh {
namespace impl {

template <typename T, Operation OP>
struct DoOp;

template <typename T>
struct DoOp<T, Operation::SUM>
{
  KOKKOS_INLINE_FUNCTION
  T operator()(T lhs, T rhs) const
  { return lhs + rhs; }

  KOKKOS_INLINE_FUNCTION
  T initial_value() const
  { return T(0); }
};

template<>
struct DoOp<double,Operation::SUM>
{
  KOKKOS_INLINE_FUNCTION
  double operator()(double lhs, double rhs) const
  {
    KOKKOS_IF_ON_DEVICE((return lhs + rhs;));
    KOKKOS_IF_ON_HOST((return static_cast<long double>(lhs) + static_cast<long double>(rhs);));
  }

  KOKKOS_INLINE_FUNCTION
  double initial_value() const
  { return 0.; }
};

template <typename T>
struct DoOp<T, Operation::MIN>
{
  KOKKOS_INLINE_FUNCTION
  T operator()(T lhs, T rhs) const
  { return lhs < rhs ? lhs : rhs; }

  KOKKOS_INLINE_FUNCTION
  T initial_value() const
  { return std::numeric_limits<T>::max(); }
};

template <>
struct DoOp<std::complex<double>, Operation::MIN>
{
  std::complex<double> operator()(std::complex<double> lhs, std::complex<double> rhs) const
  {
    const auto lhsMag = std::abs(lhs);
    const auto rhsMag = std::abs(rhs);
    return (lhsMag < rhsMag) ? lhs : rhs;
  }

  std::complex<double> initial_value() const
  { return {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()}; }
};

template <>
struct DoOp<std::complex<float>, Operation::MIN>
{
  std::complex<float> operator()(std::complex<float> lhs, std::complex<float> rhs) const
  {
    const auto lhsMag = std::abs(lhs);
    const auto rhsMag = std::abs(rhs);
    return (lhsMag < rhsMag) ? lhs : rhs;
  }

  std::complex<float> initial_value() const
  { return {std::numeric_limits<float>::max(), std::numeric_limits<float>::max()}; }
};

template <>
struct DoOp<Kokkos::complex<double>, Operation::MIN>
{
  KOKKOS_INLINE_FUNCTION
  Kokkos::complex<double> operator()(Kokkos::complex<double> lhs, Kokkos::complex<double> rhs) const
  {
    const auto lhsMag = Kokkos::abs(lhs);
    const auto rhsMag = Kokkos::abs(rhs);
    return (lhsMag < rhsMag) ? lhs : rhs;
  }

  KOKKOS_INLINE_FUNCTION
  Kokkos::complex<double> initial_value() const
  { return {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()}; }
};

template <>
struct DoOp<Kokkos::complex<float>, Operation::MIN>
{
  KOKKOS_INLINE_FUNCTION
  Kokkos::complex<float> operator()(Kokkos::complex<float> lhs, Kokkos::complex<float> rhs) const
  {
    const auto lhsMag = Kokkos::abs(lhs);
    const auto rhsMag = Kokkos::abs(rhs);
    return (lhsMag < rhsMag) ? lhs : rhs;
  }

  KOKKOS_INLINE_FUNCTION
  Kokkos::complex<float> initial_value() const
  { return {std::numeric_limits<float>::max(), std::numeric_limits<float>::max()}; }
};


template <typename T>
struct DoOp<T, Operation::MAX>
{
  KOKKOS_INLINE_FUNCTION
  T operator()(T lhs, T rhs) const
  { return lhs > rhs ? lhs : rhs; }

  KOKKOS_INLINE_FUNCTION
  T initial_value() const
  { return std::numeric_limits<T>::lowest(); }
};

template <>
struct DoOp<std::complex<double>, Operation::MAX>
{
  std::complex<double> operator()(std::complex<double> lhs, std::complex<double> rhs) const
  {
    const auto lhsMag = std::abs(lhs);
    const auto rhsMag = std::abs(rhs);
    return (lhsMag > rhsMag) ? lhs : rhs;
  }

  std::complex<double> initial_value() const
  { return {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()}; }
};

template <>
struct DoOp<std::complex<float>, Operation::MAX>
{
  std::complex<float> operator()(std::complex<float> lhs, std::complex<float> rhs) const
  {
    const auto lhsMag = std::abs(lhs);
    const auto rhsMag = std::abs(rhs);
    return (lhsMag > rhsMag) ? lhs : rhs;
  }

  std::complex<float> initial_value() const
  { return {std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()}; }
};

template <>
struct DoOp<Kokkos::complex<double>, Operation::MAX>
{
  KOKKOS_INLINE_FUNCTION
  Kokkos::complex<double> operator()(Kokkos::complex<double> lhs, Kokkos::complex<double> rhs) const
  {
    const auto lhsMag = Kokkos::abs(lhs);
    const auto rhsMag = Kokkos::abs(rhs);
    return (lhsMag > rhsMag) ? lhs : rhs;
  }

  KOKKOS_INLINE_FUNCTION
  Kokkos::complex<double> initial_value() const
  { return {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()}; }
};

template <>
struct DoOp<Kokkos::complex<float>, Operation::MAX>
{
  KOKKOS_INLINE_FUNCTION
  Kokkos::complex<float> operator()(Kokkos::complex<float> lhs, Kokkos::complex<float> rhs) const
  {
    const auto lhsMag = Kokkos::abs(lhs);
    const auto rhsMag = Kokkos::abs(rhs);
    return (lhsMag > rhsMag) ? lhs : rhs;
  }

  KOKKOS_INLINE_FUNCTION
  Kokkos::complex<float> initial_value() const
  { return {std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()}; }
};

}
}
}

#endif
