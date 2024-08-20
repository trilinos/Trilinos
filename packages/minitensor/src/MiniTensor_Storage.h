// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Kokkos_Macros.hpp"
#if !defined(MiniTensor_Storage_h)
#define MiniTensor_Storage_h

#include "MiniTensor_Definitions.h"

namespace minitensor {

/// Set to constant value if not dynamic
template<Index N, Index C>
struct dimension_const {
  static constexpr Index value = C;
};

template<Index C>
struct dimension_const<DYNAMIC, C> {
  static constexpr Index value = DYNAMIC;
};

/// Validate dimension
template<Index D>
struct check_static {

#if defined(KOKKOS_ENABLE_CUDA)
    // Empty
#else
  static constexpr Index maximum_dimension =
      static_cast<Index>(std::numeric_limits<Index>::digits);

  static_assert(D > maximum_dimension, "Dimension is too large");
#endif

    static constexpr Index value = D;
};

template<typename Store>
inline
void
check_dynamic(Index const dimension)
{
  Index const
  maximum_dimension = static_cast<Index>(std::numeric_limits<Index>::digits);

  assert(Store::IS_DYNAMIC == true);

  if (dimension > maximum_dimension) {
    MT_ERROR_EXIT("Requested dimension exceeds maximum allowed: %d", dimension);
  }
}

/// Integer power template restricted to orders defined below
template<Index D, Index O>
struct dimension_power {
  static constexpr Index value = 0;
};

template<Index D>
struct dimension_power<D, 1> {
  static constexpr Index value = D;
};

template<Index D>
struct dimension_power<D, 2> {
  static constexpr Index value = D * D;
};

template<Index D>
struct dimension_power<D, 3> {
  static constexpr Index value = D * D * D;
};

template<Index D>
struct dimension_power<D, 4> {
  static constexpr Index value = D * D * D * D;
};

/// Integer square for manipulations between 2nd and 4rd-order tensors.
template<Index N>
struct dimension_square {
  static constexpr Index value = 0;
};

template<>
struct dimension_square<DYNAMIC> {
  static constexpr Index value = DYNAMIC;
};

template <> struct dimension_square<1> { static constexpr Index value = 1; };

template <> struct dimension_square<2> { static constexpr Index value = 4; };

template <> struct dimension_square<3> { static constexpr Index value = 9; };

template <> struct dimension_square<4> { static constexpr Index value = 16; };

/// Integer square root template restricted to dimensions defined below.
/// Useful for constructing a 2nd-order tensor from a 4th-order
/// tensor with static storage.
template <Index N> struct dimension_sqrt { static constexpr Index value = 0; };

template<>
struct dimension_sqrt<DYNAMIC> {
  static constexpr Index value = DYNAMIC;
};

template <> struct dimension_sqrt<1> { static constexpr Index value = 1; };

template <> struct dimension_sqrt<4> { static constexpr Index value = 2; };

template <> struct dimension_sqrt<9> { static constexpr Index value = 3; };

template <> struct dimension_sqrt<16> { static constexpr Index value = 4; };

/// Manipulation of static and dynamic dimensions.
template<Index N, Index P>
struct dimension_add {
  static constexpr Index value = N + P;
};

template<Index P>
struct dimension_add<DYNAMIC, P> {
  static constexpr Index value = DYNAMIC;
};

template<Index N, Index P>
struct dimension_subtract {
  static constexpr Index value = N - P;
};

template<Index P>
struct dimension_subtract<DYNAMIC, P> {
  static constexpr Index value = DYNAMIC;
};

template<Index N, Index P>
struct dimension_product {
  static constexpr Index value = N * P;
};

template<Index N>
struct dimension_product<N, DYNAMIC> {
  static constexpr Index value = DYNAMIC;
};

template<Index P>
struct dimension_product<DYNAMIC, P> {
  static constexpr Index value = DYNAMIC;
};

template<>
struct dimension_product<DYNAMIC, DYNAMIC> {
  static constexpr Index value = DYNAMIC;
};

///
/// Base static storage class. Simple linear access memory model.
///
template<typename T, Index N>
class Storage
{
public:
  using value_type = T;
  using pointer_type = T *;
  using reference_type = T &;
  using const_pointer_type = T const *;
  using const_reference_type = T const &;

  static constexpr
  bool
  IS_STATIC = true;

  static constexpr
  bool
  IS_DYNAMIC = false;

  KOKKOS_INLINE_FUNCTION
  Storage()
  {
  }

  explicit KOKKOS_INLINE_FUNCTION Storage(Index const number_entries) {
    resize(number_entries);
  }

  Storage(Storage<T, N> const & s) = delete;

  Storage<T, N> &
  operator=(Storage<T, N> const & s) = delete;

  KOKKOS_INLINE_FUNCTION
  ~Storage()
  {
  }

  KOKKOS_INLINE_FUNCTION
  T const &
  operator[](Index const i) const
  {
    assert(i < size());
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
    return storage_[i];
#pragma GCC diagnostic pop
  }

  KOKKOS_INLINE_FUNCTION
  T &
  operator[](Index const i)
  {
    assert(i < size());
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
    return storage_[i];
#pragma GCC diagnostic pop
  }

  KOKKOS_INLINE_FUNCTION
  Index
  size() const
  {
    return size_;
  }

  KOKKOS_INLINE_FUNCTION
  void
  resize(Index const number_entries)
  {
    assert(number_entries <= N);
    size_ = number_entries;
  }

  KOKKOS_INLINE_FUNCTION
  void
  clear()
  {
  }

  KOKKOS_INLINE_FUNCTION
  pointer_type
  get_pointer()
  {
    return &storage_[0];
  }

  KOKKOS_INLINE_FUNCTION
  const_pointer_type
  get_const_pointer() const
  {
    return &storage_[0];
  }

  static KOKKOS_INLINE_FUNCTION constexpr Index static_size() { return N; }

private:

  T
  storage_[N];

  Index
  size_{N};
};

///
/// Base dynamic storage class. Simple linear access memory model.
///
template<typename T>
class Storage<T, DYNAMIC>
{
public:
  using value_type = T;
  using pointer_type = T *;
  using reference_type = T &;
  using const_pointer_type = T const *;
  using const_reference_type = T const &;

  static constexpr
  bool
  IS_DYNAMIC = true;

  static constexpr
  bool
  IS_STATIC = false;

  KOKKOS_INLINE_FUNCTION
  Storage()
  {
  }

  explicit KOKKOS_INLINE_FUNCTION Storage(Index const number_entries) {
    resize(number_entries);
  }

  Storage(Storage<T, DYNAMIC> const & s) = delete;

  Storage<T, DYNAMIC> &
  operator=(Storage<T, DYNAMIC> const & s) = delete;

  KOKKOS_INLINE_FUNCTION
  ~Storage()
  {
    clear();
  }

  KOKKOS_INLINE_FUNCTION
  T const &
  operator[](Index const i) const
  {
    assert(i < size());
    return storage_[i];
  }

  KOKKOS_INLINE_FUNCTION
  T &
  operator[](Index const i)
  {
    assert(i < size());
    return storage_[i];
  }

  KOKKOS_INLINE_FUNCTION
  Index
  size() const
  {
    return size_;
  }

  KOKKOS_INLINE_FUNCTION
  void
  resize(Index const number_entries)
  {
    if (number_entries != size_) {
      clear();
      storage_ = new T[number_entries];
      size_ = number_entries;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void
  clear()
  {
    if (storage_ != nullptr) {
      delete[] storage_;
      storage_ = nullptr;
      size_ = 0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  pointer_type
  get_pointer()
  {
    return storage_;
  }

  KOKKOS_INLINE_FUNCTION
  const_pointer_type
  get_const_pointer() const
  {
    return storage_;
  }

  static KOKKOS_INLINE_FUNCTION constexpr Index static_size() { return 0; }

private:

  T *
  storage_{nullptr};

  Index
  size_{0};
};

} // namespace minitensor

#include "MiniTensor_Storage.i.h"

#endif // MiniTensor_Storage_h
