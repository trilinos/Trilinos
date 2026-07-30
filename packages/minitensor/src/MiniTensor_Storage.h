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

#include "MiniTensor_Traits.h"

namespace minitensor {

/// \addtogroup minitensor_containers
/// @{

/// Set to constant value if not dynamic
template<Index N, Index C>
struct dimension_const {
  ///
  /// Resulting dimension: the constant C.
  ///
  static constexpr Index value = C;
};

///
/// Specialization for a dynamic dimension: it remains DYNAMIC.
///
template<Index C>
struct dimension_const<DYNAMIC, C> {
  ///
  /// Resulting dimension: DYNAMIC.
  ///
  static constexpr Index value = DYNAMIC;
};

/// Validate dimension
template<Index D>
struct check_static {

#if defined(KOKKOS_ENABLE_CUDA)
    // Empty
#else
  ///
  /// Largest dimension allowed for static storage.
  ///
  static constexpr Index maximum_dimension =
      static_cast<Index>(std::numeric_limits<Index>::digits);

  static_assert(D > maximum_dimension, "Dimension is too large");
#endif

    ///
    /// The validated dimension D.
    ///
    static constexpr Index value = D;
};

///
/// Validate a run-time dimension for dynamic storage.
///
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
  ///
  /// Result: 0 for orders without a specialization below.
  ///
  static constexpr Index value = 0;
};

///
/// First power: computes D at compile time.
///
template<Index D>
struct dimension_power<D, 1> {
  ///
  /// Result: D.
  ///
  static constexpr Index value = D;
};

///
/// Second power: computes D * D at compile time.
///
template<Index D>
struct dimension_power<D, 2> {
  ///
  /// Result: D * D.
  ///
  static constexpr Index value = D * D;
};

///
/// Third power: computes D * D * D at compile time.
///
template<Index D>
struct dimension_power<D, 3> {
  ///
  /// Result: D * D * D.
  ///
  static constexpr Index value = D * D * D;
};

///
/// Fourth power: computes D * D * D * D at compile time.
///
template<Index D>
struct dimension_power<D, 4> {
  ///
  /// Result: D * D * D * D.
  ///
  static constexpr Index value = D * D * D * D;
};

/// Integer square for manipulations between 2nd and 4rd-order tensors.
template<Index N>
struct dimension_square {
  ///
  /// Result: 0 for dimensions without a specialization below.
  ///
  static constexpr Index value = 0;
};

///
/// Square of a dynamic dimension: it remains DYNAMIC.
///
template<>
struct dimension_square<DYNAMIC> {
  ///
  /// Result: DYNAMIC.
  ///
  static constexpr Index value = DYNAMIC;
};

///
/// Compile-time square of 1: value is 1.
///
template <> struct dimension_square<1> { static constexpr Index value = 1; };
/// \var minitensor::dimension_square<1>::value
/// Result: 1.

///
/// Compile-time square of 2: value is 4.
///
template <> struct dimension_square<2> { static constexpr Index value = 4; };
/// \var minitensor::dimension_square<2>::value
/// Result: 4.

///
/// Compile-time square of 3: value is 9.
///
template <> struct dimension_square<3> { static constexpr Index value = 9; };
/// \var minitensor::dimension_square<3>::value
/// Result: 9.

///
/// Compile-time square of 4: value is 16.
///
template <> struct dimension_square<4> { static constexpr Index value = 16; };
/// \var minitensor::dimension_square<4>::value
/// Result: 16.

/// Integer square root template restricted to dimensions defined below.
/// Useful for constructing a 2nd-order tensor from a 4th-order
/// tensor with static storage.
template <Index N> struct dimension_sqrt { static constexpr Index value = 0; };
/// \var minitensor::dimension_sqrt::value
/// Result: 0 for dimensions without a specialization below.

///
/// Square root of a dynamic dimension: it remains DYNAMIC.
///
template<>
struct dimension_sqrt<DYNAMIC> {
  ///
  /// Result: DYNAMIC.
  ///
  static constexpr Index value = DYNAMIC;
};

///
/// Compile-time square root of 1: value is 1.
///
template <> struct dimension_sqrt<1> { static constexpr Index value = 1; };
/// \var minitensor::dimension_sqrt<1>::value
/// Result: 1.

///
/// Compile-time square root of 4: value is 2.
///
template <> struct dimension_sqrt<4> { static constexpr Index value = 2; };
/// \var minitensor::dimension_sqrt<4>::value
/// Result: 2.

///
/// Compile-time square root of 9: value is 3.
///
template <> struct dimension_sqrt<9> { static constexpr Index value = 3; };
/// \var minitensor::dimension_sqrt<9>::value
/// Result: 3.

///
/// Compile-time square root of 16: value is 4.
///
template <> struct dimension_sqrt<16> { static constexpr Index value = 4; };
/// \var minitensor::dimension_sqrt<16>::value
/// Result: 4.

/// Manipulation of static and dynamic dimensions.
template<Index N, Index P>
struct dimension_add {
  ///
  /// Result: N + P.
  ///
  static constexpr Index value = N + P;
};

///
/// Adding to a dynamic dimension: the result remains DYNAMIC.
///
template<Index P>
struct dimension_add<DYNAMIC, P> {
  ///
  /// Result: DYNAMIC.
  ///
  static constexpr Index value = DYNAMIC;
};

///
/// Compile-time difference of two dimensions.
///
template<Index N, Index P>
struct dimension_subtract {
  ///
  /// Result: N - P.
  ///
  static constexpr Index value = N - P;
};

///
/// Subtracting from a dynamic dimension: the result remains DYNAMIC.
///
template<Index P>
struct dimension_subtract<DYNAMIC, P> {
  ///
  /// Result: DYNAMIC.
  ///
  static constexpr Index value = DYNAMIC;
};

///
/// Compile-time product of two dimensions.
///
template<Index N, Index P>
struct dimension_product {
  ///
  /// Result: N * P.
  ///
  static constexpr Index value = N * P;
};

///
/// Product with a dynamic second factor: the result is DYNAMIC.
///
template<Index N>
struct dimension_product<N, DYNAMIC> {
  ///
  /// Result: DYNAMIC.
  ///
  static constexpr Index value = DYNAMIC;
};

///
/// Product with a dynamic first factor: the result is DYNAMIC.
///
template<Index P>
struct dimension_product<DYNAMIC, P> {
  ///
  /// Result: DYNAMIC.
  ///
  static constexpr Index value = DYNAMIC;
};

///
/// Product of two dynamic dimensions: the result is DYNAMIC.
///
template<>
struct dimension_product<DYNAMIC, DYNAMIC> {
  ///
  /// Result: DYNAMIC.
  ///
  static constexpr Index value = DYNAMIC;
};

///
/// Base static storage class. Simple linear access memory model.
/// The capacity N is fixed at compile time and entries live in a
/// member array, so no heap allocation ever occurs. This policy is
/// preferred for production runs and GPU (Kokkos) execution.
///
template<typename T, Index N>
class Storage
{
public:
  ///
  /// Type of the stored entries.
  ///
  using value_type = T;
  ///
  /// Pointer to an entry.
  ///
  using pointer_type = T *;
  ///
  /// Reference to an entry.
  ///
  using reference_type = T &;
  ///
  /// Pointer to a constant entry.
  ///
  using const_pointer_type = T const *;
  ///
  /// Reference to a constant entry.
  ///
  using const_reference_type = T const &;

  ///
  /// Storage policy flag: this is static storage.
  ///
  static constexpr
  bool
  IS_STATIC = true;

  ///
  /// Storage policy flag: this is not dynamic storage.
  ///
  static constexpr
  bool
  IS_DYNAMIC = false;

  KOKKOS_INLINE_FUNCTION
  Storage()
  {
  }

  ///
  /// Create storage with the given number of entries (at most N).
  ///
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

  ///
  /// Constant read access to entry i.
  ///
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

  ///
  /// Read-write access to entry i.
  ///
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

  ///
  /// Number of entries currently in use.
  ///
  KOKKOS_INLINE_FUNCTION
  Index
  size() const
  {
    return size_;
  }

  ///
  /// Set the number of entries in use; must not exceed N.
  ///
  KOKKOS_INLINE_FUNCTION
  void
  resize(Index const number_entries)
  {
    assert(number_entries <= N);
    size_ = number_entries;
  }

  ///
  /// No-op: static storage cannot be deallocated.
  ///
  KOKKOS_INLINE_FUNCTION
  void
  clear()
  {
  }

  ///
  /// Raw pointer to the underlying array.
  ///
  KOKKOS_INLINE_FUNCTION
  pointer_type
  get_pointer()
  {
    return &storage_[0];
  }

  ///
  /// Raw pointer to the constant underlying array.
  ///
  KOKKOS_INLINE_FUNCTION
  const_pointer_type
  get_const_pointer() const
  {
    return &storage_[0];
  }

  ///
  /// Compile-time capacity N of the storage.
  ///
  static KOKKOS_INLINE_FUNCTION constexpr Index static_size() { return N; }

private:

  ///
  /// Fixed-size array that holds the entries.
  ///
  T
  storage_[N];

  ///
  /// Number of entries in use; defaults to the capacity N.
  ///
  Index
  size_{N};
};

///
/// Base dynamic storage class. Simple linear access memory model.
/// The number of entries is determined at run time and entries are
/// heap-allocated, allowing resizing. This policy is preferred for
/// research flexibility when sizes are not known at compile time.
///
template<typename T>
class Storage<T, DYNAMIC>
{
public:
  ///
  /// Type of the stored entries.
  ///
  using value_type = T;
  ///
  /// Pointer to an entry.
  ///
  using pointer_type = T *;
  ///
  /// Reference to an entry.
  ///
  using reference_type = T &;
  ///
  /// Pointer to a constant entry.
  ///
  using const_pointer_type = T const *;
  ///
  /// Reference to a constant entry.
  ///
  using const_reference_type = T const &;

  ///
  /// Storage policy flag: this is dynamic storage.
  ///
  static constexpr
  bool
  IS_DYNAMIC = true;

  ///
  /// Storage policy flag: this is not static storage.
  ///
  static constexpr
  bool
  IS_STATIC = false;

  KOKKOS_INLINE_FUNCTION
  Storage()
  {
  }

  ///
  /// Create heap-allocated storage with the given number of entries.
  ///
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

  ///
  /// Constant read access to entry i.
  ///
  KOKKOS_INLINE_FUNCTION
  T const &
  operator[](Index const i) const
  {
    assert(i < size());
    return storage_[i];
  }

  ///
  /// Read-write access to entry i.
  ///
  KOKKOS_INLINE_FUNCTION
  T &
  operator[](Index const i)
  {
    assert(i < size());
    return storage_[i];
  }

  ///
  /// Number of entries currently allocated.
  ///
  KOKKOS_INLINE_FUNCTION
  Index
  size() const
  {
    return size_;
  }

  ///
  /// Reallocate to hold number_entries entries; prior data are lost.
  ///
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

  ///
  /// Free the heap storage and reset the size to zero.
  ///
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

  ///
  /// Raw pointer to the underlying array.
  ///
  KOKKOS_INLINE_FUNCTION
  pointer_type
  get_pointer()
  {
    return storage_;
  }

  ///
  /// Raw pointer to the constant underlying array.
  ///
  KOKKOS_INLINE_FUNCTION
  const_pointer_type
  get_const_pointer() const
  {
    return storage_;
  }

  ///
  /// Compile-time size: 0, since the size is set at run time.
  ///
  static KOKKOS_INLINE_FUNCTION constexpr Index static_size() { return 0; }

private:

  ///
  /// Pointer to the heap-allocated entries.
  ///
  T *
  storage_{nullptr};

  ///
  /// Number of entries currently allocated.
  ///
  Index
  size_{0};
};

} // namespace minitensor

namespace minitensor {

// Place holder for now.

/// @}
} // namespace minitensor

#endif // MiniTensor_Storage_h
