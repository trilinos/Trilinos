// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(Intrepid_MiniTensor_Storage_h)
#define Intrepid_MiniTensor_Storage_h

#include "Intrepid2_MiniTensor_Definitions.h"

namespace Intrepid2 {

/// Set to constant value if not dynamic
template<Index N, Index C>
struct dimension_const {
  static Index const value = C;
};

template<Index C>
struct dimension_const<DYNAMIC, C> {
  static Index const value = DYNAMIC;
};

/// Validate dimension
template<Index D>
struct check_static {
  static Index const
  maximum_dimension = static_cast<Index>(std::numeric_limits<Index>::digits);
#if defined(KOKKOS_HAVE_CUDA)
    // if(maximum_dimension<D) {}
    //      Kokkos::abort("Dimension is too large");}
#else
  static_assert(D < maximum_dimension, "Dimension is too large");
#endif
  static Index const value = D;
};

template<typename Store>
inline
void
check_dynamic(Index const dimension)
{
  Index const
  maximum_dimension = static_cast<Index>(std::numeric_limits<Index>::digits);

#if defined(KOKKOS_HAVE_CUDA)
    if (Store::IS_DYNAMIC == true){
      if (dimension > maximum_dimension) {Kokkos::abort("Requested dimension exceeds maximum allowed");}
   }
#else
  assert(Store::IS_DYNAMIC == true);

  if (dimension > maximum_dimension) {
    std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
    std::cerr << std::endl;
    std::cerr << "Requested dimension (" << dimension;
    std::cerr << ") exceeds maximum allowed: " << maximum_dimension;
    std::cerr << std::endl;
    exit(1);
  }
#endif
}

/// Integer power template restricted to orders defined below
template<Index D, Index O>
struct dimension_power {
  static Index const value = 0;
};

template<Index D>
struct dimension_power<D, 1> {
  static Index const value = D;
};

template<Index D>
struct dimension_power<D, 2> {
  static Index const value = D * D;
};

template<Index D>
struct dimension_power<D, 3> {
  static Index const value = D * D * D;
};

template<Index D>
struct dimension_power<D, 4> {
  static Index const value = D * D * D * D;
};

/// Integer square for manipulations between 2nd and 4rd-order tensors.
template<Index N>
struct dimension_square {
  static Index const value = 0;
};

template<>
struct dimension_square<DYNAMIC> {
  static Index const value = DYNAMIC;
};

template<>
struct dimension_square<1> {
  static Index const value = 1;
};

template<>
struct dimension_square<2> {
  static Index const value = 4;
};

template<>
struct dimension_square<3> {
  static Index const value = 9;
};

template<>
struct dimension_square<4> {
  static Index const value = 16;
};

/// Integer square root template restricted to dimensions defined below.
/// Useful for constructing a 2nd-order tensor from a 4th-order
/// tensor with static storage.
template<Index N>
struct dimension_sqrt {
  static Index const value = 0;
};

template<>
struct dimension_sqrt<DYNAMIC> {
  static Index const value = DYNAMIC;
};

template<>
struct dimension_sqrt<1> {
  static Index const value = 1;
};

template<>
struct dimension_sqrt<4> {
  static Index const value = 2;
};

template<>
struct dimension_sqrt<9> {
  static Index const value = 3;
};

template<>
struct dimension_sqrt<16> {
  static Index const value = 4;
};

/// Manipulation of static and dynamic dimensions.
template<Index N, Index P>
struct dimension_add {
  static Index const value = N + P;
};

template<Index P>
struct dimension_add<DYNAMIC, P> {
  static Index const value = DYNAMIC;
};

template<Index N, Index P>
struct dimension_subtract {
  static Index const value = N - P;
};

template<Index P>
struct dimension_subtract<DYNAMIC, P> {
  static Index const value = DYNAMIC;
};

template<Index N, Index P>
struct dimension_product {
  static Index const value = N * P;
};

template<Index N>
struct dimension_product<N, DYNAMIC> {
  static Index const value = DYNAMIC;
};

template<Index P>
struct dimension_product<DYNAMIC, P> {
  static Index const value = DYNAMIC;
};

template<>
struct dimension_product<DYNAMIC, DYNAMIC> {
  static Index const value = DYNAMIC;
};

///
/// Base static storage class. Simple linear access memory model.
///
template<typename T, Index N,typename ES>
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

  Storage()
  {
  }

  explicit
  Storage(Index const number_entries)
  {
    resize(number_entries);
  }

  ~Storage()
  {
  }

  T const &
  operator[](Index const i) const
  {
#if defined(KOKKOS_HAVE_CUDA)
   if (i>=size()) Kokkos::abort("index i in perator[] >= than size of the array");
#else
    assert(i < size());
#endif
    return storage_[i];
  }

  T &
  operator[](Index const i)
  {
#if defined(KOKKOS_HAVE_CUDA)
    if (i>=size()) Kokkos::abort("index i in perator[] >= than size of the array");
#else
    assert(i < size());
#endif
    return storage_[i];
  }

  Index
  size() const
  {
    return size_;
  }

  void
  resize(Index const number_entries)
  {
#if defined(KOKKOS_HAVE_CUDA)
     if (number_entries>N) Kokkos::abort (" IntrepidMiniTensor resize: number_entries>N");
#else
    assert(number_entries <= N);
#endif
    size_ = number_entries;
  }

  void
  clear()
  {
  }

  pointer_type
  get_pointer()
  {
    return &storage_[0];
  }

  const_pointer_type
  get_const_pointer() const
  {
    return &storage_[0];
  }

  static constexpr
  Index
  static_size()
  {
    return N;
  }

private:

  Storage(Storage<T, N, ES> const & s);

  Storage<T, N, ES> &
  operator=(Storage<T, N, ES> const & s);

  T
  storage_[N];

  Index
  size_{N};
};

///
/// Base dynamic storage class. Simple linear access memory model.
///
template<typename T, typename ES>
class Storage<T, DYNAMIC, ES>
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

  Storage()
  {
  }

  explicit
  Storage(Index const number_entries)
  {
    resize(number_entries);
  }

  ~Storage()
  {
    clear();
  }

  T const &
  operator[](Index const i) const
  {
#if defined(KOKKOS_HAVE_CUDA)
    if (i>=size()) Kokkos::abort("index i in perator[] >= than size of the array");
#else
    assert(i < size());
#endif
    return storage_[i];
  }

  T &
  operator[](Index const i)
  {
#if defined(KOKKOS_HAVE_CUDA)
    if (i>=size()) Kokkos::abort("index i in perator[] >= than size of the array");
#else
    assert(i < size());
#endif
    return storage_[i];
  }

  Index
  size() const
  {
    return size_;
  }

  void
  resize(Index const number_entries)
  {
    if (number_entries != size_) {
      clear();
      storage_ = new T[number_entries];
      size_ = number_entries;
    }
  }

  void
  clear()
  {
    if (storage_ != nullptr) {
      delete[] storage_;
      storage_ = nullptr;
      size_ = 0;
    }
  }

  pointer_type
  get_pointer()
  {
    return storage_;
  }

  const_pointer_type
  get_const_pointer() const
  {
    return storage_;
  }

  static constexpr
  Index
  static_size()
  {
    return 0;
  }

private:

  Storage(Storage<T, DYNAMIC, ES> const & s);

  Storage<T, DYNAMIC, ES> &
  operator=(Storage<T, DYNAMIC, ES> const & s);

  T *
  storage_{nullptr};

  Index
  size_{0};
};

} // namespace Intrepid2

#include "Intrepid2_MiniTensor_Storage.i.h"

#endif // Intrepid_MiniTensor_Storage_h
