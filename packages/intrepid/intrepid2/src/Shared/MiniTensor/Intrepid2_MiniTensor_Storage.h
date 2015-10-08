// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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

#if !defined(Intrepid2_MiniTensor_Storage_h)
#define Intrepid2_MiniTensor_Storage_h

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
#if defined(HAVE_INTREPID_KOKKOSCORE) && defined(KOKKOS_HAVE_CUDA)
  static Index const
  maximum_dimension =  NPP_MAX_32U;
#else
  static Index const
  maximum_dimension = static_cast<Index>(std::numeric_limits<Index>::digits);
#endif

  static_assert(D < maximum_dimension, "Dimension is too large");
  static Index const value = D;
};

template <typename Store>
KOKKOS_INLINE_FUNCTION
void
check_dynamic(Index const dimension)
{
#if defined(HAVE_INTREPID_KOKKOSCORE) && defined(KOKKOS_HAVE_CUDA)
  static Index const
  maximum_dimension =  NPP_MAX_32U;
#else
  Index const
  maximum_dimension = static_cast<Index>(std::numeric_limits<Index>::digits);
#endif

  assert(Store::IS_DYNAMIC == true);

  if (dimension > maximum_dimension) {
#if defined(HAVE_INTREPID_KOKKOSCORE)
    Kokkos::abort("ERROR (check_dynamic) : Requested dimension exceeds maximum allowed");
#else
    std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
    std::cerr << std::endl;
    std::cerr << "Requested dimension (" << dimension;
    std::cerr << ") exceeds maximum allowed: " << maximum_dimension;
    std::cerr << std::endl;
    exit(1);
#endif
  }
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
template<typename T, Index N, class ES=NOKOKKOS>
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

  static
  bool const
  IS_KOKKOS = false;

KOKKOS_INLINE_FUNCTION
  Storage() {}

KOKKOS_INLINE_FUNCTION
  explicit
  Storage(Index const number_entries)
  {
    resize(number_entries);
  }

KOKKOS_INLINE_FUNCTION
  ~Storage() {}

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
  resize(Index const number_entries) {
#if defined (HAVE_INTREPID_KOKKOSCORE)
   if(number_entries == N) ;
   else Kokkos::abort("ERROR (MiniTensor): wrong # of entries in MiniTensor resize()");
#else
    assert(number_entries == N);
#endif
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

  static constexpr
  KOKKOS_INLINE_FUNCTION
  Index
  static_size()
  {
    return N;
  }

private:
KOKKOS_INLINE_FUNCTION
  Storage(Storage<T, N, ES> const & s);

KOKKOS_INLINE_FUNCTION
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
//Kokkos specialization of Storage
template<typename T, class ES>
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

  static 
  bool const
  IS_KOKKOS = true;

KOKKOS_INLINE_FUNCTION
  Storage()
  {
  }

KOKKOS_INLINE_FUNCTION
  explicit
  Storage(Index const number_entries)
  {
    resize(number_entries);
  }

KOKKOS_INLINE_FUNCTION
  ~Storage() {clear();}

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
  get_pointer() {return storage_;}

KOKKOS_INLINE_FUNCTION
  const_pointer_type
  get_const_pointer() const {return storage_;}

 static constexpr
  KOKKOS_INLINE_FUNCTION
  Index
  static_size()
  {
    return 0;
  }


private:
  Storage(Storage<T, DYNAMIC,ES> const & s);
KOKKOS_INLINE_FUNCTION
  Storage<T, DYNAMIC, ES> &
  operator=(Storage<T, DYNAMIC, ES> const & s);

  T *
  storage_;

  Index
  size_;
};

#if defined(HAVE_INTREPID_KOKKOSCORE) && defined(KOKKOS_HAVE_CUDA)
//CUDA specialization for the DYNAMIC storage
template<typename T>
class Storage<T, DYNAMIC, Kokkos::Cuda>
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

  static
  bool const
  IS_KOKKOS = true;

KOKKOS_INLINE_FUNCTION
  Storage() : storage_(NULL), size_(0) {}

KOKKOS_INLINE_FUNCTION
  explicit
  Storage(Index const number_entries) : storage_(NULL), size_(0)
  {resize(number_entries);}

KOKKOS_INLINE_FUNCTION
  ~Storage() {clear();}

KOKKOS_INLINE_FUNCTION
  T const &
  operator[](Index const i) const
  {assert(i < size()); return storage_[i];}

KOKKOS_INLINE_FUNCTION
  T &
  operator[](Index const i)
  {assert(i < size()); return storage_[i];}

KOKKOS_INLINE_FUNCTION
  Index
  size() const
  {return size_;}

KOKKOS_INLINE_FUNCTION
  void
  resize(Index const number_entries)
  {
    if (number_entries != size_) {
      clear(); storage_ = (T*)malloc(number_entries); size_ = number_entries;
    }
  }

KOKKOS_INLINE_FUNCTION
  void
  clear()
  {
    if (storage_ != NULL) {
      free(storage_); storage_ = NULL; size_ = 0;
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

  static constexpr
  KOKKOS_INLINE_FUNCTION
  Index
  static_size()
  {
    return 0;
  }

private:

KOKKOS_INLINE_FUNCTION
  Storage(Storage<T, DYNAMIC, Kokkos::Cuda> const & s);

KOKKOS_INLINE_FUNCTION
  Storage<T, DYNAMIC, Kokkos::Cuda> &
  operator=(Storage<T, DYNAMIC, Kokkos::Cuda> const & s);

  T *
  storage_{nullptr};

  Index
  size_{0};
};
#endif

} // namespace Intrepid

#include "Intrepid2_MiniTensor_Storage.i.h"

#endif // Intrepid2_MiniTensor_Storage_h
