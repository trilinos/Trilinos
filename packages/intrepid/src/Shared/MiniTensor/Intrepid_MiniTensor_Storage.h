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

#include "Intrepid_MiniTensor_Definitions.h"
#include "Teuchos_ArrayRCP.hpp"

namespace Intrepid {

/// Integer power template restricted to orders defined below
template <Index D, Index O>
struct dimension_order {
  static const Index value = 0;
};

template <Index D>
struct dimension_order<D, 1> {
  static const Index value = D;
};

template <Index D>
struct dimension_order<D, 2> {
  static const Index value = D * D;
};

template <Index D>
struct dimension_order<D, 3> {
  static const Index value = D * D * D;
};

template <Index D>
struct dimension_order<D, 4> {
  static const Index value = D * D * D * D;
};

/// Integer square root template restricted to dimensions defined below.
/// Useful for constructing a 2nd-order tensor from a 4th-order
/// tensor with static storage.
template <Index N>
struct second_to_fourth_dimension {
  static const Index value = 0;
};

template <>
struct second_to_fourth_dimension<DYNAMIC> {
  static const Index value = DYNAMIC;
};

template <>
struct second_to_fourth_dimension<1> {
  static const Index value = 1;
};

template <>
struct second_to_fourth_dimension<4> {
  static const Index value = 2;
};

template <>
struct second_to_fourth_dimension<9> {
  static const Index value = 3;
};

template <>
struct second_to_fourth_dimension<16> {
  static const Index value = 4;
};

/// Manipulation of static and dynamic dimensions.
template <Index N, Index P>
struct dimension_add {
  static const Index value = N + P;
};

template <Index P>
struct dimension_add<DYNAMIC, P> {
  static const Index value = DYNAMIC;
};

template <Index N, Index P>
struct dimension_subtract {
  static const Index value = N - P;
};

template <Index P>
struct dimension_subtract<DYNAMIC, P> {
  static const Index value = DYNAMIC;
};

///
/// Base static storage class. Simple linear access memory model.
///
template<typename T, Index N>
class Storage
{
public:
  typedef T type;

  static
  bool const
  IS_STATIC = true;

  static
  bool const
  IS_DYNAMIC = false;

  Storage() {}

  explicit
  Storage(Index const number_entries) {resize(number_entries);}

  ~Storage() {}

  T const &
  operator[](Index const i) const
  {assert(i < N); return storage_[i];}

  T &
  operator[](Index const i)
  {assert(i < N); return storage_[i];}

  Index
  size() const {return N;}

  void
  resize(Index const number_entries) {assert(number_entries == N);}

  void
  clear() {}

private:

  Storage(Storage<T, N> const & s);

  Storage<T, N> &
  operator=(Storage<T, N> const & s);

  T
  storage_[N];

};

///
/// Base dynamic storage class. Simple linear access memory model.
///
template<typename T>
class Storage<T, DYNAMIC>
{
public:
  typedef T type;

  static
  bool const
  IS_DYNAMIC = true;

  static
  bool const
  IS_STATIC = false;

  Storage();

  explicit
  Storage(Index const number_entries);

  ~Storage();

  T const &
  operator[](Index const i) const;

  T &
  operator[](Index const i);

  Index
  size() const;

  void
  resize(Index const number_entries);

  void
  clear();

private:

  Storage(Storage<T, DYNAMIC> const & s);

  Storage<T, DYNAMIC> &
  operator=(Storage<T, DYNAMIC> const & s);

  T *
  storage_;

  Index
  size_;
};

} // namespace Intrepid

#include "Intrepid_MiniTensor_Storage.i.h"

#endif // Intrepid_MiniTensor_Storage_h
