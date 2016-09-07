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

#if !defined(Intrepid2_MiniTensor_TensorBase_h)
#define Intrepid2_MiniTensor_TensorBase_h

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "Intrepid2_MiniTensor_Storage.h"
#include "Intrepid2_MiniTensor_Utilities.h"

namespace Intrepid2 {

///
/// Type for setting components all at once
///
enum ComponentValue {
  ZEROS,
  ONES,
  SEQUENCE,
  RANDOM,
  RANDOM_UNIFORM,
  RANDOM_NORMAL,
  NANS
};

///
/// Base class for all vector and tensor types.
/// R, S, T: Component types
/// SR, SS, ST: Corresponding storage types.
///
template<typename T, typename ST>
class TensorBase
{
public:

  ///
  /// Component type
  ///
  using value_type = T;

  ///
  /// Storage type
  ///
  using storage_type = ST;

  ///
  /// Default constructor
  ///
  KOKKOS_INLINE_FUNCTION
  TensorBase();

  ///
  /// Constructor that initializes to NaNs
  /// \param dimension the space dimension
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  TensorBase(Index const dimension, Index const order);

  ///
  /// Create with specified value
  /// \param dimension the space dimension
  /// \param value all components are set equal to this
  ///
  KOKKOS_INLINE_FUNCTION
  TensorBase(Index const dimension, Index const order,
      ComponentValue const value);

  ///
  /// Create from a scalar
  /// \param dimension the space dimension
  /// \param s all components are set equal to this value
  ///
  KOKKOS_INLINE_FUNCTION
  TensorBase(Index const dimension, Index const order, T const & s);

  ///
  /// Create from array
  /// \param dimension the space dimension
  /// \param data_ptr pointer into the array
  ///
  // TensorBase for Kokkos Data Types (we can 't use pointers with Kokkos::View)
  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  TensorBase(
      Index const dimension,
      Index const order,
      ArrayT & data,
      iType index1);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  TensorBase(
      Index const dimension,
      Index const order,
      ArrayT & data,
      iType index1,
      iType index2);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  TensorBase(
      Index const dimension,
      Index const order,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  TensorBase(
      Index const dimension,
      Index const order,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  TensorBase(
      Index const dimension,
      Index const order,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  TensorBase(
      Index const dimension,
      Index const order,
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5,
      iType index6);

  //TensorBase for Shards and other data Types
  KOKKOS_INLINE_FUNCTION
  TensorBase(Index const dimension, Index const order, T const * data_ptr);
  ///
  /// Copy constructor
  /// \param X the values of its components are copied to the new tensor
  ///
  KOKKOS_INLINE_FUNCTION
  TensorBase(TensorBase<T, ST> const & X);

  ///
  /// Copy assignment
  /// \param X the values of its components are copied to this tensor
  ///
  KOKKOS_INLINE_FUNCTION
  TensorBase<T, ST> &
  operator=(TensorBase<T, ST> const & X);

  ///
  /// Simple destructor
  ///
  virtual
  KOKKOS_INLINE_FUNCTION
  ~TensorBase() = 0;

  ///
  /// Linear access to components
  /// \param i the index
  ///
  KOKKOS_INLINE_FUNCTION
  T const &
  operator[](Index const i) const;

  ///
  /// Linear access to components
  /// \param i the index
  ///
  KOKKOS_INLINE_FUNCTION
  T &
  operator[](Index const i);

  ///
  /// \return total number of components
  ///
  KOKKOS_INLINE_FUNCTION
  Index
  get_number_components() const;

  ///
  /// Fill components with value
  /// \param value all components are set equal to this specification
  ///
  KOKKOS_INLINE_FUNCTION
  void
  fill(ComponentValue const value);

  ///
  /// Fill components with value
  /// \param value all components are set equal to this parameter
  ///
  KOKKOS_INLINE_FUNCTION
  void
  fill(T const & s);

  ///
  /// Fill components from array defined by pointer.
  /// \param data_ptr pointer into array for filling components
  ///
  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  void fill(
      ArrayT & data,
      iType index1);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  void fill(
      ArrayT & data,
      iType index1,
      iType index2);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  void fill(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  void fill(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  void fill(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5);

  template<class ArrayT, typename iType>
  KOKKOS_INLINE_FUNCTION
  void fill(
      ArrayT & data,
      iType index1,
      iType index2,
      iType index3,
      iType index4,
      iType index5,
      iType index6);

  KOKKOS_INLINE_FUNCTION
  void fill(T const * data_ptr);

  ///
  /// Component increment
  /// \param X components are added to current components
  ///
  template<typename S, typename SS>
  KOKKOS_INLINE_FUNCTION
  TensorBase<T, ST> &
  operator+=(TensorBase<S, SS> const & X);

  ///
  /// Component decrement
  /// \param X components are subtracted to current components
  ///
  template<typename S, typename SS>
  KOKKOS_INLINE_FUNCTION
  TensorBase<T, ST> &
  operator-=(TensorBase<S, SS> const & X);

  ///
  /// Component scale
  /// \param X components are scaled by the parameter
  ///
  template<typename S>
  KOKKOS_INLINE_FUNCTION
  TensorBase<T, ST> &
  operator*=(S const & X);

  ///
  /// Component divide
  /// \param X components are divided by the parameter
  ///
  template<typename S>
  KOKKOS_INLINE_FUNCTION
  TensorBase<T, ST> &
  operator/=(S const & X);

  ///
  /// Fill with zeros
  ///
  KOKKOS_INLINE_FUNCTION
  void
  clear();

protected:

  ///
  /// \param number_components set storage for this number of components
  ///
  KOKKOS_INLINE_FUNCTION
  void
  set_number_components(Index const number_components);

  ///
  /// \return dimension
  ///
  KOKKOS_INLINE_FUNCTION
  Index
  get_dimension() const;

  ///
  /// \param dimension
  /// \param order
  ///
  KOKKOS_INLINE_FUNCTION
  void
  set_dimension(Index const dimension, Index const order);

  ///
  /// TensorBase components
  ///
  ST
  components_;

  ///
  /// dimension
  ///
  Index
  dimension_{0};

};

///
/// Frobenius norm
///
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
T
norm_f(TensorBase<T, ST> const & X);

///
/// Square of Frobenius norm
///
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
T
norm_f_square(TensorBase<T, ST> const & X);

///
/// Base addition
///
template<typename R, typename S, typename T, typename SR, typename SS,
    typename ST>
KOKKOS_INLINE_FUNCTION
void
add(
    TensorBase<R, SR> const & A,
    TensorBase<S, SS> const & B,
    TensorBase<T, ST> & C);

///
/// Base subtraction
///
template<typename R, typename S, typename T, typename SR, typename SS,
    typename ST>
KOKKOS_INLINE_FUNCTION
void
subtract(
    TensorBase<R, SR> const & A,
    TensorBase<S, SS> const & B,
    TensorBase<T, ST> & C);

///
/// Base minus
///
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
void
minus(TensorBase<T, ST> const & A, TensorBase<T, ST> & B);

///
/// Base equality
///
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
bool
equal(TensorBase<T, ST> const & A, TensorBase<T, ST> const & B);

///
/// Base not equality
///
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
bool
not_equal(TensorBase<T, ST> const & A, TensorBase<T, ST> const & B);

///
/// Base scaling
///
template<typename R, typename S, typename T, typename SR, typename ST>
KOKKOS_INLINE_FUNCTION
void
scale(TensorBase<R, SR> const & A, S const & s, TensorBase<T, ST> & B);

///
/// Base division
///
template<typename R, typename S, typename T, typename SR, typename ST>
KOKKOS_INLINE_FUNCTION
void
divide(TensorBase<R, SR> const & A, S const & s, TensorBase<T, ST> & B);

///
/// Base split (scalar divided by tensor)
///
template<typename R, typename S, typename T, typename SR, typename ST>
KOKKOS_INLINE_FUNCTION
void
split(TensorBase<R, SR> const & A, S const & s, TensorBase<T, ST> & B);

} // namespace Intrepid

#include "Intrepid2_MiniTensor_TensorBase.i.h"
#include "Intrepid2_MiniTensor_TensorBase.t.h"

#endif //Intrepid2_MiniTensor_TensorBase_h
