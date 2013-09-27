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

#if !defined(Intrepid_MiniTensor_TensorBase_h)
#define Intrepid_MiniTensor_TensorBase_h

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "Intrepid_MiniTensor_Storage.h"
#include "Intrepid_MiniTensor_Utilities.h"

namespace Intrepid {

///
/// Type for setting components all at once
///
enum ComponentValue {ZEROS, ONES, SEQUENCE, RANDOM, NANS};

///
/// Base class for all vector and tensor types.
///
template<typename T, typename Store>
class TensorBase
{
public:

  ///
  /// Component type
  ///
  typedef T type;

  ///
  /// Default constructor
  ///
  TensorBase();

  ///
  /// Constructor that initializes to NaNs
  /// \param dimension the space dimension
  ///
  explicit
  TensorBase(Index const dimension, Index const order);

  ///
  /// Create with specified value
  /// \param dimension the space dimension
  /// \param value all components are set equal to this
  ///
  TensorBase(Index const dimension, Index const order,
      ComponentValue const value);

  ///
  /// Create from a scalar
  /// \param dimension the space dimension
  /// \param s all components are set equal to this value
  ///
  TensorBase(Index const dimension, Index const order, T const & s);

  ///
  /// Create from array
  /// \param dimension the space dimension
  /// \param data_ptr pointer into the array
  ///
  TensorBase(Index const dimension, Index const order, T const * data_ptr);

  ///
  /// Copy constructor
  /// \param X the values of its components are copied to the new tensor
  ///
  TensorBase(TensorBase<T, Store> const & X);

  ///
  /// Copy assignment
  /// \param X the values of its components are copied to this tensor
  ///
  TensorBase<T, Store> &
  operator=(TensorBase<T, Store> const & X);

  ///
  /// Simple destructor
  ///
  virtual
  ~TensorBase() = 0;

  ///
  /// Linear access to components
  /// \param i the index
  ///
  T const &
  operator[](Index const i) const;

  ///
  /// Linear access to components
  /// \param i the index
  ///
  T &
  operator[](Index const i);

  ///
  /// \return total number of components
  ///
  Index
  get_number_components() const;

  ///
  /// Fill components with value
  /// \param value all components are set equal to this specification
  ///
  void
  fill(ComponentValue const value);

  ///
  /// Fill components with value
  /// \param value all components are set equal to this parameter
  ///
  void
  fill(T const & s);

  ///
  /// Fill components from array defined by pointer.
  /// \param data_ptr pointer into array for filling components
  ///
  void
  fill(T const * data_ptr);

  ///
  /// Component increment
  /// \param X components are added to current components
  ///
  TensorBase<T, Store> &
  operator+=(TensorBase<T, Store> const & X);

  ///
  /// Component decrement
  /// \param X components are subtracted to current components
  ///
  TensorBase<T, Store> &
  operator-=(TensorBase<T, Store> const & X);

  ///
  /// Fill with zeros
  ///
  void
  clear();

protected:

  ///
  /// \param number_components set storage for this number of components
  ///
  void
  set_number_components(Index const number_components);

  ///
  /// \return dimension
  ///
  Index
  get_dimension() const;

  ///
  /// \param dimension
  /// \param order
  ///
  void
  set_dimension(Index const dimension, Index const order);

  ///
  /// TensorBase components
  ///
  Store
  components_;

  ///
  /// dimension
  ///
  Index
  dimension_;

};

///
/// Frobenius norm
///
template<typename T, typename Store>
T
norm_f(TensorBase<T, Store> const & X);

///
/// Square of Frobenius norm
///
template<typename T, typename Store>
T
norm_f_square(TensorBase<T, Store> const & X);

///
/// Base addition
///
template<typename R, typename S, typename T, typename Store>
void
add(
    TensorBase<R, Store> const & A,
    TensorBase<S, Store> const & B,
    TensorBase<T, Store> & C);

///
/// Base subtraction
///
template<typename R, typename S, typename T, typename Store>
void
subtract(
    TensorBase<R, Store> const & A,
    TensorBase<S, Store> const & B,
    TensorBase<T, Store> & C);

///
/// Base minus
///
template<typename T, typename Store>
void
minus(TensorBase<T, Store> const & A, TensorBase<T, Store> & B);

///
/// Base equality
///
template<typename T, typename Store>
bool
equal(TensorBase<T, Store> const & A, TensorBase<T, Store> const & B);

///
/// Base not equality
///
template<typename T, typename Store>
bool
not_equal(TensorBase<T, Store> const & A, TensorBase<T, Store> const & B);

///
/// Base scaling
///
template<typename R, typename S, typename T, typename Store>
void
scale(TensorBase<R, Store> const & A, S const & s, TensorBase<T, Store> & B);

///
/// Base division
///
template<typename R, typename S, typename T, typename Store>
void
divide(TensorBase<R, Store> const & A, S const & s, TensorBase<T, Store> & B);

///
/// Base split (scalar divided by tensor)
///
template<typename R, typename S, typename T, typename Store>
void
split(TensorBase<R, Store> const & A, S const & s, TensorBase<T, Store> & B);

} // namespace Intrepid

#include "Intrepid_MiniTensor_TensorBase.i.h"
#include "Intrepid_MiniTensor_TensorBase.t.h"

#endif //Intrepid_MiniTensor_TensorBase_h
