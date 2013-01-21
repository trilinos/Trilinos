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

#if !defined(Intrepid_MiniTensor_Utilities_h)
#define Intrepid_MiniTensor_Utilities_h

#include "Sacado.hpp"

namespace Intrepid {

  ///
  /// Sign function
  ///
  template<typename T>
  int
  sgn(T const & s);

  ///
  /// Copysign function
  ///
  template<typename T>
  T
  copysign(T const & a, T const & b);

  ///
  /// NaN function. Necessary to choose the proper underlying NaN
  /// for non-floating-point types.
  /// Assumption: non-floating-point types have a typedef that
  /// determines the underlying floating-point type.
  ///
  template<typename T>
  typename Sacado::ScalarType<T>::type
  not_a_number();

  ///
  /// Machine epsilon function. Necessary to choose the proper underlying
  /// machine epsilon for non-floating-point types.
  /// Assumption: non-floating-point types have a typedef that
  /// determines the underlying floating-point type.
  ///
  template<typename T>
  typename Sacado::ScalarType<T>::type
  machine_epsilon();

  ///
  /// Duet type. Holder of two objects of the same type.
  /// Useful as return type for functions that need to return two objects.
  ///
  template<typename T>
  struct Duet
  {
    T first;
    T second;
  };

  ///
  /// Triplet type. Holder of three objects of the same type.
  /// Useful as return type for functions that need to return three objects.
  ///
  template<typename T>
  struct Triplet
  {
    T first;
    T second;
    T third;
  };

  ///
  /// Create a Duet structure.
  ///
  template<typename T>
  inline
  Duet<T>
  make_duet(T const & a, T const & b);

  ///
  /// Create a Triplet structure.
  ///
  template<typename T>
  inline
  Triplet<T>
  make_triplet(T const & a, T const & b, T const & c);

  ///
  /// Tie function template to hold values of functions
  /// that return a Duet.
  ///
  template<typename T>
  inline
  Duet<T>
  tie(T & a, T& b);

  ///
  /// Tie function template to hold values of functions
  /// that return a Duet.
  ///
  template<typename T>
  inline
  Triplet<T>
  tie(T & a, T & b, T & c);

} // namespace Intrepid

#include "Intrepid_MiniTensor_Utilities.i.cc"

#endif // Intrepid_MiniTensor_Utilities_h
