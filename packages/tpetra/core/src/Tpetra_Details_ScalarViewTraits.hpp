// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_SCALARVIEWTRAITS_HPP
#define TPETRA_DETAILS_SCALARVIEWTRAITS_HPP

///
/// \file Tpetra_Details_ScalarViewTraits.hpp
/// \brief Declaration and generic definition of traits class that
///   tells Tpetra::CrsMatrix how to pack and unpack data.
///

#include "Tpetra_Details_PackTraits.hpp"

namespace Tpetra {
namespace Details {

/// \brief Traits class for allocating a Kokkos::View<T*, D>.
///
/// \tparam T The type of the data to pack / unpack.
/// \tparam D The Kokkos "device" type; where the data live.
template<class T, class D>
struct ScalarViewTraits {
  using value_type = T;
  using device_type = D;

  /// \brief Given an instance of \c value_type allocated with the
  ///   right size, allocate and return a one-dimensional array of
  ///   \c value_type.
  ///
  /// This function lets the pack and unpack code that uses
  /// ScalarViewTraits correctly handle types that have a size
  /// specified at run time.  In particular, it's helpful if that code
  /// needs to allocate temporary buffers of \c value_type.
  /// ScalarViewTraits still assumes that all instances of \c
  /// value_type in an input or output array have the same run-time
  /// size.
  ///
  /// \param x [in] Instance of \c value_type with the correct
  ///  (run-time) size.
  ///
  /// \param numEnt [in] Number of entries in the returned
  ///   Kokkos::View.
  ///
  /// \param label [in] Optional string label of the returned
  ///   Kokkos::View.  (Kokkos::View's constructor takes a string
  ///   label, which Kokkos uses for debugging output.)
  ///
  /// \return One-dimensional array of \c value_type, all instances of
  ///   which have the same (run-time) size as \c x.
  ///
  /// \note To implementers of specializations: If the number of bytes
  ///   to pack or unpack your type may be determined at run time, you
  ///   might be able just to use this implementation as-is, and just
  ///   reimplement numValuesPerScalar().
  static Kokkos::View<value_type*, device_type>
  allocateArray (const value_type& x,
                 const size_t numEnt,
                 const std::string& label = "")
  {
    using view_type = Kokkos::View<value_type*, device_type>;

    // When the traits::specialize type is non-void this exploits
    // the fact that Kokkos::View's constructor ignores
    // size arguments beyond what the View's type specifies.  For
    // value_type = Stokhos::UQ::PCE<S>, numValuesPerScalar returns
    // something other than 1, and the constructor will actually use
    // that value.
    // Otherwise, the number of arguments must match the dynamic rank
    // (i.e. number *'s with the value_type of the View)
    const size_t numVals = PackTraits<value_type>::numValuesPerScalar (x);
    if (std::is_same<typename view_type::traits::specialize, void>::value ) {
      return view_type (label, numEnt);
    }
    else {
      return view_type (label, numEnt, numVals);
    }
  }
}; // struct ScalarViewTraits

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_SCALARVIEWTRAITS_HPP
