// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
