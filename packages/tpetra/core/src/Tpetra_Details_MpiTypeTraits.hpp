// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_MPITYPETRAITS_HPP
#define TPETRA_DETAILS_MPITYPETRAITS_HPP

/// \file Tpetra_Details_MpiTypeTraits.hpp
/// \brief Add specializations of Teuchos::Details::MpiTypeTraits for
///   Kokkos::complex<float> and Kokkos::complex<double>, and import
///   Teuchos::Details::MpiTypeTraits into the Tpetra::Details
///   namespace.
///
/// \warning This file, and its contents are implementation details of
///   Tpetra.  DO NOT DEPEND ON THEM!
/// \warning This file's contents only exist if building with MPI.
///
/// This file exists mainly to let Tpetra developers add
/// specializations to Teuchos::Details::MpiTypeTraits.  Tpetra in
/// particular needs specializations for the Kokkos versions of
/// std::complex, since std::complex does not work with Kokkos::View.
///
/// If you want to add a new Scalar type to Tpetra, and Kokkos already
/// supports that type, you should add a specialization of
/// MpiTypeTraits for your Scalar types.  If you want to add a new
/// Scalar type to Tpetra, and Kokkos does not support that type,
/// first begin by adding a specialization of
/// Kokkos::ArithTraits for your type, that maps it to a type
/// that Kokkos does support.  Then, add a specialization of
/// MpiTypeTraits for your Scalar types.

#include "TpetraCore_config.h"
#ifdef HAVE_TPETRACORE_MPI

#include "Kokkos_Complex.hpp"
#include "Teuchos_Details_MpiTypeTraits.hpp"
// Include this file for uintptr_t.
// Windows (Visual Studio 2015) claims to have that type too:
//
// https://msdn.microsoft.com/en-us/library/323b6b3k.aspx
#include <cstdint>

namespace Teuchos {
namespace Details {

//! Specialization of MpiTypeTraits for Kokkos::complex<double>.
template<>
class MpiTypeTraits< ::Kokkos::complex<double> > {
private:
#if MPI_VERSION >= 3
  static const bool hasMpi3 = true;
#else
  static const bool hasMpi3 = false;
#endif // MPI_VERSION >= 3

public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = ! hasMpi3;

  //! MPI_Datatype corresponding to the given Kokkos::complex<double> instance.
  static MPI_Datatype getType (const ::Kokkos::complex<double>& z);

  //! MPI_Datatype corresponding to all Kokkos::complex<double> instances.
  static MPI_Datatype getType ();
};

template<>
class MpiTypeTraits< ::Kokkos::complex<float> > {
private:
#if MPI_VERSION >= 3
  static const bool hasMpi3 = true;
#else
  static const bool hasMpi3 = false;
#endif // MPI_VERSION >= 3

public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = ! hasMpi3;

  //! MPI_Datatype corresponding to the given Kokkos::complex<float> instance.
  static MPI_Datatype getType (const ::Kokkos::complex<float>& z);

  //! MPI_Datatype corresponding to all Kokkos::complex<float> instances.
  static MPI_Datatype getType ();
};

} // namespace Details
} // namespace Teuchos

namespace Tpetra {
namespace Details {

// Import MpiTypeTraits into the Tpetra::Details namespace.
using ::Teuchos::Details::MpiTypeTraits;

} // namespace Details
} // namespace Tpetra

#endif // HAVE_TPETRACORE_MPI
#endif // TPETRA_DETAILS_MPITYPETRAITS_HPP
