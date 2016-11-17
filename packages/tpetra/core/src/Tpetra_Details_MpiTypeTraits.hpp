/*
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

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
/// Kokkos::Details::ArithTraits for your type, that maps it to a type
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

namespace Impl {

/// \brief Struct for use in computeKokkosComplexMpiDatatype (see below).
///
/// The actual re_ and im_ fields in Kokkos::complex<T> are private.
/// While the instance methods real() and imag() return references to
/// those fields, it's not legal to take the addresses of the return
/// values of those methods ("invalid initialization" errors, etc.).
/// Thus, we construct another struct here which should have exactly
/// the same layout, and use it in the function below.
template<class T>
struct MyComplex {
  T re;
  T im;
};

/// \brief Compute MPI_Datatype for instance of Kokkos::complex<T>.
///
/// This function assumes the following:
/// <ul>
/// <li> <tt> MpiTypeTraits<T>::isSpecialized </tt> </li>
/// <li> <tt> ! MpiTypeTraits<T>::needsFree </tt>
/// <li> Kokkos::complex<T> has the same layout as
///      <tt>struct { T re; T im; };</tt> </li>
/// <li> Every instance of T has the same MPI_Datatype </li>
/// </ul>
template<class T>
MPI_Datatype
computeKokkosComplexMpiDatatype (const ::Kokkos::complex<T>& z)
{
  static_assert (MpiTypeTraits<T>::isSpecialized, "This function only "
                 "works if MpiTypeTraits<T>::isSpecialized.");
  static_assert (! MpiTypeTraits<T>::needsFree, "This function requires "
                 "! MpiTypeTraits<T>::needsFree, since otherwise it would "
                 "leak memory.");
  // We assume here that every instance of T has the same
  // MPI_Datatype, i.e., has the same binary representation.
  MPI_Datatype innerDatatype = MpiTypeTraits<T>::getType (z.real ());
  MPI_Datatype outerDatatype; // return value

  // If Kokkos::complex<T> has the same layout as T[2], then we can
  // use a contiguous derived MPI_Datatype.  This is likely the only
  // code path that will execute.  Contiguous types are likely more
  // efficient for MPI to execute, and almost certainly more efficient
  // for MPI to set up.
  if (sizeof ( ::Kokkos::complex<T>) == 2 * sizeof (T)) {
    (void) MPI_Type_contiguous (2, innerDatatype, &outerDatatype);
  }
  else { // must use the general struct approach
    // I borrowed and adapted the code below from the MPICH
    // documentation:
    //
    // www.mpich.org/static/docs/v3.1/www3/MPI_Type_struct.html
    int blockLengths[3];
    MPI_Aint arrayOfDisplacements[3];
    MPI_Datatype arrayOfTypes[3];
    MPI_Datatype outerDatatype;

    // See documentation of MyComplex (above) for explanation.
    static_assert (sizeof (MyComplex<T>) == sizeof ( ::Kokkos::complex<T>),
                   "Attempt to construct a struct of the same size and layout "
                   "as Kokkos::complex<T> failed.");
    MyComplex<T> z2;

    // First entry in the struct.
    blockLengths[0] = 1;
    // Normally, &z2.re would equal &z2, but I'll be conservative and
    // actually compute the offset, even though it's probably just 0.
    //
    // Need the cast to prevent the compiler complaining about
    // subtracting addresses of different types.
    arrayOfDisplacements[0] = reinterpret_cast<uintptr_t> (&z2.re) - reinterpret_cast<uintptr_t> (&z2);
    arrayOfTypes[0] = innerDatatype;

    // Second entry in the struct.
    blockLengths[1] = 1;
    arrayOfDisplacements[1] = reinterpret_cast<uintptr_t> (&z2.im) - reinterpret_cast<uintptr_t> (&z2);
    arrayOfTypes[1] = innerDatatype;

#if MPI_VERSION < 2
    // Upper bound of the struct.
    blockLengths[2] = 1;
    arrayOfDisplacements[2] = sizeof (MyComplex<T>);
    arrayOfTypes[2] = MPI_UB; // "upper bound type"; signals end of struct
#endif // MPI_VERSION < 2

    // Define the MPI_Datatype.
#if MPI_VERSION < 2
    (void) MPI_Type_struct (3, blockLengths, arrayOfDisplacements,
                            arrayOfTypes, &outerDatatype);
#else
    // Don't include the upper bound with MPI_Type_create_struct.
    (void) MPI_Type_create_struct (2, blockLengths, arrayOfDisplacements,
                                   arrayOfTypes, &outerDatatype);
#endif // MPI_VERSION < 2
  }

  MPI_Type_commit (&outerDatatype);
  return outerDatatype;
}

} // namespace Impl

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
  static MPI_Datatype getType (const ::Kokkos::complex<double>& z) {
    if (hasMpi3) {
      return MPI_C_DOUBLE_COMPLEX; // requires MPI 2.?
    }
    else { // ! hasMpi3
      return Impl::computeKokkosComplexMpiDatatype<double> (z);
    }
  }

  //! MPI_Datatype corresponding to all Kokkos::complex<double> instances.
  static MPI_Datatype getType () {
    if (hasMpi3) {
      return MPI_C_DOUBLE_COMPLEX; // requires MPI 2.?
    }
    else { // ! hasMpi3
      // Values are arbitrary.  The function just looks at the address
      // offsets of the class fields, not their contents.
      ::Kokkos::complex<double> z (3.0, 4.0);
      return Impl::computeKokkosComplexMpiDatatype<double> (z);
    }
  }
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
  static MPI_Datatype getType (const ::Kokkos::complex<float>& z) {
    if (hasMpi3) {
      return MPI_C_FLOAT_COMPLEX; // requires MPI 2.?
    }
    else { // ! hasMpi3
      return Impl::computeKokkosComplexMpiDatatype<float> (z);
    }
  }

  //! MPI_Datatype corresponding to all Kokkos::complex<float> instances.
  static MPI_Datatype getType () {
    if (hasMpi3) {
      return MPI_C_FLOAT_COMPLEX; // requires MPI 2.?
    }
    else { // ! hasMpi3
      // Values are arbitrary.  The function just looks at the address
      // offsets of the class fields, not their contents.
      ::Kokkos::complex<float> z (3.0, 4.0);
      return Impl::computeKokkosComplexMpiDatatype<float> (z);
    }
  }
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
