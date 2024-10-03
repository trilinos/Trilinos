// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_MpiTypeTraits.hpp"

#ifdef HAVE_TPETRACORE_MPI

namespace Teuchos {
namespace Details {
namespace Impl {

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
computeKokkosComplexMpiDatatypeImpl (const ::Kokkos::complex<T>& z)
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

    // See documentation of MyComplex (above) for explanation.
    static_assert (sizeof (MyComplex<T>) == sizeof ( ::Kokkos::complex<T>),
                   "Attempt to construct a struct of the same size and layout "
                   "as Kokkos::complex<T> failed.");
    ::Teuchos::Details::Impl::MyComplex<T> z2;

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

//! Compute MPI_Datatype for instance of Kokkos::complex<double>.
MPI_Datatype
computeKokkosComplexMpiDatatype (const ::Kokkos::complex<double>& z)
{

  return computeKokkosComplexMpiDatatypeImpl<double> (z);
}

//! Compute MPI_Datatype for instance of Kokkos::complex<float>.
MPI_Datatype
computeKokkosComplexMpiDatatype (const ::Kokkos::complex<float>& z)
{
  return computeKokkosComplexMpiDatatypeImpl<float> (z);
}

} // namespace Impl

MPI_Datatype
MpiTypeTraits< ::Kokkos::complex<double> >::
getType (const ::Kokkos::complex<double>& z)
{
  if (hasMpi3) {
#if MPI_VERSION >= 3
    return MPI_C_DOUBLE_COMPLEX; // requires MPI 2.?
#else
    return MPI_DATATYPE_NULL; // FIXME (mfh 17 Nov 2016) Better to throw?
#endif // MPI_VERSION >= 3
  }
  else { // ! hasMpi3
    return Impl::computeKokkosComplexMpiDatatype (z);
  }
}

MPI_Datatype
MpiTypeTraits< ::Kokkos::complex<double> >::
getType ()
{
  if (hasMpi3) {
#if MPI_VERSION >= 3
    return MPI_C_DOUBLE_COMPLEX; // requires MPI 2.?
#else
    return MPI_DATATYPE_NULL; // FIXME (mfh 17 Nov 2016) Better to throw?
#endif // MPI_VERSION >= 3
  }
  else { // ! hasMpi3
    // Values are arbitrary.  The function just looks at the address
    // offsets of the class fields, not their contents.
    ::Kokkos::complex<double> z (3.0, 4.0);
    return Impl::computeKokkosComplexMpiDatatype (z);
  }
}

MPI_Datatype
MpiTypeTraits< ::Kokkos::complex<float> >::
getType (const ::Kokkos::complex<float>& z)
{
  if (hasMpi3) {
#if MPI_VERSION >= 3
    return MPI_C_FLOAT_COMPLEX; // requires MPI 2.?
#else
    return MPI_DATATYPE_NULL; // FIXME (mfh 17 Nov 2016) Better to throw?
#endif // MPI_VERSION >= 3
  }
  else { // ! hasMpi3
    return Impl::computeKokkosComplexMpiDatatype (z);
  }
}

MPI_Datatype
MpiTypeTraits< ::Kokkos::complex<float> >::
getType ()
{
  if (hasMpi3) {
#if MPI_VERSION >= 3
    return MPI_C_FLOAT_COMPLEX; // requires MPI 2.?
#else
    return MPI_DATATYPE_NULL; // FIXME (mfh 17 Nov 2016) Better to throw?
#endif // MPI_VERSION >= 3
  }
  else { // ! hasMpi3
    // Values are arbitrary.  The function just looks at the address
    // offsets of the class fields, not their contents.
    ::Kokkos::complex<float> z (3.0, 4.0);
    return Impl::computeKokkosComplexMpiDatatype (z);
  }
}

} // namespace Details
} // namespace Teuchos


#endif // HAVE_TPETRACORE_MPI
