// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_DETAILS_MPITYPETRAITS_HPP
#define TEUCHOS_DETAILS_MPITYPETRAITS_HPP

/// \file Teuchos_Details_MpiTypeTraits.hpp
/// \brief Declaration of Teuchos::Details::MpiTypeTraits (only if
///   building with MPI)

#include "Teuchos_config.h"

#ifdef HAVE_TEUCHOS_MPI

#include <mpi.h>
#include <complex>

namespace Teuchos {
namespace Details {

/// \class MpiTypeTraits
/// \brief Traits class mapping from type T to its MPI_Datatype
/// \tparam T The type being sent or received.  T must be default
///   constructible.  It must also be either one of C++'s built-in
///   types (like \c int or \c double), or a struct or "struct-like"
///   type like <tt>std::complex<double</tt>, for which sizeof(T)
///   correctly conveys the amount of data to send or receive.
template<class T>
class MpiTypeTraits {
public:
  /// \brief Whether this class is specialized for T.
  ///
  /// If this class has <i>not</i> been specialized for T, then the
  /// return value of getType (either the one-argument or
  /// zero-argument version) is undefined.
  static const bool isSpecialized = false;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  ///
  /// It is illegal to call MPI_Type_free on a built-in MPI_Datatype.
  /// It is required to call MPI_Type_free on a non-built-in ("custom"
  /// or "derived") MPI_Datatype after use.  In the latter case, the
  /// MPI standard says that you may call MPI_Type_free on an
  /// MPI_Datatype as soon as you are done using the MPI_Datatype in
  /// your code on that process, even if there is an outstanding
  /// asynchronous operation on that process that uses the
  /// MPI_Datatype.
  ///
  /// This applies to both the one-argument and the zero-argument
  /// version of getType.  If the return value of one needs freeing,
  /// so must the return value of the other one.  (IMPLEMENTERS:
  /// Please make note of the previous sentence.)
  static const bool needsFree = false;

  /// \brief The MPI_Datatype corresponding to the given T instance.
  ///
  /// For more generality, this method requires passing in a T
  /// instance.  The method may or may not ignore this instance,
  /// depending on the type T.  The reason for passing in an instance
  /// is that some MPI_Datatype constructors, e.g., MPI_Type_struct,
  /// need actual offsets of the fields in an actual instance of T, in
  /// order to construct the MPI_Datatype safely and portably.  If T
  /// has no default constructor, we have no way of doing so without
  /// accepting a T instance.
  ///
  /// Specializations for T that do not need an instance of T in order
  /// to construct the MPI_Datatype safely, may overload this method
  /// not to require an instance of T.  However, all specializations
  /// must retain the overload that takes a T instance.  This lets
  /// users invoke this method in the same way for all types T.
  static MPI_Datatype getType (const T&) {
    // In this default implementation, isSpecialized == false, so the
    // return value of getType is undefined.  We have to return
    // something, so we return the predefined "invalid" MPI_Datatype,
    // MPI_DATATYPE_NULL.  Specializations of getType need to redefine
    // this method to return something other than MPI_DATATYPE_NULL.
    return MPI_DATATYPE_NULL;
  }
};

//
// Specializations of MpiTypeTraits.
//

/// \brief Specialization for T = char.
///
/// This requires MPI 1.2.
template<>
class MpiTypeTraits<char> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const char&) {
    return MPI_CHAR;
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_CHAR;
  }
};

/// \brief Specialization for T = unsigned char.
///
/// This requires MPI 1.2.
template<>
class MpiTypeTraits<unsigned char> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const unsigned char&) {
    return MPI_UNSIGNED_CHAR;
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_UNSIGNED_CHAR;
  }
};

#if MPI_VERSION >= 2
/// \brief Specialization for T = signed char.
///
/// This requires MPI 2.0.
template<>
class MpiTypeTraits<signed char> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const signed char&) {
    return MPI_SIGNED_CHAR;
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_SIGNED_CHAR;
  }
};
#endif // MPI_VERSION >= 2

// mfh 09 Nov 2016: amb reports on 11 Nov 2014: "I am disabling these
// specializations for now. MPI_C_DOUBLE_COMPLEX is causing a problem
// in some builds. This code was effectively turned on only yesterday
// (10 Nov 2014) when TEUCHOS_HAVE_COMPLEX was corrected to be
// HAVE_TEUCHOS_COMPLEX, so evidently there are no users of these
// specializations."
//
// mfh 14 Nov 2016: my work-around for the above issue, is
// conservatively to assume that I need MPI_VERSION >= 3 for these
// types to exist.  If I don't have MPI_VERSION >= 3, I create custom
// MPI_Datatype for these types.

namespace Impl {

/// \brief Struct for use in computeStdComplexMpiDatatype (see below).
///
/// The actual real and imaginary fields in std::complex<T> are
/// private.  While the instance methods real() and imag() return
/// references to those fields, it's not legal to take the addresses
/// of the return values of those methods ("invalid initialization"
/// errors, etc.).  Thus, we construct another struct here which
/// should have exactly the same layout, and use it in the function
/// below.
template<class T>
struct MyComplex {
  T re;
  T im;
};

/// \brief Compute MPI_Datatype for instance of std::complex<T>.
///
/// This function assumes the following:
/// <ul>
/// <li> <tt> MpiTypeTraits<T>::isSpecialized </tt> </li>
/// <li> <tt> ! MpiTypeTraits<T>::needsFree </tt>
/// <li> std::complex<T> has the same layout as
///      <tt>struct { T re; T im; };</tt> </li>
/// <li> Every instance of T has the same MPI_Datatype </li>
/// </ul>
template<class T>
MPI_Datatype
computeStdComplexMpiDatatype (const std::complex<T>& z)
{
#ifdef HAVE_TEUCHOSCORE_CXX11
  static_assert (MpiTypeTraits<T>::isSpecialized, "This function only "
                 "works if MpiTypeTraits<T>::isSpecialized.");
  static_assert (! MpiTypeTraits<T>::needsFree, "This function requires "
                 "! MpiTypeTraits<T>::needsFree, since otherwise it would "
                 "leak memory.");
#endif // HAVE_TEUCHOSCORE_CXX11

  // We assume here that every instance of T has the same
  // MPI_Datatype, i.e., has the same binary representation.
  MPI_Datatype innerDatatype = MpiTypeTraits<T>::getType (z.real ());
  MPI_Datatype outerDatatype; // return value

  // If std::complex<T> has the same layout as T[2], then we can use a
  // contiguous derived MPI_Datatype.  This is likely the only code
  // path that will execute.  Contiguous types are likely more
  // efficient for MPI to execute, and almost certainly more efficient
  // for MPI to set up.
  if (sizeof (std::complex<T>) == 2 * sizeof (T)) {
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

#ifdef HAVE_TEUCHOSCORE_CXX11
    // See documentation of MyComplex (above) for explanation.
    static_assert (sizeof (MyComplex<T>) == sizeof (std::complex<T>),
                   "Attempt to construct a struct of the same size and layout "
                   "as std::complex<T> failed.");
#endif // HAVE_TEUCHOSCORE_CXX11
    MyComplex<T> z2;

    // First entry in the struct.
    blockLengths[0] = 1;
    // Normally, &z2.re would equal &z2, but I'll be conservative and
    // actually compute the offset, even though it's probably just 0.
    //
    // Need the cast to prevent the compiler complaining about
    // subtracting addresses of different types.
    arrayOfDisplacements[0] = reinterpret_cast<char*>(&z2.re) - reinterpret_cast<char*>(&z2);
    arrayOfTypes[0] = innerDatatype;

    // Second entry in the struct.
    blockLengths[1] = 1;
    arrayOfDisplacements[1] = reinterpret_cast<char*>(&z2.im) - reinterpret_cast<char*>(&z2);
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


//! Specialization of MpiTypeTraits for std::complex<double>.
template<>
class MpiTypeTraits< std::complex<double> > {
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

  //! MPI_Datatype corresponding to the given std::complex<double> instance.
  static MPI_Datatype getType (const std::complex<double>& z) {
    if (hasMpi3) {
#if MPI_VERSION >= 3
      return MPI_C_DOUBLE_COMPLEX; // requires MPI 2.?
#else
      return MPI_DATATYPE_NULL; // FIXME (mfh 17 Nov 2016) Better to throw?
#endif // MPI_VERSION >= 3
    }
    else { // ! hasMpi3
      return Impl::computeStdComplexMpiDatatype<double> (z);
    }
  }

  //! MPI_Datatype corresponding to all std::complex<double> instances.
  static MPI_Datatype getType () {
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
      std::complex<double> z (3.0, 4.0);
      return Impl::computeStdComplexMpiDatatype<double> (z);
    }
  }
};

template<>
class MpiTypeTraits< std::complex<float> > {
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

  //! MPI_Datatype corresponding to the given std::complex<float> instance.
  static MPI_Datatype getType (const std::complex<float>& z) {
    if (hasMpi3) {
#if MPI_VERSION >= 3
      return MPI_C_FLOAT_COMPLEX; // requires MPI 2.?
#else
      return MPI_DATATYPE_NULL; // FIXME (mfh 17 Nov 2016) Better to throw?
#endif // MPI_VERSION >= 3
    }
    else { // ! hasMpi3
      return Impl::computeStdComplexMpiDatatype<float> (z);
    }
  }

  //! MPI_Datatype corresponding to all std::complex<float> instances.
  static MPI_Datatype getType () {
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
      std::complex<float> z (3.0, 4.0);
      return Impl::computeStdComplexMpiDatatype<float> (z);
    }
  }
};

//! Specialization for T = double.
template<>
class MpiTypeTraits<double> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const double&) {
    return MPI_DOUBLE;
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_DOUBLE;
  }
};

//! Specialization for T = float.
template<>
class MpiTypeTraits<float> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const float&) {
    return MPI_FLOAT;
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_FLOAT;
  }
};

//! Specialization for T = long long.
template<>
class MpiTypeTraits<long long> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const long long&) {
    return MPI_LONG_LONG; // synonym for MPI_LONG_LONG_INT in MPI 2.2
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_LONG_LONG;
  }
};


#if MPI_VERSION >= 2
/// \brief Specialization for T = long long.
///
/// This requires MPI 2.0.
template<>
class MpiTypeTraits<unsigned long long> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const unsigned long long&) {
    return MPI_UNSIGNED_LONG_LONG;
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_UNSIGNED_LONG_LONG;
  }
};
#endif // MPI_VERSION >= 2

//! Specialization for T = long.
template<>
class MpiTypeTraits<long> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const long&) {
    return MPI_LONG;
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_LONG;
  }
};

//! Specialization for T = unsigned long.
template<>
class MpiTypeTraits<unsigned long> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const unsigned long&) {
    return MPI_UNSIGNED_LONG;
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_UNSIGNED_LONG;
  }
};

//! Specialization for T = int.
template<>
class MpiTypeTraits<int> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const int&) {
    return MPI_INT;
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_INT;
  }
};

//! Specialization for T = unsigned int.
template<>
class MpiTypeTraits<unsigned int> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const unsigned int&) {
    return MPI_UNSIGNED;
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_UNSIGNED;
  }
};

//! Specialization for T = short.
template<>
class MpiTypeTraits<short> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const short&) {
    return MPI_SHORT;
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_SHORT;
  }
};

//! Specialization for T = unsigned short.
template<>
class MpiTypeTraits<unsigned short> {
public:
  //! Whether this is a defined specialization of MpiTypeTraits (it is).
  static const bool isSpecialized = true;

  /// \brief Whether you must call MPI_Type_free on the return value
  ///   of getType (both versions) after use.
  static const bool needsFree = false;

  //! MPI_Datatype corresponding to the given T instance.
  static MPI_Datatype getType (const unsigned short&) {
    return MPI_UNSIGNED_SHORT;
  }

  //! MPI_Datatype corresponding to the type T.
  static MPI_Datatype getType () {
    return MPI_UNSIGNED_SHORT;
  }
};

} // namespace Details
} // namespace Teuchos

#endif // HAVE_TEUCHOS_MPI

#endif // TEUCHOS_DETAILS_MPITYPETRAITS_HPP

