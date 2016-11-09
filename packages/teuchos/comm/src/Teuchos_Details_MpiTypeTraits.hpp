// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_DETAILS_MPITYPETRAITS_HPP
#define TEUCHOS_DETAILS_MPITYPETRAITS_HPP

/// \file Teuchos_Details_MpiTypeTraits.hpp
/// \brief Declaration of Teuchos::Details::MpiTypeTraits (only if
///   building with MPI)

#include "Teuchos_config.h"

#ifdef HAVE_TEUCHOS_MPI

#include "mpi.h"
#ifdef HAVE_TEUCHOS_COMPLEX
#  include <complex>
#endif // HAVE_TEUCHOS_COMPLEX

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
  /// \brief The MPI_Datatype corresponding to the given T instance.
  ///
  /// For more generality, this method requires passing in a T
  /// instance.  The method may or may not ignore this instance,
  /// depending on the type T.
  static MPI_Datatype getType (const T&);
};

//
// Specializations of MpiTypeTraits.
//

// mfh 10 Nov 2016: Requires MPI 1.2.
template<>
class MpiTypeTraits<char> {
public:
  static MPI_Datatype getType (const char&) {
    return MPI_CHAR;
  }
};

// mfh 10 Nov 2016: Requires MPI 1.2.
template<>
class MpiTypeTraits<unsigned char> {
public:
  static MPI_Datatype getType (const unsigned char&) {
    return MPI_UNSIGNED_CHAR;
  }
};

// mfh 10 Nov 2016: Requires MPI 2.0.
#if MPI_VERSION >= 2
template<>
class MpiTypeTraits<signed char> {
public:
  static MPI_Datatype getType (const signed char&) {
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
#if 0
#ifdef HAVE_TEUCHOS_COMPLEX
template<>
class MpiTypeTraits<std::complex<double> > {
public:
  static MPI_Datatype getType (const std::complex<double>&) {
    return MPI_C_DOUBLE_COMPLEX; // requires MPI 2.?
  }
};

template<>
class MpiTypeTraits<std::complex<float> > {
public:
  static MPI_Datatype getType (const std::complex<float>&) {
    return MPI_C_FLOAT_COMPLEX; // requires MPI 2.?
  }
};
#endif // HAVE_TEUCHOS_COMPLEX
#endif // if 0

template<>
class MpiTypeTraits<double> {
public:
  static MPI_Datatype getType (const double&) {
    return MPI_DOUBLE;
  }
};

template<>
class MpiTypeTraits<float> {
public:
  static MPI_Datatype getType (const float&) {
    return MPI_FLOAT;
  }
};

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
template<>
class MpiTypeTraits<long long> {
public:
  static MPI_Datatype getType (const long long&) {
    return MPI_LONG_LONG; // synonym for MPI_LONG_LONG_INT in MPI 2.2
  }
};

// mfh 10 Nov 2016: This is new in MPI 2.0.
#if MPI_VERSION >= 2
template<>
class MpiTypeTraits<unsigned long long> {
public:
  static MPI_Datatype getType (const unsigned long long&) {
    return MPI_UNSIGNED_LONG_LONG;
  }
};
#endif // MPI_VERSION >= 2
#endif // HAVE_TEUCHOS_LONG_LONG_INT

template<>
class MpiTypeTraits<long> {
public:
  static MPI_Datatype getType (const long&) {
    return MPI_LONG;
  }
};

template<>
class MpiTypeTraits<unsigned long> {
public:
  static MPI_Datatype getType (const unsigned long&) {
    return MPI_UNSIGNED_LONG;
  }
};

template<>
class MpiTypeTraits<int> {
public:
  static MPI_Datatype getType (const int&) {
    return MPI_INT;
  }
};

template<>
class MpiTypeTraits<unsigned int> {
public:
  static MPI_Datatype getType (const unsigned int&) {
    return MPI_UNSIGNED;
  }
};

template<>
class MpiTypeTraits<short> {
public:
  static MPI_Datatype getType (const short&) {
    return MPI_SHORT;
  }
};

template<>
class MpiTypeTraits<unsigned short> {
public:
  static MPI_Datatype getType (const unsigned short&) {
    return MPI_UNSIGNED_SHORT;
  }
};

} // namespace Details
} // namespace Teuchos

#endif // HAVE_TEUCHOS_MPI

#endif // TEUCHOS_DETAILS_MPITYPETRAITS_HPP

