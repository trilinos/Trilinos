// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tsqr_Util.hpp
/// \brief Utilities for TSQR (the Tall Skinny QR factorization)

#ifndef TSQR_UTIL_HPP
#define TSQR_UTIL_HPP

#include "Teuchos_ScalarTraits.hpp"
#include "Tsqr_MatView.hpp"

#ifdef HAVE_TPETRATSQR_COMPLEX
#  include <complex>
#endif // HAVE_TPETRATSQR_COMPLEX

#include <algorithm>
#include <ostream>

namespace TSQR {

  /// \class ScalarPrinter
  /// \brief Print a Scalar value to the given output stream
  ///
  /// \tparam Scalar The type of the value to print.
  /// \tparam isComplex Whether Scalar represents a complex number
  ///   type (such as std::complex<T>).
  ///
  /// C++ (before C++0x) doesn't let me do partial template
  /// specialization of functions.  Because of that, I can't use a
  /// template function; instead, I have to reify the function into a
  /// class ("function object").  This is typical Java style, where
  /// everything is a noun with a "run()" method; not my favorite, but
  /// it's the only way to do it.
  ///
  template<class Scalar, bool isComplex>
  class ScalarPrinter {
  public:
    ///
    /// Print elt to out
    void operator() (std::ostream& out, const Scalar& elt) const;
  };

  // Partial specialization for real Scalar
  template< class Scalar >
  class ScalarPrinter< Scalar, false > {
  public:
    void operator() (std::ostream& out, const Scalar& elt) const {
      out << elt;
    }
  };

  // Partial specialization for complex Scalar
  template< class Scalar >
  class ScalarPrinter< Scalar, true > {
  public:
    void operator() (std::ostream& out, const Scalar& elt) const {
      typedef Teuchos::ScalarTraits<Scalar> STS;
      typedef typename STS::magnitudeType magnitude_type;
      typedef Teuchos::ScalarTraits<magnitude_type> STM;

      const magnitude_type ZERO (0);
      const magnitude_type& realPart = std::real (elt);
      const magnitude_type& imagPart = std::imag (elt);

      out << realPart;
      if (imagPart < ZERO) {
        out << "-" << STM::magnitude (imagPart) << "*i";
      } else if (imagPart > ZERO) {
        out << "+" << imagPart << "*i";
      }
    }
  };

  template< class LocalOrdinal, class Scalar >
  void
  print_local_matrix (std::ostream& out,
                      const LocalOrdinal nrows_local,
                      const LocalOrdinal ncols,
                      const Scalar A[],
                      const LocalOrdinal lda)
  {
    typedef Teuchos::ScalarTraits<Scalar> STS;
    ScalarPrinter<Scalar, STS::isComplex> printer;
    for (LocalOrdinal i = 0; i < nrows_local; ++i) {
      for (LocalOrdinal j = 0; j < ncols; ++j) {
        const Scalar& curElt = A[i + j*lda];
        printer (out, curElt);
        if (j < ncols - 1) {
          out << ", ";
        }
      }
      out << ";" << std::endl;
    }
  }

  template< class Ordinal, class Scalar, class Generator >
  void
  generate_matrix (const Ordinal nrows,
                   const Ordinal ncols,
                   Scalar* const A,
                   const Ordinal lda,
                   Generator gen)
  {
    for (Ordinal j = 0; j < ncols; ++j) {
      Scalar* const A_j = &A[j*lda];
      std::generate (A_j, A_j + nrows, gen);
    }
  }

  template< class Scalar >
  class SumSquare {
  public:
    Scalar operator() (const Scalar& result, const Scalar& x) const {
      return result + x*x;
    }
  };

#ifdef HAVE_TPETRATSQR_COMPLEX
  // Specialization for complex numbers
  template<class Scalar>
  class SumSquare<std::complex<Scalar> >  {
  public:
    Scalar operator() (const std::complex<Scalar>& result,
                       const std::complex<Scalar>& x) const {
      const Scalar absval = std::norm (x);
      return result + absval * absval;
    }
  };
#endif // HAVE_TPETRATSQR_COMPLEX

} // namespace TSQR

#endif // TSQR_UTIL_HPP
