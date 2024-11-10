// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teuchos_MatrixMarket_assignScalar_hpp
#define __Teuchos_MatrixMarket_assignScalar_hpp

#include <Teuchos_as.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <string>


namespace Teuchos {
  namespace MatrixMarket {
    namespace details {

      /// \class ScalarAssigner
      /// \brief Implementation detail of assignScalar().
      ///
      /// We define a class because template functions can't (in the
      /// current C++ standard) have default template parameters.
      template<class Scalar, bool isComplex=Teuchos::ScalarTraits<Scalar>::isComplex>
      class ScalarAssigner {
      public:
        static void
        assign (Scalar& val,
                const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& real,
                const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& imag);
      };

      // Partial specialization for real Scalar types.
      template<class RealType>
      class ScalarAssigner<RealType, false> {
      public:
        static void
        assign (RealType& val,
                const typename Teuchos::ScalarTraits<RealType>::magnitudeType& real,
                const typename Teuchos::ScalarTraits<RealType>::magnitudeType& imag)
        {
          // imag had better be zero.  We're ignoring it regardless.
          (void) imag;
          val = real;
        }
      };

#ifdef HAVE_TEUCHOS_COMPLEX
      // Partial specialization for complex Scalar types.
      // MagType is the template parameter for std::complex.
      template<class MagType>
      class ScalarAssigner<std::complex<MagType>, true> {
      public:
        static void
        assign (std::complex<MagType>& val,
                const typename Teuchos::ScalarTraits<std::complex<MagType> >::magnitudeType& real,
                const typename Teuchos::ScalarTraits<std::complex<MagType> >::magnitudeType& imag)
        {
          val = std::complex<MagType> (real, imag);
        }
      };
#endif // HAVE_TEUCHOS_COMPLEX

      /// \fn assignScalar
      /// \brief Assign to a Scalar value: val = S(real, imag).
      ///
      /// We have to template it because we don't know that S is a
      /// complex type; if we write S(real,imag), the compiler will
      /// complain if S is a real type.
      template<class Scalar>
      void
      assignScalar (Scalar& val,
                    const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& real,
                    const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& imag)
      {
        ScalarAssigner<Scalar>::assign (val, real, imag);
      }

    } // namespace details
  } // namespace MatrixMarket
} // namespace Teuchos

#endif // __Teuchos_MatrixMarket_assignScalar_hpp
