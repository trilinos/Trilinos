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
