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

#ifndef __Teuchos_MatrixMarket_SetScientific_hpp
#define __Teuchos_MatrixMarket_SetScientific_hpp

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <string>


namespace Teuchos {
  namespace MatrixMarket {
    namespace details {
      /// \class SetScientific
      /// \brief Temporarily make an output stream use scientific
      ///   notation with sufficient precision.
      /// \author Mark Hoemmen
      ///
      /// On construction, apply the necessary flags to the given
      /// output stream so that floating-point numbers are written in
      /// scientific notation with precision (dependent on the Scalar
      /// type) sufficient to ensure that they can be read in with the
      /// same value.  On destruction, restore the original
      /// (pre-construction) flags to the output stream.
      ///
      /// This makes SetScientific good for scope-protected alteration
      /// of the output stream's flags; no matter how the scope exits
      /// (normally or by a thrown exception), the original flags will
      /// be restored.  Hence, "temporarily" (or even "politely"): we
      /// restore the original flags on scope exit.
      ///
      /// \tparam Scalar A floating-point type, either real or
      ///   complex, for which Teuchos::ScalarTraits<Scalar> has a
      ///   specialization.  Currently we also require that
      ///   std::log10() take arguments of type Scalar, but this may
      ///   be relaxed in the future if Teuchos::ScalarTraits gets its
      ///   own log10() class method.
      template<class Scalar>
      class SetScientific {
      public:
        //! The Scalar type with which SetScientific was specialized.
        typedef Scalar scalar_type;

        /// \brief Constructor.
        ///
        /// \param out [in/out] Output stream to which to apply the
        ///   scientific notation flags.
        SetScientific (std::ostream& out) :
          out_ (out), originalFlags_ (out.flags())
        {
          typedef Teuchos::ScalarTraits<scalar_type> STS;
          typedef typename STS::magnitudeType magnitude_type;
          typedef Teuchos::ScalarTraits<magnitude_type> STM;

          // Print floating-point values in scientific notation.
          out << std::scientific;

          // We're writing decimal digits, so compute the number of
          // digits we need to get reasonable accuracy when reading
          // values back in.
          //
          // There is actually an algorithm, due to Guy Steele (yes,
          // Java's Guy Steele) et al., for idempotent printing of
          // finite-length floating-point values.  We should actually
          // implement that algorithm, but I don't have time for that
          // now.  Currently, I just print no more than (one decimal
          // digit more than (the number of decimal digits justified
          // by the precision of magnitude_type)).
          //
          // We need to use STM's log10() rather than (say) std::log10
          // here, because STM::base() returns a magnitude_type, not
          // one of C++'s standard integer types.
          const magnitude_type numDecDigits = STM::t() * STM::log10 (STM::base());

          // Round and add one.  The cast to int should not overflow
          // unless STM::t() is _extremely_ large, so we don't need to
          // check for that case here.
          const magnitude_type one = STM::one();
          const magnitude_type two = one + one;
          // Cast from magnitude_type to int, since std::ostream's
          // precision() method expects an int input.
          const int prec = 1 + Teuchos::as<int> ((two*numDecDigits + one) / two);

          // Set the number of (decimal) digits after the decimal
          // point to print.
          out.precision (prec);
        }

        /// \brief Destructor.
        ///
        /// The destructor sets the output stream's flags back to
        /// their original state, that is, the state before the
        /// constructor of this object was called.
        ~SetScientific () {
          out_.flags (originalFlags_);
        }

      private:
        //! The output stream to which to apply flags.
        std::ostream& out_;

        //! The output stream's original flags.
        std::ios_base::fmtflags originalFlags_;
      };

      // Specialization for Scalar=int.  Specializations for built-in
      // integer types let Tpetra::MatrixMarket::Reader or
      // Tpetra::MatrixMarket::Writer work with a Tpetra::CrsMatrix
      // whose Scalar template parameter is a built-in integer type.
      // These specializations are trivial because there are no
      // decimal digits to print.
      template<>
      class SetScientific<int> {
      public:
        typedef int scalar_type;
        SetScientific (std::ostream& out) {}
        ~SetScientific () {}
      };

      // Specialization for Scalar=unsigned int.
      template<>
      class SetScientific<unsigned int> {
      public:
        typedef unsigned int scalar_type;
        SetScientific (std::ostream& out) {}
        ~SetScientific () {}
      };

      // Specialization for Scalar=long.
      template<>
      class SetScientific<long> {
      public:
        typedef long scalar_type;
        SetScientific (std::ostream& out) {}
        ~SetScientific () {}
      };

      // Specialization for Scalar=unsigned long.
      template<>
      class SetScientific<unsigned long> {
      public:
        typedef unsigned long scalar_type;
        SetScientific (std::ostream& out) {}
        ~SetScientific () {}
      };

#ifdef HAVE_TEUCHOS_LONG_LONG_INT

      // Specialization for Scalar=long long.
      template<>
      class SetScientific<long long> {
      public:
        typedef long long scalar_type;
        SetScientific (std::ostream& out) {}
        ~SetScientific () {}
      };

      // Specialization for Scalar=unsigned long long.
      template<>
      class SetScientific<unsigned long long> {
      public:
        typedef unsigned long long scalar_type;
        SetScientific (std::ostream& out) {}
        ~SetScientific () {}
      };

#endif // HAVE_TEUCHOS_LONG_LONG_INT


    } // namespace details
  } // namespace MatrixMarket
} // namespace Teuchos

#endif // __Teuchos_MatrixMarket_SetScientific_hpp
