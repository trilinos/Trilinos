//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER

/// \file Tsqr_ApplyType.hpp
/// \brief NoTranspose, Transpose, or ConjugateTranspose
///
#ifndef __TSQR_TsqrApplyType_hpp
#define __TSQR_TsqrApplyType_hpp

#include <Tsqr_ConfigDefs.hpp>
#include <string>


namespace TSQR {

  /// \class ApplyType
  /// \brief NoTranspose, Transpose, or ConjugateTranspose
  ///
  /// An ApplyType instance represents one of three ways one can apply
  /// an implicitly stored Q factor from a QR factorization to the
  /// left of a matrix C: either as Q (NoTranspose), as Q^T
  /// (Transpose), or as Q^H (ConjugateTranspose).  Transpose and
  /// ConjugateTranspose only mean different things in complex
  /// arithmetic.  This class is a kind of "checked enum" that only
  /// allows these three values.  It knows how to go from a length-one
  /// string to the appropriate ApplyType object: "N" -> NoTranspose,
  /// "T" -> Transpose, and "C" or "H" -> ConjugateTranspose (both "C"
  /// and "H" mean the same thing).
  ///
  /// std::invalid_argument is thrown if an invalid input is given.
  class ApplyType {
  public:
    /// \brief Constructor
    ///
    /// \param op [in] One of "N", "T", "C", or "H".  Only the first
    ///   character of op is read, in a case-insensitive way.
    ApplyType (const std::string& op);

    /// \brief Whether this corresponds to (Q^T or Q^H)
    ///
    /// If op corresponds to applying Q, return false, else if op
    /// corresponds to applying Q^T or Q^H, return true.
    ///
    /// \note We lump Q^T and Q^H together because they both involve
    /// applying the Q factor pieces in the same order as they were
    /// computed in factor(), whereas applying Q involves applying
    /// those pieces in the reverse order of their computation in
    /// factor().
    bool transposed () const { return type_ != NoTranspose_; }

    //! Copy constructor
    ApplyType (const ApplyType& rhs);

    //! Assignment operator
    ApplyType& operator= (const ApplyType& rhs);

    //! Does rhs equal this?
    bool operator== (const ApplyType& rhs) const {
      return type_ == rhs.type_;
    }

    //! Does rhs not equal this?
    bool operator!= (const ApplyType& rhs) const {
      return !(type_ == rhs.type_);
    }

    //! Represents applying Q to a matrix.
    static const ApplyType NoTranspose;

    //! Represents applying Q^T (transpose of Q) to a matrix.
    static const ApplyType Transpose;

    //! Represents applying Q^H (conjugate transpose of Q) to a matrix.
    static const ApplyType ConjugateTranspose;

    /// Return a reference to the canonical LAPACK string representing
    /// the apply type.  Different for each of NoTranspose, Transpose,
    /// or ConjugateTranspose.
    ///
    /// \note This is useful for e.g., calling into LAPACK's
    ///   Householder QR routines.  this->toString().c_str() will
    ///   return a character array which LAPACK routines such as
    ///   DORMQR and ZORMQR will understand.
    const std::string& toString () const { return lapackString_; }

  private:
    /// Enum representing the underlying three possibilities for an
    /// ApplyType: Not transposed, Transposed, or Conjugate
    /// transposed.
    enum ApplyType_ { NoTranspose_, Transpose_, ConjugateTranspose_ };

    //! The state of this ApplyType.
    ApplyType_ type_;

    /// For a given ApplyType_ enum value, return the corresponding
    /// canonical LAPACK string.
    static std::string 
    enumToLapackString (const ApplyType::ApplyType_ theType);

    //! Return true if op[0] == 'T', 'C', or 'H', false otherwise.
    bool 
    decide_transposed (const std::string& op) const;

    /// Return the ApplyType_ enum value corresponding to the given
    /// string: NoTranspose_ for 'N', Transpose_ for 'T', and
    /// ConjugateTranspose_ for 'C' or 'H' (both of which mean the
    /// same thing).
    ApplyType_
    decide_apply_type (const std::string& op) const;

    /// \brief Canonical LAPACK string representing the apply type.
    ///
    /// \note We keep this around in the expectation that ApplyType
    ///   objects won't get copied much, but the underlying string
    ///   will need to be passed into LAPACK often for a single
    ///   ApplyType object.
    std::string lapackString_;
  };

} // namespace TSQR

#endif // __TSQR_TsqrApplyType_hpp
