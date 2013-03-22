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

#ifndef __Teuchos_MatrixMarket_SymmetrizingAdder_hpp
#define __Teuchos_MatrixMarket_SymmetrizingAdder_hpp

#include <Teuchos_as.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <string>


namespace Teuchos {
  namespace MatrixMarket {
    // Anonymous namespace for helper functions for SymmetrizingAdder.
    namespace {
      bool isSkew (const std::string& symmType) {
        return symmType.size() >= 4 && symmType.substr(0,4) == "skew";
      }

      bool isConj (const std::string& symmType) {
        return std::string::npos != symmType.find ("hermitian");
      }

      bool needsSymmetrization (const std::string& symmType) {
        return symmType != "general";
      }
    } // namespace (anonymous)

    /// \class SymmetrizingAdder
    /// \author Mark Hoemmen
    /// \brief Adds entries with optional symmetry to a sparse matrix.
    ///
    /// This class wraps any existing class (AdderType) with the
    /// interface shown below.  Given the Matrix Market symmetry type,
    /// this class' corresponding operator() may invoke AdderType's
    /// operator() twice, in order to add entry (j,i) if entry (i,j)
    /// is to be added.
    ///
    /// \tparam AdderType A class with at least the following interface:
    ///   \code
    ///   class AdderType {
    ///   public:
    ///     typedef ... index_type; // Ellipsis represents the actual type
    ///     typedef ... value_type; // Ellipsis represents the actual type
    ///     void operator() (const index_type, const index_type, const value_type&);
    ///   };
    ///   \endcode
    template<class AdderType>
    class SymmetrizingAdder {
    public:
      //! The type of indices of the sparse matrix
      typedef typename AdderType::index_type index_type;
      //! The type of entries of the sparse matrix
      typedef typename AdderType::value_type value_type;

      /// \brief Constructor
      ///
      /// \param adder [in/out] The wrapped AdderType instance
      ///
      /// \param symmType [in] Canonical Matrix Market string
      ///   representing the symmetry storage type of the matrix data.
      SymmetrizingAdder (const Teuchos::RCP<AdderType>& adder,
                         const std::string& symmType) :
        adder_ (adder),
        symmetrize_ (needsSymmetrization (symmType)),
        conjugate_ (isConj (symmType)),
        skew_ (isSkew (symmType))
      {}

      //! Add value A_ij to entry (i,j), and optionally symmetrize.
      void
      operator() (const index_type i,
                  const index_type j,
                  const value_type& Aij)
      {
        AdderType& theAdder = *adder_;

        theAdder (i, j, Aij);
        if (symmetrize_ && i != j) {
          typedef Teuchos::ScalarTraits<value_type> STS;
          const value_type Aji = skew_ ?
            -(conjugate_ ? STS::conjugate(Aij) : Aij) :
            (conjugate_ ? STS::conjugate(Aij) : Aij);
          // The optional fourth argument (which defaults to true)
          // specifies whether or not to count the entry against the
          // total expected number of entries.  We don't want to count
          // this entry because it wasn't part of the original data;
          // we inserted it because the caller doesn't want symmetric
          // storage.  The original data's total expected number of
          // entries only counts the entries that are in the original
          // data, not those that we insert.
          theAdder (j, i, Aji, false);
        }
      }

      /// \brief Persisting non-const view of the underlying adder object.
      ///
      /// This violates encapsulation, so please be careful with this.
      Teuchos::RCP<AdderType> getAdder() const {
        return adder_;
      }

    private:
      //! The wrapped AdderType instance.
      Teuchos::RCP<AdderType> adder_;
      //! Whether to do symmetrization at all.
      bool symmetrize_;
      //! Whether to conjugate when symmetrizing.
      bool conjugate_;
      //! Whether to negate when symmetrizing.
      bool skew_;
    };

  } // namespace MatrixMarket
} // namespace Teuchos

#endif // __Teuchos_MatrixMarket_SymmetrizingAdder_hpp
