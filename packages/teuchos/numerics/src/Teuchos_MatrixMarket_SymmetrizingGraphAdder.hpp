// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teuchos_MatrixMarket_SymmetrizingGraphAdder_hpp
#define __Teuchos_MatrixMarket_SymmetrizingGraphAdder_hpp

#include <Teuchos_as.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <string>


// Macro that marks a function as "possibly unused," in order to
// suppress build warnings.
#if ! defined(TRILINOS_UNUSED_FUNCTION)
#  if defined(__GNUC__) || (defined(__INTEL_COMPILER)  && !defined(_MSC_VER))
#    define TRILINOS_UNUSED_FUNCTION __attribute__((__unused__))
#  elif defined(__clang__)
#    if __has_attribute(unused)
#      define TRILINOS_UNUSED_FUNCTION __attribute__((__unused__))
#    else
#      define TRILINOS_UNUSED_FUNCTION
#    endif // Clang has 'unused' attribute
#  elif defined(__IBMCPP__)
// IBM's C++ compiler for Blue Gene/Q (V12.1) implements 'used' but not 'unused'.
//
// http://pic.dhe.ibm.com/infocenter/compbg/v121v141/index.jsp
#    define TRILINOS_UNUSED_FUNCTION
#  else // some other compiler
#    define TRILINOS_UNUSED_FUNCTION
#  endif
#endif // ! defined(TRILINOS_UNUSED_FUNCTION)


namespace Teuchos {
  namespace MatrixMarket {

    /// \class SymmetrizingGraphAdder
    /// \author Alicia Klinvex
    /// \brief Adds entries with optional symmetry to a sparse graph.
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
    class SymmetrizingGraphAdder {
    public:
      //! The type of indices of the sparse graph
      typedef typename AdderType::index_type index_type;

      /// \brief Constructor
      ///
      /// \param adder [in/out] The wrapped AdderType instance
      ///
      /// \param symmType [in] Canonical Matrix Market string
      ///   representing the symmetry storage type of the matrix data.
      SymmetrizingGraphAdder (const Teuchos::RCP<AdderType>& adder,
                         const std::string& symmType) :
        adder_ (adder),
        symmetrize_ (needsSymmetrization (symmType)),
        skew_ (isSkew (symmType))
      {}

      //! Add entry (i,j), and optionally symmetrize.
      void
      operator() (const index_type i,
                  const index_type j)
      {
        AdderType& theAdder = *adder_;

        theAdder (i, j);
        if (symmetrize_ && i != j) {
          // The optional third argument (which defaults to true)
          // specifies whether or not to count the entry against the
          // total expected number of entries.  We don't want to count
          // this entry because it wasn't part of the original data;
          // we inserted it because the caller doesn't want symmetric
          // storage.  The original data's total expected number of
          // entries only counts the entries that are in the original
          // data, not those that we insert.
          theAdder (j, i, false);
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
      //! Whether to negate when symmetrizing.
      bool skew_;
    };

  } // namespace MatrixMarket
} // namespace Teuchos

#endif // __Teuchos_MatrixMarket_SymmetrizingGraphAdder_hpp
