// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teuchos_MatrixMarket_split_hpp
#define __Teuchos_MatrixMarket_split_hpp

#include <string>
#include <vector>


namespace Teuchos {
  namespace MatrixMarket {
    namespace details {

      //! Trim whitespace from both sides of the given string.
      std::string
      trim (const std::string& in);

      //! Return lowercase version of the given string.
      std::string
      lowercase (const std::string& in);

      //! Trim whitespace from both sides, and make lowercase.
      std::string
      trim_and_lowercase (const std::string& in);

      /// \brief Split the given string using the given set of delimiters.
      ///
      /// Split the string \c str, optionally starting at position \c
      /// start, into zero or more tokens separated by one or more of the
      /// given delimiter characters in \c delimiters.
      ///
      /// \param str [in] String to split into tokens
      /// \param delimiters [in] Array of one or more delimiter character(s)
      /// \param start [in] Position in \c str where the search should begin.
      ///   Defaults to zero.
      ///
      /// \return Vector of zero or more tokens, none of which contain any
      ///   of the delimiter character(s)
      std::vector<std::string>
      split (const std::string& str,
             const std::string& delimiters,
             const size_t start=0);

    } // namespace details
  } // namespace MatrixMarket
} // namespace Teuchos

#endif // __Teuchos_MatrixMarket_split_hpp
