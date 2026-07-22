// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_MatrixMarket_split.hpp"
#include <algorithm>
#include <cctype> // tolower
#include <utility> // std::pair


namespace Teuchos {
  namespace MatrixMarket {
    namespace details {

      std::string
      trim (const std::string& in)
      {
        size_t start = in.find_first_not_of(" \t");
        size_t end = in.find_last_not_of(" \t");
        if (start == std::string::npos)
          return std::string("");
        else
          return in.substr (start, end-start+1);
      }

      std::string
      lowercase (const std::string& in)
      {
        std::string out (in);
        std::transform (in.begin(), in.end(), out.begin(), tolower);
        return out;
      }

      std::string
      trim_and_lowercase (const std::string& in)
      {
        std::string out = trim (in);
        std::string out2 (out);
        std::transform (out.begin(), out.end(), out2.begin(), tolower);
        return out2;
      }

      //! Does \c range signify "no more tokens"?
      static bool
      endToken (const std::pair<size_t, size_t>& range)
      {
        return range.first == std::string::npos && range.second == std::string::npos;
      }

      /// Get the next token, if there is one
      ///
      /// \return If there is a next token, a (startpos,length) pair
      ///   suitable for string::substr(); otherwise, npos,npos), which
      ///   alone satisfies the endToken() predicate.
      static std::pair<size_t, size_t>
      nextToken (const std::string& str,
                 const std::string& delimiters,
                 const size_t start,
                 const size_t size)
      {
        using std::make_pair;
        using std::string;

        if (start >= size)
          return make_pair (string::npos, string::npos);

        // First index of a non-delimiter character
        const size_t first = str.find_first_not_of (delimiters, start);
        if (first == string::npos)
          // There are only delimiter characters left
          return make_pair (string::npos, string::npos);
        else if (first == size-1)
          // There's only one non-delimiter character left
          return make_pair (first, 1);
        else
          { // Next index of a delimiter character
            // Search for the next delimiting token from "first" (instead of
            // "start+1" as originally coded and commented out below.)
            // const size_t next = str.find_first_of (delimiters, start+1);
            const size_t next = str.find_first_of (delimiters, first);
            return make_pair (first, next - first);
          }
      }

      std::vector<std::string>
      split (const std::string& str, const std::string& delimiters, const size_t start)
      {
        size_t curStart = start;
        size_t size = str.size();
        std::vector<std::string> tokens;
        while (true) {
          std::pair<size_t, size_t> token = nextToken (str, delimiters, curStart, size);
          if (endToken (token)) {
            break;
          }
          tokens.push_back (str.substr (token.first, token.second));
          curStart = token.first + token.second;
        }
        return tokens;
      }
    } // namespace details
  } // namespace MatrixMarket
} // namespace Teuchos
