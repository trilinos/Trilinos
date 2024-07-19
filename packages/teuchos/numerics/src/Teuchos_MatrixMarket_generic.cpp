// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_MatrixMarket_generic.hpp"
#include "Teuchos_MatrixMarket_split.hpp"
#include <algorithm>

namespace Teuchos {
  namespace MatrixMarket {

    int maxLineLength() { return 1024; }

    bool
    checkCommentLine (const std::string& line,
                      size_t& start,
                      size_t& size,
                      const size_t lineNumber,
                      const bool tolerant,
                      const bool maybeBannerLine)
    {
      // In tolerant mode, empty lines are considered comment lines.
      if (line.empty ()) {
        if (tolerant) {
          return true;
        }
        else {
          std::ostringstream os;
          os << "Line " << lineNumber << " contains no characters";
          throw std::invalid_argument (os.str());
        }
      }
      // The line of comments or data "starts" after any whitespace
      // characters.  Whitespace-only lines are considered "empty."
      start = line.find_first_not_of (" \t");
      if (start == std::string::npos) {
        // It's a whitespace-only line.  We consider those comments
        // in tolerant mode, syntax errors otherwise.
        if (tolerant) {
          return true;
        }
        else {
          std::ostringstream os;
          os << "Line " << lineNumber << " contains only whitespace";
          throw std::invalid_argument (os.str());
        }
      }
      // Position of the first comment character (if any), relative to
      // the first non-whitespace character in the line.  (If we got
      // this far, then the line has at least one non-whitespace
      // character.)
      const size_t commentPos = line.find_first_of("%#", start);
      if (commentPos == std::string::npos) {
        // There are no comment characters in the line.
        // line.substr(start,npos) gives the substring of line
        // containing valid data.
        size = std::string::npos;
        return false;
      }
      else if (commentPos == start) {
        // The line has 0 or more whitespace characters, followed by a
        // start-of-comment character.  However, the Matrix Market
        // banner line starts with "%%MatrixMarket", so we have to
        // look for this, if the caller allows this.
        if (maybeBannerLine) {
          const size_t bannerStart =
            line.substr (commentPos).find ("%%MatrixMarket");
          if (bannerStart != std::string::npos) { // It's a banner line!
            size = line.size() - commentPos;
            return false;
          }
          else { // It's a comment line.  Ah well.
            size = 0;
            return true;
          }
        }
        else {
          size = 0;
          return true;
        }
      }
      else {
        // [start, start+size-1] is the (inclusive) range of
        // characters (if any) between the first non-whitespace
        // character, and the first comment character.  That range
        // could contain valid data, so we don't consider this a
        // "comment line."
        size = commentPos - start;
        return false;
      }
    }

  } // namespace MatrixMarket
} // namespace Teuchos
