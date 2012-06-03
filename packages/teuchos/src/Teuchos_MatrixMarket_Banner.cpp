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

#include "Teuchos_MatrixMarket_Banner.hpp"
#include "Teuchos_MatrixMarket_split.hpp"
#include "Teuchos_TestForException.hpp"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>


namespace Teuchos {
  namespace MatrixMarket {

    using details::split;
    using details::trim_and_lowercase;

    std::string
    Banner::validateObjectType (const std::string& objectType, const bool tolerant)
    {
      // Canonical representation is lowercase
      std::string out = trim_and_lowercase (objectType);

      const char* const validValues[] = {"matrix"};
      const int numValidValues = 1;
      if (tolerant)
        // This is the only value currently defined for this token in
        // the Matrix Market format, so we just return it.
        return std::string (validValues[0]);
      else if (validValues + numValidValues ==
               std::find (validValues, validValues + numValidValues, out))
        throw std::invalid_argument("Object type \"" + out + "\" is "
                                    "not one of the valid values");
      else
        return out;
    }

    std::string
    Banner::validateMatrixType (const std::string& matrixType, const bool tolerant)
    {
      // Canonical representation is lowercase
      std::string out = trim_and_lowercase (matrixType);

      const char* const validValues[] = {"coordinate", "array"};
      const int numValidValues = 2;
      if (validValues + numValidValues == std::find (validValues, validValues + numValidValues, out))
        throw std::invalid_argument("Matrix type \"" + out + "\" is not one of the valid values");
      else
        return out;
    }

    std::string
    Banner::validateDataType (const std::string& dataType, const bool tolerant)
    {
      // Canonical representation is lowercase
      std::string out = trim_and_lowercase (dataType);

      const char* const validValues[] = {"real", "complex", "integer", "pattern"};
      const int numValidValues = 4;
      if (validValues + numValidValues == std::find (validValues, validValues + numValidValues, out))
        throw std::invalid_argument("Data type \"" + out + "\" is not one of the valid values");
      else
        return out;
    }

    std::string
    Banner::validateSymmType (const std::string& symmType, const bool tolerant)
    {
      // Canonical representation is lowercase
      std::string out = trim_and_lowercase (symmType);

      if (tolerant)
        {
          const char* const validValues[] =
            {"general", "nonsymmetric", "unsymmetric", "symmetric",
             "skew-symmetric", "skew", "hermitian"};
          const int numValidValues = 7;
          if (validValues + numValidValues == std::find (validValues, validValues + numValidValues, out))
            throw std::invalid_argument("Symmetry type \"" + out + "\" is not one of the valid values");
          else
            {
              if (out == "nonsymmetric" || out == "unsymmetric")
                return std::string("general");
              else if (out == "skew")
                return std::string("skew-symmetric");
              else
                return out;
            }
        }
      else
        {
          const char* const validValues[] = {"general", "symmetric", "skew-symmetric", "hermitian"};
          const int numValidValues = 4;
          if (validValues + numValidValues == std::find (validValues, validValues + numValidValues, out))
            throw std::invalid_argument("Symmetry type \"" + out + "\" is not one of the valid values");
          else
            return out;
        }
    }


    void
    Banner::setDefaults (const int howMany)
    {
      if (howMany >= 4)
        objectType_ = "matrix";
      if (howMany >= 3)
        matrixType_ = "coordinate";
      if (howMany >= 2)
        dataType_ = "real";
      if (howMany >= 1)
        symmType_ = "general";
    }

    Banner::Banner (const std::string& line, const bool tolerant)
    {
      size_t start;

      if (line.empty()) {
        if (tolerant) {
          setDefaults (4);
          return;
        }
        else {
          throw std::invalid_argument ("The banner line is empty");
        }
      }
      start = line.find_first_not_of (" \t");
      if (start == std::string::npos) {
        if (tolerant) {
          setDefaults (4);
          return;
        }
        else {
          throw std::invalid_argument ("The banner line contains only "
                                       "whitespace characters");
        }
      }
      else if (start != 0 && ! tolerant) {
        // If tolerant, we allow the banner line to start with
        // whitespace characters, and keep reading.
        throw std::invalid_argument ("The banner line is not allowed to start "
                                     "with whitespace characters");
      }

      // Find "%%MatrixMarket" -- it should be the first thing in the
      // banner line, and all other data should come after it.
      // Optionally relax to allow any case, and possibly a space
      // between "%%" and "MatrixMarket".
      size_t ppStart = line.find ("%%", start);
      size_t tokenStart;
      if (ppStart == std::string::npos) {
        if (tolerant) {
          tokenStart = start; // Just ignore the missing %%
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "The Matrix "
            "Market file's banner line should always start with \"%%\".  Here "
            "is the offending line: " << std::endl << line);
        }
      }
      else {
        tokenStart = ppStart + 2;
        if (tokenStart >= line.size()) {
          // There's no banner information after the %%.
          if (tolerant) {
            setDefaults (4);
            return;
          }
          else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "The Matrix "
              "Market file's banner line needs to contain information after the "
              "\"%%\" marker.  Here is the offending line: " << std::endl << line);
          }
        }
      }
      //
      // In tolerant mode, fill in missing tokens with their default
      // values.
      //
      // After extracting the %%, search for the five tokens.
      std::vector<std::string> tokens = split (line, " \t", 2);
      const int numTokens = tokens.size();
      if (numTokens < 1) {
        if (tolerant) {
          setDefaults (4);
          return;
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "The Matrix "
            "Market file's banner line must always begin with the \"Matrix"
            "Market\" keyword.  Here is the offending line: " << std::endl
            << line);
        }
      }
      // In tolerant mode, just ignore the first token.
      if (! tolerant && tokens[0] != "MatrixMarket") {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "The Matrix "
          "Market file's banner line must always begin with the \"Matrix"
          "Market\" keyword.  Here is the offending line: " << std::endl
          << line);
      }
      if (numTokens < 5) {
        if (tolerant) {
          setDefaults (5 - numTokens); // how many defaults to set
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "The Matrix "
            "Market file's banner line must always have 5 tokens, but yours "
            "only has " << numTokens << "token" << (numTokens != 1 ? "s" : "")
            << ".  Here is the offending line: " << std::endl << line);
        }
      }
      if (numTokens >= 2) {
        objectType_ = validateObjectType (tokens[1], tolerant);
      }
      if (numTokens >= 3) {
        matrixType_ = validateMatrixType (tokens[2], tolerant);
      }
      if (numTokens >= 4) {
        dataType_ = validateDataType (tokens[3], tolerant);
      }
      if (numTokens >= 5) {
        symmType_ = validateSymmType (tokens[4], tolerant);
      }
    }

    std::ostream&
    operator<< (std::ostream& out, const Banner& banner)
    {
      out << "%%MatrixMarket"
          << " " << banner.objectType()
          << " " << banner.matrixType()
          << " " << banner.dataType()
          << " " << banner.symmType();
      return out;
    }
  } // namespace MatrixMarket
} // namespace Teuchos

