// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STRUTILS_H
#define ROL_STRUTILS_H

#include "Teuchos_ConfigDefs.hpp"

#include <string>
#include <vector>

namespace ROL {


  /**
   * \brief Provides std::string manipulation utilities that are not provided in the
   * standard C++ std::string class.
   */
class TEUCHOSCORE_LIB_DLL_EXPORT StrUtils
{
public:
  /** \brief Read a file, putting each line into a std::string */
  static std::vector<std::string> readFile(std::istream& is, char comment);

  /** \brief Split an input std::string that contains newlines into an array
      of strings, one for each line */
  static std::vector<std::string> splitIntoLines(const std::string& input);

  /** \brief Tokenize a file into whitespace-delimited tokens */
  static std::vector<std::vector<std::string> > tokenizeFile(std::istream& is, char comment);

  /** \brief Read a single line into a std::string */
  static bool readLine(std::istream& is, std::string& line);

  /** \brief . */
  static std::vector<std::string> stringTokenizer(const std::string& str);

  /** \brief . */
  static std::vector<std::string> getTokensPlusWhitespace(const std::string& str);

  /** \brief . */
  static std::string reassembleFromTokens(const std::vector<std::string>& tokens, int iStart=0);

  /** \brief . */
  static void splitList(const std::string& bigstring, std::vector<std::string>& elements);

  /** \brief . */
  static int findNextWhitespace(const std::string& str, int offset);

  /** \brief . */
  static int findNextNonWhitespace(const std::string& str, int offset);

  /** \brief . */
  static  std::string varSubstitute(const std::string& rawLine,
    const std::string& varName,
    const std::string& varValue);

  /** \brief . */
  static std::string varTableSubstitute(const std::string& rawLine,
    const std::vector<std::string>& varNames,
    const std::vector<std::string>& varValues);

  /** \brief . */
  static std::string envSubstitute(const std::string& line);

  /** \brief Find the substring before a specified substring. For example,
   * before("abcdefghij", "gh") returns "abcdef". */
  static std::string before(const std::string& str, const std::string& sub);

  /** \brief Find the substring before a specified character. For example,
   * before("abcdefghij", 'g') returns "abcdef". */
  static std::string before(const std::string& str, char sub);

  /** \brief Find the substring after a specified substring. For example,
   * before("abcdefghij", "gh") returns "ij". */
  static std::string after(const std::string& str, const std::string& sub);

  /** \brief Find the position at which a substring first occurs. For example,
   * find("abcdefghij", "gh") returns 6. */
  static int find(const std::string& str, const std::string& sub);

  /** \brief Returns true if a std::string consists entirely of whitespace */
  static bool isWhite(const std::string& str);

  /** \brief Convert unprintable non-null characters to whitespace */
  static std::string fixUnprintableCharacters(const std::string& str);

  /** \brief Returns true if a std::string has any non-whitespace */
  static bool isNonWhite(const std::string& str) {return !isWhite(str);}

  /** \brief Returns the std::string between two delimiting strings, and returns
   * by reference the strings before and after the delimiters.
   *
   * For example, between("abcdefghij", "c", "g", front, back)
   * returns "def" and sets front to "ab", back to "hij". */
  static std::string between(const std::string& str, const std::string& begin,
    const std::string& end, std::string& front, std::string& back);

  /** \brief Returns the substring between two positions.
   *
   *  For example, subString("abcdefghij", 2, 5) returns "cde". */
  static std::string subString(const std::string& str, int begin, int end);

  static std::string readFromStream(std::istream& is);

  /** \brief Converts a std::string to all upper case */
  static std::string allCaps(const std::string& str);

  /** \brief Returns the double value of a std::string. */
  static double atof(const std::string& str);

  /** \brief Returns the int value of a std::string. */
  static int atoi(const std::string& str);

  /** \brief Print lines with prefix first. */
  static std::ostream& printLines(
    std::ostream             &os
    ,const std::string       &linePrefix
    ,const std::string       &lines
    );
	
  /** \brief Removes all the spaces in a string. */
  static std::string removeAllSpaces(std::string stringToClean);

};


} // namespace ROL


#endif
