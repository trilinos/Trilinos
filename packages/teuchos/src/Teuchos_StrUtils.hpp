// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_STRUTILS_H
#define TEUCHOS_STRUTILS_H

/*! \file Teuchos_StrUtils.hpp
    \brief A string utilities class for Teuchos
*/

#include "Teuchos_ConfigDefs.hpp"

#include "Teuchos_Utils.hpp"
#include "Teuchos_Array.hpp"

namespace Teuchos
{
  /** 
   * \brief Provides string manipulation utilities that are not provided in the
   * standard C++ string class.
   */
  class StrUtils
    {
    public:
      /** \brief Read a file, putting each line into a string */
      static Array<string> readFile(istream& is, char comment);

      /** \brief Split an input string that contains newlines into an array
 	 of strings, one for each line */
      static Array<string> splitIntoLines(const string& input);

      /** \brief Tokenize a file into whitespace-delimited tokens */
      static Array<Array<string> > tokenizeFile(istream& is, char comment);

      /** \brief Read a single line into a string */
      static bool readLine(istream& is, string& line);

      static Array<string> stringTokenizer(const string& str);

      static Array<string> getTokensPlusWhitespace(const string& str);

      static string reassembleFromTokens(const Array<string>& tokens, int iStart=0);

      static void splitList(const string& bigstring, Array<string>& elements);

      static int findNextWhitespace(const string& str, int offset);

      static int findNextNonWhitespace(const string& str, int offset);


      static  string varSubstitute(const string& rawLine,
                                   const string& varName,
                                   const string& varValue);

      static string varTableSubstitute(const string& rawLine,
                                       const Array<string>& varNames,
                                       const Array<string>& varValues);

      static string envSubstitute(const string& line);

      /** \brief Find the substring before a specified substring. For example,
       * before("abcdefghij", "gh") returns "abcdef". */
      static string before(const string& str, const string& sub);

      /** \brief Find the substring before a specified character. For example,
       * before("abcdefghij", 'g') returns "abcdef". */
      static string before(const string& str, char sub);

      /** \brief Find the substring after a specified substring. For example,
       * before("abcdefghij", "gh") returns "ij". */
      static string after(const string& str, const string& sub);

      /** \brief Find the position at which a substring first occurs. For example,
       * find("abcdefghij", "gh") returns 6. */
      static int find(const string& str, const string& sub);

      /** \brief Returns true if a string consists entirely of whitespace */
      static bool isWhite(const string& str);

      /** \brief Convert unprintable non-null characters to whitespace */
      static string fixUnprintableCharacters(const string& str);

      /** \brief Returns true if a string has any non-whitespace */
      static bool isNonWhite(const string& str) {return !isWhite(str);}

      /** \brief Returns the string between two delimiting strings, and returns
       * by reference the strings before and after the delimiters.
       *
       * For example, between("abcdefghij", "c", "g", front, back)
       * returns "def" and sets front to "ab", back to "hij". */
      static string between(const string& str, const string& begin,
                            const string& end, string& front, string& back);

      /** \brief Returns the substring between two positions. 
       *
       *  For example, subString("abcdefghij", 2, 5) returns "cde". */
      static string subString(const string& str, int begin, int end);

      static string readFromStream(istream& is);

      /** \brief Converts a string to all upper case */
      static string allCaps(const string& str);

      /** \brief Returns the double value of a string. */
      static double atof(const string& str);

      /** \brief Returns the int value of a string. */
      static int atoi(const string& str);
    };



}

#endif
