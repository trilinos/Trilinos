#ifndef STRUTILS_H
#define STRUTILS_H

#include "Teuchos_ConfigDefs.hpp"

#include "Teuchos_Utils.hpp"
#include "Teuchos_Array.hpp"

namespace Teuchos
{
  /** \ingroup General
   * Some string manipulation utilities that are not provided in the
   * standard C++ string class.
   */

  class StrUtils
    {
    public:
      /** read a file, putting each line into a string */
      static Array<string> readFile(istream& is, char comment);

      /** split an input string that contains newlines into an array
       * of strings, one for each line */
      static Array<string> splitIntoLines(const string& input);

      /** tokenize a file into whitespace-delimited tokens */
      static Array<Array<string> > tokenizeFile(istream& is, char comment);

      /** read a single line into a string */
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

      /** find the substring before a specified substring. For example,
       * before("abcdefghij", "gh") returns "abcdef". */
      static string before(const string& str, const string& sub);

      /** find the substring before a specified character. For example,
       * before("abcdefghij", 'g') returns "abcdef". */
      static string before(const string& str, char sub);

      /** find the substring after a specified substring. For example,
       * before("abcdefghij", "gh") returns "ij". */
      static string after(const string& str, const string& sub);

      /** find the position at which a substring first occurs. For example,
       * find("abcdefghij", "gh") returns 6. */
      static int find(const string& str, const string& sub);

      /** returns true if a string consists entirely of whitespace */
      static bool isWhite(const string& str);

      /** convert unprintable non-null characters to whitespace */
      static string fixUnprintableCharacters(const string& str);

      /** returns true if a string has any non-whitespace */
      static bool isNonWhite(const string& str) {return !isWhite(str);}

      /** returns the string between two delimiting strings, and returns
       * by reference the strings before and after the delimiters.
       * For example, between("abcdefghij", "c", "g", front, back)
       * returns "def" and sets front to "ab", back to "hij". */
      static string between(const string& str, const string& begin,
                            const string& end, string& front, string& back);

      /** returns the substring between two positions. For example,
          subString("abcdefghij", 2, 5) returns "cde". */
      static string subString(const string& str, int begin, int end);

      static string readFromStream(istream& is);

      /** converts a string to all upper case */
      static string allCaps(const string& str);

      /** returns the double value of a string. */
      static double atof(const string& str);

      /** returns the int value of a string. */
      static int atoi(const string& str);
    };



}

#endif
