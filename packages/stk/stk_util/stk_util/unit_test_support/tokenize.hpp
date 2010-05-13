/**********************************************************************
 * tokenize() - Tokenizes a string into parts by separators.          *
 * Copyright (C) 2000  C. Brandon Forehand                            *
 * <b4hand@users.sourceforge.net>                                     *
 *                                                                    *
 * This code is free software; you can redistribute it and/or         *
 * modify it under the terms of the GNU General Public License        *
 * as published by the Free Software Foundation; either version 2     *
 * of the License, or (at your option) any later version.             *
 *                                                                    *
 * This program is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      *
 * GNU General Public License for more details.                       *
 **********************************************************************/

#ifndef STK_IO_UTIL_TOKENIZE_H
#define STK_IO_UTIL_TOKENIZE_H

#include <algorithm>
#include <functional>
#include <string>
#include <vector>

namespace stk {
  namespace io {
    namespace util {
      //-----------------------------------------------------------------------------
      /** This is the default predicate for the tokenize function.  It is necessary
	  because the prototype for isspace() is int isspace(int). */
      //----------------------------------------------------------------------------
      struct is_space : public std::unary_function<char,bool>
      {
	bool operator() (char c) const
	{ return static_cast<bool>(isspace(c)); }
      };

      struct is_punct : public std::unary_function<char,bool>
      {
	bool operator() (char c) const
	{ return static_cast<bool>(isspace(c) || ispunct(c)); }
      };

      //----------------------------------------------------------------------------
      /** Powerful general predicate for tokenize().  Predicate returns true
	  for all arguments that are contained in the string passed in at
	  construction. For example, recognize(" \v\f\n\r\t") is the same as
	  is_space() above. */
      //----------------------------------------------------------------------------
      class recognize : public std::unary_function<char,bool>
      {
      public:
	explicit recognize(const std::string &str) : mStr (str) {}
	bool operator() (char c) const
	{ return (mStr.end() != std::find(mStr.begin(),mStr.end(),c)); }

      private:
	std::string mStr;
      };

      //----------------------------------------------------------------------------
      /** Takes a string and a predicate, and tokenizes the string according to the
	  separators which are defined by the provided predicate.  The predicate
	  should evaluate to true when applied to a separator. */
      //----------------------------------------------------------------------------
      template <class Pred>
      inline std::vector<std::string> tokenize (const std::string &s, Pred p)
      {
	using namespace std;

	vector<string> result;
	string::const_iterator i = s.begin();
	string::const_iterator tokenEnd = s.begin();

	while (i != s.end())
	  {
	    // Eat seperators
	    while (p(*i))
	      i++;

	    // Find next token
	    tokenEnd = find_if(i,s.end(),p);

	    // Append token to result
	    if (i != tokenEnd)
	      result.push_back(string(i,tokenEnd));

	    i = tokenEnd;
	  }

	return result;
      }

      //----------------------------------------------------------------------------
      /** Tokenizes a string.  Same as above, but uses is_space() as default
	  predicate.  For some reason, g++ won't recognize a default parameter. */
      //----------------------------------------------------------------------------
      inline std::vector<std::string> tokenize (const std::string &s)
      {
	return tokenize(s,is_space());
      }
    }
  }
}
#endif

#if 0
#include <iostream>
using std::cout;
using std::cin;

typedef std::vector<std::string> TokenList;

int main()
{
  char s[128];
  while(!cin.eof()) {
    cout << "Enter a string: ";
    cin.getline(s,128);
    std::string input_line(s);
    if (input_line != "quit") {
      std::vector<std::string> tokens = tokenize(input_line);
      cout << "There were " << tokens.size() << " tokens in the line\n";
      TokenList::const_iterator I = tokens.begin();
      while (I != tokens.end()) {
	cout << "'" << *I++ << "'\t";
      }
      cout << '\n';
    } else {
      exit(0);
    }
  }
}
#endif
