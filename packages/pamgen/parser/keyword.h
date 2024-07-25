// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef keywordH
#define keywordH

#include "pamgen_code_types.h"
#include "token.h"
#include <strings.h>
#include "token_stream.h"
#include <cstring>

// This file supports the keyword model of a token stream.  Under this
// model, the token stream is a sequence of keywords (of type TK_IDENTIFIER)
// each of which is optionally followed by additional tokens.  Each valid
// keyword has associated with it a parsing function that parses the
// additional tokens, if any, and which takes an optional integer argument.
// The latter can be used to enable the same parse function to parse closely 
// related keywords and their associated tokens.

// Standard signature for a keyword-associated parse function.
/* typedef Token (*Keyword_Function)(PAMGEN_NEVADA::Token_Stream *, int); */

// Structure representing a keyword.  The name must be a valid TK_IDENTIFIER 
// token.  The argument may have any value.  It is passed into the keyword
// function as its second argument and can be used to distinguish related
// keywords that are parsed by the same keyword function.  The module is
// a unique identifier indicating where the keyword came from.  It is 
// usually a class name but could be a developer or institution name.




extern "C" int PAMGEN_Keyword_Compare(const void *data1, const void *data2);
  // This is a qsort-compatible comparison function for keywords.
  // The two arguments must point to valid const Keyword structs.

extern "C" int PAMGEN_Cstring_Keyword_Compare(const void *token, const void *key);
  // This is a bsearch-compatible comparison function.  The first pointer
  // must point to a valid const char* and the second must point to a valid
  // const Keyword struct.


namespace PAMGEN_NEVADA{

typedef Token (*Keyword_Function)(Token_Stream *, int);


struct  Keyword {
  const char *name;
  int argument;
  /*   Token (*func) (PAMGEN_NEVADA::Token_Stream *,int); */
  Keyword_Function func; 
};

// Compare two keywords for strict equality.  In fact, a keyword table will
// not behave properly if it has two keyword with identical names, so for
// many purposes comparing the names is sufficient. 
inline bool operator==(const Keyword &a, const Keyword &b)
{
  assert(a.name != 0);
  assert(a.func != 0);
  assert(b.name != 0);
  assert(b.func != 0);
  return !strcmp(a.name, b.name) && a.argument==b.argument && a.func==b.func;
}
inline bool operator!=(const Keyword &a, const Keyword &b)
{
  assert(a.name != 0);
  assert(a.func != 0);
  assert(b.name != 0);
  assert(b.func != 0);
  return strcmp(a.name, b.name) || a.argument!=b.argument || a.func!=b.func;
}

inline bool operator>(const Keyword &a, const Keyword &b)
  // This is an STL-compatible ordering function for keywords.
{
  return PAMGEN_Keyword_Compare(&a, &b)>0;
}

inline bool operator<(const Keyword &a, const Keyword &b)
  // This is an STL-compatible ordering function for keywords.
{
  return PAMGEN_Keyword_Compare(&a, &b)<0;
}

}



#endif
