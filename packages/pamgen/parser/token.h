// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef tokenH
#define tokenH
 
#include <stdio.h>
#include <setjmp.h>

#include "pamgen_code_types.h"
#include "token_enum.h"
#include "token_value.h"

#include <string>
#include <assert.h>

namespace PAMGEN_NEVADA{

/*****************************************************************************/
class Token
/*****************************************************************************/
// Represents a token -- the atom of input recognized by the scanner and
// passed to the parser.
 { 
  public:
    // The token type defaults to the exit token.
    Token();

    // In general, a token consists of a type and an associated 
    // semantic value
    Token(Token_Type t, Token_Value& sv); 

    // Not all tokens have an associated semantic value.
    // However, we do *not* want an implicit conversion from integral
    // types to Token, so we add a dummy argument.
    enum No_Value { no_value };
    Token(Token_Type t, No_Value);

    Token(const Token &src);

    ~Token();

    Token& operator=(const Token&);

    Token_Type  Type()          const {return type;}

    void Convert_Real(Real& x)    {s.fval = x;}

    int         As_Int()    const 
    {
      assert(Type()==TK_INTEGER);
      return s.ival;
    }

    Real        As_Real()   const 
    {
      assert(Type()==TK_REAL || Type()==TK_INTEGER);
      if(Type()==TK_INTEGER) return s.ival;
      //Type()==TK_REAL
      return s.fval;
    }

    const char* As_String() const 
    {
      assert(Type()==TK_IDENTIFIER || Type()==TK_STRING);
      return s.sval;
    }

    // Return the value of a quoted string token with its 
    // quotes stripped away. 
    std::string      As_Stripped_String() const;

    friend bool operator==(const Token &, const Token &);
    friend bool operator!=(const Token &, const Token &);

    friend bool operator==(const Token& tk, const char *c){ 
      return (int)tk.type == (int)TK_IDENTIFIER && !Token_Match(tk.s.sval, c);
    }
    friend bool operator==(const char *c, const Token& tk){ 
      return (int)tk.type == (int)TK_IDENTIFIER && !Token_Match(tk.s.sval, c);
    }
    friend bool operator!=(const Token& tk, const char *c){ 
      return (int)tk.type != (int)TK_IDENTIFIER || Token_Match(tk.s.sval, c);
    }
    friend bool operator!=(const char *c, const Token& tk){ 
      return (int)tk.type != (int)TK_IDENTIFIER || Token_Match(tk.s.sval, c);
    }

    // Matches two character strings according to the following rules:
    // Two character strings t and k match if each word of t matches
    // each word of k.
    // Two words match if word from t matches the first part of word 
    // from k.
    static int Token_Match(const char *t, const char *k);

  private:
    Token_Type type;
    Token_Value s;

    void copy(const Token&);
    void free_resources();
};
}//end namespace
#endif
