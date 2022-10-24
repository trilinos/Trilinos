// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef stk_expreval_lexer_hpp
#define stk_expreval_lexer_hpp

#include <vector>
#include <string>

namespace stk {
namespace expreval {

template <class T>
T convert_cast(const std::string & s);

enum Token {
  TOKEN_PLUS,
  TOKEN_MINUS,
  TOKEN_MULTIPLY,
  TOKEN_DIVIDE,
  TOKEN_PERCENT,
  TOKEN_QUESTION,
  TOKEN_COMMA,
  TOKEN_COLON,
  TOKEN_SEMI,
  TOKEN_LPAREN,
  TOKEN_RPAREN,
  TOKEN_LBRACK,
  TOKEN_RBRACK,
  TOKEN_ASSIGN,
  TOKEN_LESS,
  TOKEN_GREATER,
  TOKEN_ARITHMETIC_OR,
  TOKEN_ARITHMETIC_AND,
  TOKEN_NOT,
  TOKEN_EQUAL,
  TOKEN_NOT_EQUAL,
  TOKEN_GREATER_EQUAL,
  TOKEN_LESS_EQUAL,
  TOKEN_LOGICAL_OR,
  TOKEN_LOGICAL_AND,
  TOKEN_IDENTIFIER,
  TOKEN_REAL_CONSTANT,
  TOKEN_INTEGER_CONSTANT,
  TOKEN_LITERAL,
  TOKEN_EXPONENTIATION,
  TOKEN_END
};

class Lexem
{
public:
  Lexem(Token token, const char * from, const char * to)
    : m_token(token),
      m_value(from, to)
  {}

  Lexem(Token token, const char * value)
    : m_token(token),
      m_value(value)
  {}

  Token getToken() const { return m_token; }

  const std::string &getString() const { return m_value; }

  template<class T>
  T getValue() const { return convert_cast<T>(m_value); }

private:
  Token m_token;
  std::string m_value;
};

typedef std::vector<Lexem> LexemVector;

LexemVector tokenize(const std::string & expression);

} // namespace expreval
} // namespace stk

#endif // stk_expreval_lexer_hpp
