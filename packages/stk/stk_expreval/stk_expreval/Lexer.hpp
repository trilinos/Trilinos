// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
T convert_cast(const std::string &s);

/**
 * Enumeration <b>Tokens</b> defines the tokens returned by the lexer.
 *
 */
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


/**
 * Class <b>Lexem</b> associates a token with a simple expression of the type
 * described by the token.
 *
 */
class Lexem
{
public:
  /**
   * Creates a new <b>Lexem</b> instance.
   *
   * @param token		an <b>int</b> value of the token which describes the
   *				value.
   *
   * @param from		a <b>char</b> const pointer to the beginning of the
   *				value.
   *
   * @param to			a <b>char</b> const pointer to the end of the
   *				value.
   */
  Lexem(Token token, const char *from, const char *to)
    : m_token(token),
      m_value(from, to)
  {}

  /**
   * Creates a new <b>Lexem</b> instance.
   *
   * @param token		an <b>int</b> value of the token which describes the
   *				value.
   *
   * @param value		a <b>char</b> const pointer to a c-style string of
   *				the value.
   */
  Lexem(Token token, const char *value)
    : m_token(token),
      m_value(value)
  {}

  /**
   * Member function <b>getToken</b> returns the token.
   *
   * @return			a <b>Token</b>.
   */
  Token getToken() const {
    return m_token;
  }

  /**
   * Member function <b>getValue</b> returns const reference to the value associated
   * with the token.
   *
   * @return			a <b>std::string</b> const reference to the value
   *				associated with the token.
   */
  const std::string &getString() const {
    return m_value;
  }

  /**
   * Member function <b>getValue</b> returns const reference to the value associated
   * with the token.
   *
   * @return			a <b>std::string</b> const reference to the value
   *				associated with the token.
   */
  template<class T>
  T getValue() const {
    return convert_cast<T>(m_value);
  }

private:
  Token			m_token;		///< Token which describes the value
  std::string		m_value;		///< Value
};

typedef std::vector<Lexem> LexemVector;

LexemVector tokenize(const std::string &expression);

} // namespace expreval
} // namespace stk

#endif // stk_expreval_lexer_hpp
