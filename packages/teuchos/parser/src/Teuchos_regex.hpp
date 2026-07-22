// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_REGEX_HPP
#define TEUCHOS_REGEX_HPP

#include <Teuchos_Language.hpp>
#include <Teuchos_FiniteAutomaton.hpp>
#include <Teuchos_ReaderTables.hpp>
#include <Teuchos_Reader.hpp>

namespace Teuchos {
namespace regex {

enum {
  PROD_REGEX,
  PROD_UNION_DECAY,
  PROD_UNION,
  PROD_CONCAT_DECAY,
  PROD_CONCAT,
  PROD_QUAL_DECAY,
  PROD_STAR,
  PROD_PLUS,
  PROD_MAYBE,
  PROD_SINGLE_CHAR,
  PROD_ANY,
  PROD_SINGLE_SET,
  PROD_PARENS_UNION,
  PROD_SET_POSITIVE,
  PROD_SET_NEGATIVE,
  PROD_POSITIVE_SET,
  PROD_NEGATIVE_SET,
  PROD_SET_ITEMS_DECAY,
  PROD_SET_ITEMS_ADD,
  PROD_SET_ITEM_CHAR,
  PROD_SET_ITEM_RANGE,
  PROD_RANGE
};

enum { NPRODS = PROD_RANGE + 1 };

enum {
  TOK_CHAR,
  TOK_DOT,
  TOK_LRANGE,
  TOK_RRANGE,
  TOK_LPAREN,
  TOK_RPAREN,
  TOK_UNION,
  TOK_RANGE,
  TOK_NEGATE,
  TOK_STAR,
  TOK_PLUS,
  TOK_MAYBE
};

enum { NTOKS = TOK_MAYBE + 1 };

Language make_language();
LanguagePtr ask_language();

void make_lexer(FiniteAutomaton& result);

ReaderTablesPtr ask_reader_tables();

void make_dfa(FiniteAutomaton& result, std::string const& name, std::string const& regex, int token);

class Reader : public Teuchos::Reader {
 public:
  Reader(int result_token_in);
  virtual ~Reader() {}
 protected:
  virtual void at_shift(any& result, int token, std::string& text);
  virtual void at_reduce(any& result, int token, std::vector<any>& rhs);
 private:
  int result_token;
};

}  // end namespace regex
}  // end namespace Teuchos

#endif
