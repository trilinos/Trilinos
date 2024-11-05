// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_MATHEXPR_HPP
#define TEUCHOS_MATHEXPR_HPP

#include <set>

#include <Teuchos_Language.hpp>
#include <Teuchos_Reader.hpp>

namespace Teuchos {

namespace MathExpr {

enum {
  PROD_PROGRAM,
  PROD_NO_STATEMENTS,
  PROD_NEXT_STATEMENT,
  PROD_ASSIGN,
  PROD_NO_EXPR,
  PROD_YES_EXPR,
  PROD_EXPR,
  PROD_TERNARY_DECAY,
  PROD_OR_DECAY,
  PROD_AND_DECAY,
  PROD_ADD_SUB_DECAY,
  PROD_MUL_DIV_DECAY,
  PROD_NEG_DECAY,
  PROD_POW_DECAY,
  PROD_TERNARY,
  PROD_OR,
  PROD_AND,
  PROD_GT,
  PROD_LT,
  PROD_GEQ,
  PROD_LEQ,
  PROD_EQ,
  PROD_BOOL_PARENS,
  PROD_ADD,
  PROD_SUB,
  PROD_MUL,
  PROD_DIV,
  PROD_POW,
  PROD_CALL,
  PROD_NO_ARGS,
  PROD_SOME_ARGS,
  PROD_FIRST_ARG,
  PROD_NEXT_ARG,
  PROD_NEG,
  PROD_VAL_PARENS,
  PROD_CONST,
  PROD_VAR,
  PROD_NO_SPACES,
  PROD_SPACES
};

enum { NPRODS = PROD_SPACES + 1 };

enum {
  TOK_SPACE,
  TOK_NAME,
  TOK_ADD,
  TOK_SUB,
  TOK_MUL,
  TOK_DIV,
  TOK_POW,
  TOK_LPAREN,
  TOK_RPAREN,
  TOK_COMMA,
  TOK_CHECK,
  TOK_CHOOSE,
  TOK_GT,
  TOK_LT,
  TOK_GEQ,
  TOK_LEQ,
  TOK_EQ,
  TOK_AND,
  TOK_OR,
  TOK_CONST,
  TOK_SEMICOLON,
  TOK_ASSIGN
};

enum { NTOKS = TOK_ASSIGN + 1 };

Language make_language();

LanguagePtr ask_language();

ReaderTablesPtr ask_reader_tables();

class SymbolSetReader : public Reader {
 public:
  SymbolSetReader();
  virtual ~SymbolSetReader();
 public:
  std::set<std::string> variable_names;
  std::set<std::string> function_names;
 private:
  virtual void at_shift(any& result, int token, std::string& text);
  virtual void at_reduce(any& result, int prod, std::vector<any>& rhs);
};

std::set<std::string> get_variables_used(std::string const& expr);
std::set<std::string> get_symbols_used(std::string const& expr);

Reader* new_calc_reader();

}  // end namespace MathExpr

}  // end namespace Teuchos

#endif
