#ifndef TEUCHOS_READER_TABLES_HPP
#define TEUCHOS_READER_TABLES_HPP

#include <memory>

#include <Teuchos_finite_automaton.hpp>
#include <Teuchos_parser.hpp>

namespace Teuchos {

struct IndentInfo {
  bool is_sensitive;
  int indent_token;
  int dedent_token;
  int eqdent_token;
  int nodent_token;
};

struct ReaderTables {
  Parser parser;
  FiniteAutomaton lexer;
  IndentInfo indent_info;
};

using ReaderTablesPtr = std::shared_ptr<ReaderTables const>;

}

#endif
