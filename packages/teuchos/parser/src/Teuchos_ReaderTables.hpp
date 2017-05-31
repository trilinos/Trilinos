#ifndef TEUCHOS_READER_TABLES_HPP
#define TEUCHOS_READER_TABLES_HPP

#include <Teuchos_FiniteAutomaton.hpp>
#include <Teuchos_Parser.hpp>
#include <Teuchos_RCP.hpp>

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

typedef RCP<const ReaderTables> ReaderTablesPtr;

}

#endif
