// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_READER_TABLES_HPP
#define TEUCHOS_READER_TABLES_HPP

/*! \file Teuchos_ReaderTables.hpp
  \brief Declares Teuchos::ReaderTables.
*/

#include <Teuchos_FiniteAutomaton.hpp>
#include <Teuchos_Parser.hpp>
#include <Teuchos_RCP.hpp>

namespace Teuchos {

struct IndentInfo {
  bool is_sensitive;
  int indent_token;
  int dedent_token;
  int newline_token;
};

/** \brief Parser and lexer tables specifying how to read a Language. */
struct ReaderTables {
  /** \brief parser. */
  Parser parser;
  /** \brief lexer. */
  FiniteAutomaton lexer;
  IndentInfo indent_info;
};

/** \brief an RCP to a const ReaderTables */
typedef RCP<const ReaderTables> ReaderTablesPtr;

}

#endif
