// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_PARSER_HPP
#define TEUCHOS_PARSER_HPP

/*! \file Teuchos_Parser.hpp
  \brief Declares Teuchos::Parser, ParserFail and make_lalr1_parser.
*/

#include <Teuchos_TableDecl.hpp>
#include <Teuchos_Grammar.hpp>

namespace Teuchos {

enum ActionKind {
  ACTION_NONE,
  ACTION_SHIFT,
  ACTION_REDUCE
};

struct Action {
  ActionKind kind;
  union {
    int production;
    int next_state;
  };
};

#ifdef HAVE_TEUCHOSCORE_CXX11
extern template struct Table<Action>;
#endif

struct Parser {
  GrammarPtr grammar;
  /* (state x terminal) -> action */
  Table<Action> terminal_table;
  /* (state x non-terminal) -> new state */
  Table<int> nonterminal_table;
  Parser() {}
  Parser(GrammarPtr g, int nstates_reserve);
};

int add_state(Parser& p);
int get_nstates(Parser const& p);
void add_terminal_action(Parser& p, int state, int terminal, Action action);
void add_nonterminal_action(Parser& p, int state, int nonterminal, int next_state);
Action const& get_action(Parser const& p, int state, int terminal);
int execute_action(Parser const& p, std::vector<int>& stack, Action const& action);
GrammarPtr const& get_grammar(Parser const& p);

/** \brief Tries to create LALR(1) parser tables for a given grammar.
 *
 * The exception class used throughout the TeuchosParser sub-package.
 * This exception indicates that some input did not satisfy constaints
 * on its format.
 * The most common cases are Teuchos::Reader being given a stream which
 * has some syntax error in it, making it invalid, or
 * some part of the process to convert a Teuchos::Language into
 * Teuchos::ReaderTables has similarly encountered an issue with the
 * given Teuchos::Language.
 * Such exceptions are fully recoverable and should be caught and
 * handled by users as appropriate.
 * The exception string will usually contain rich information about
 * exactly why the exception was thrown and what can be done about it.
 */
class ParserFail: public std::invalid_argument {
 public:
  ParserFail(const std::string& msg);
};

/** \brief Tries to create LALR(1) parser tables for a given grammar.
 *
 * There are several restrictions required of said grammar:
 * <ol>
 *   <li>Unambiguous: There cannot exist two ASTs which produce the same string of terminals</li>
 *   <li>LR(1): It must be possible, reading the string from left to right, to determine its
 *       AST by only looking ahead one symbol at a time while reading.</li>
 *   <li>LALR(1): More of an implementation detail, but further constrains the set of
 *                parseable grammars beyond LR(1).</li>
 * </ol>
 * It is hard to determine whether a grammar satisfies these conditions without computer assitance,
 * which is why Teuchos::make_lalr1_parser includes verbose capabilities.
 * If the grammar is not LALR(1), Teuchos::make_lalr1_parser will throw a Teuchos::ParserFail exception
 * with detailed debugging information available in Teuchos::ParserFail::what.
 * It will also write a DOT file describing the parser as a state diagram.
 * The exception string will include instruction on using
 * <a href="http://www.graphviz.org/">GraphViz</a> to convert that DOT file into
 * a format that can be viewed.
 * Together with other information in the exception string, 
 * there should be enough information to understand why the grammar is unsatisfactory
 * and hopefully devise changes to make it LALR(1).
 *
 * Designing an LALR(1) grammar is more art than science, and is probably one of the more
 * difficult aspects of using TeuchosParser.
 */
Parser make_lalr1_parser(GrammarPtr grammar, bool verbose = false);

}

#endif
