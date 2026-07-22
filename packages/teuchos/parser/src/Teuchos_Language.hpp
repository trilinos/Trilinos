// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_LANGUAGE_HPP
#define TEUCHOS_LANGUAGE_HPP

/*! \file Teuchos_Language.hpp
  \brief Declares Teuchos::Language.
*/

#include <string>
#include <vector>
#include <iosfwd>

#include <Teuchos_Grammar.hpp>
#include <Teuchos_FiniteAutomaton.hpp>
#include <Teuchos_ReaderTables.hpp>
#include <Teuchos_RCP.hpp>

namespace Teuchos {


/** \brief The main class for users to define a language using TeuchosParser.
 *
 * This is the second most important user-facing class in TeuchosParser,
 * the first being Teuchos::Reader.
 * With Language, users can define their own language at runtime,
 * which can then be converted into ReaderTables by make_reader_tables
 * and used by Reader to read text streams.
 *
 * TeuchosParser is heavily inspired by the <tt>lex</tt> and <tt>yacc</tt>
 * UNIX programs, whose modern implementations are
 * <a href="https://github.com/westes/flex">Flex</a> and
 * <a href="https://www.gnu.org/software/bison/">Bison</a>.
 * These are programs which read a custom file format and output C or C++ source files.
 * At the time of this writing, the C++ support (particularly in Flex) leaves
 * something to be desired, which is part of the reason for creating TeuchosParser.
 * TeuchosParser does a subset of what Flex and Bison can do, but it accepts
 * a Language object as input and outputs ReaderTables, which can later
 * be used by a Reader.
 * All of these are in-memory, pure C++ constructs.
 * TeuchosParser supports LALR(1) grammars for the parser and
 * basic regular expressions in the lexer.
 *
 * A Language has two portions: the productions and the tokens.
 *
 * The Language::productions vector encodes a Context Free Grammar as a set
 * of Teuchos::Language::Production objects.
 * A Context Free Grammar consists of a set of terminal symbols (referred
 * to here as tokens), a set of nonterminal symbols, and a set of productions.
 * The productions consist of a single symbol on the left hand side, and a
 * string of symbols on the right hand side.
 *
 * A production means that the symbol on the left hand side may be substituted
 * for the string on the right hand side, or vice versa.
 * The grammar also has a root nonterminal symbol.
 * Every string acceptable by the grammar can be formed by:
 * <ol>
 *   <li>Start with a string containing only the root nonterminal</li>
 *   <li>Choose a nonterminal in the string</li>
 *   <li>Choose a production whose LHS is the chosen nonterminal</li>
 *   <li>Substitute the nonterminal with the RHS of the production</li>
 *   <li>Repeat 2-4 until there are only terminal symbols in the string</li>
 * </ol>
 * This is the top-down perspective of forming a string of terminals based
 * on the choices of substitutions made.
 * Parsing is the bottom-up equivalent, of taking a string of terminals
 * and trying to deduce what substitutions were made to form it.
 * One can define an Abstract Syntax Tree (AST), as a tree defined by
 * the substitutions made.
 * Each tree node contains a symbol, the root node contains the root nonterminal,
 * and each non-leaf tree node must contain the LHS of a production while its
 * children must contain the symbols in the RHS of that production.
 * The leaf nodes of the AST must contain terminal symbols, and reading the leaf
 * nodes of the tree from left to right spells out the string that is accepted.
 *
 * Please see the Teuchos::make_lalr1_parser documentation for information
 * on the contraints required of a grammar.
 *
 * The Language::tokens vector defines each token as a regular expression.
 * The tokens and their regular expressions together form a lexer, which
 * can be a useful pre-processing step to turn raw text into high-level
 * tokens for consumption by the parser.
 * Please see the Teuchos::make_lexer, Teuchos::make_dfa, and Teuchos::Reader documentation
 * for information about requirements and the definition of the lexer.
 *
 * All symbols in Teuchos::Language are denoted by strings.
 * There are no restrictions on what strings can be used, their contents
 * have no special interpretation.
 * Thus it is wise to choose symbol names which are as convenient to read
 * as possible.
 *
 * The idiomatic way to construct a Language object is as follows:
 * <ol>
 *   <li>Create an enumeration for your productions</li>
 *   <li>Create an enumeration for your tokens</li>
 *   <li>Define each production using the overloaded operators provided:
 *     \code
 *     lang.productions[PROD_PARENS]("expr1") >> "(", "expr2", ")";
 *     \endcode
 *   </li>
 *   <li>Define each token using the overloaded operator():
 *     \code
 *     lang.tokens[TOK_LPAREN]("(", "\\)");
 *     \endcode
 *   </li>
 * </ol>
 * If your language doesn't change, it can be good design to store it as a singleton
 * LanguagePtr object.
 *
 * Please see Teuchos_XML.cpp, Teuchos_YAML.cpp, and Calc.cpp for examples
 * of 
 */
struct Language {
  struct Token {
    std::string name;
    std::string regex;
    void operator()(std::string const& name_in, std::string const& regex_in);
  };
  typedef std::vector<Token> Tokens;
  /** \brief vector of tokens */
  Tokens tokens;
  typedef std::vector<std::string> RHS;
  struct Production;
  struct RHSBuilder {
    Production& prod;
    RHSBuilder(Production& prod_in);
    RHSBuilder& operator,(std::string const& rhs_item);
    RHSBuilder& operator>>(std::string const& rhs_item);
  };
  struct Production {
    std::string lhs;
    RHS rhs;
    RHSBuilder operator()(std::string const& lhs_in);
  };
  typedef std::vector<Production> Productions;
  /** \brief vector of productions */
  Productions productions;
};

/** \brief an RCP to a const Language */
typedef RCP<const Language> LanguagePtr;

GrammarPtr make_grammar(Language const& language);

/** \brief construct a lexer for the Language tokens.
 *
 * First, a lexer for each token will be constructed using Teuchos::make_dfa.
 * Then the tokens will be united.
 * In any case when one string would match two or different more tokens,
 * only the token with the lowest index is matched (i.e. the one that
 * is listed first in Language::tokens.
 * Note that Teuchos::Reader will accept the longest match from a given
 * starting point.
 */
void make_lexer(FiniteAutomaton& result, Language const& language);


/** \brief constructs ReaderTables for the given Language.
 *
 * See Teuchos::make_lexer and Teuchos::make_lalr1_parser.
 */
ReaderTablesPtr make_reader_tables(Language const& language);

std::ostream& operator<<(std::ostream& os, Language const& lang);

}

#endif
