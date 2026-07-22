// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_READER_HPP
#define TEUCHOS_READER_HPP

/*! \file Teuchos_Reader.hpp
  \brief Declares Teuchos::Reader.
*/

#include <iosfwd>
#include <functional>

#include <Teuchos_ReaderTables.hpp>
#include <Teuchos_any.hpp>

namespace Teuchos {

/** \brief The main class for users to read text using TeuchosParser.
 *
 * This is the most important user-facing class in TeuchosParser.
 * A Reader represents the use of a set of ReaderTables to read (parse)
 * a given stream of text.
 * In most cases, users should create their own class deriving from Teuchos::Reader,
 * and override the Reader::at_shift and Reader::at_reduce methods to implement the mechanics of
 * their custom parser.
 *
 * To use Reader, it must be given an RCP to a set of ReaderTables.
 * In most cases, ReaderTables are constructed automatically from a Language
 * using the make_reader_tables function.
 * Notice that ReaderTables only have to be constructed once, and can be used
 * simultaneously by multiple Readers.
 * A Reader can only operate on one stream of text at a time, and contains all
 * the variable parser state.
 *
 * The Reader also provides some error reporting capability.
 * If the text stream being read is valid as defined by the Language,
 * the Teuchos::ParserFail exception will be thrown with a descriptive string.
 * It is strongly recommended that code using Reader catch this exception and
 * display its what() message to the console.
 *
 * Classes deriving from Teuchos::Reader are also strongly recommended to throw
 * the Teuchos::ParserFail exception class (or a class derived from it) to
 * indicate that the text stream being read is invalid, usually due to constraints
 * that can't be expressed in the Teuchos::Language.
 */
class Reader {
 public:

  /** \brief Virtual destructor, allows polymorphism */
  virtual ~Reader() {}

  /** \brief Constructor: accepts an RCP to ReaderTables.
   *
   * The set of streams accepted by this Reader is dictated solely by
   * the contents of the ReaderTables, and any additional checking implemented
   * in derived classes.
   * ReaderTables should be constructed from a Language using make_reader_tables.
   */
  Reader(ReaderTablesPtr tables_in);

  /** \brief The main method for reading a stream of text.
   *
   * \param result A reference to the Teuchos::any result object.
   *               If the stream is accepted, this will be the result
   *               of the Reader::at_reduce call which reduced the root nonterminal
   *               of the Language.
   * \param stream The input stream of text to be read.
   * \param name_in The name of the stream. This name is only used for
   *                error reporting.
   */
  void read_stream(any& result, std::istream& stream, std::string const& stream_name_in);

  /** \brief A convenience method for reading a string.
   *
   * This method calls Teuchos::Parser::read_stream
   */
  void read_string(any& result, std::string const& string, std::string const& string_name);

  /** \brief A convenience method for reading a file.
   *
   * This method calls Teuchos::Parser::read_stream
   */
  void read_file(any& result, std::string const& file_name);

 protected:

  /** \brief User-overridable SHIFT (token) method.
   *
   * This method will be called when the lexer (ReaderTables::lexer) accepts a token
   * in the text stream.
   * The default implementation does nothing.
   *
   * \param result A reference to the Teuchos::any result object.
   *               The contents of this object are user-defined, and will
   *               be passed later to the (rhs) argument of Reader::at_reduce
   *               when this token is reduced along with other Language symbols
   *               into a nonterminal.
   * \param token  An integer index denoting which token was accepted.
   *               This is the index of the token in the Language::tokens vector.
   * \param text   The full text which matched the token's regular expression.
   *               This is a non-const reference so that users may do things
   *               like move or swap it directly into the result object.
   *               Teuchos::Reader will find the longest substring accepted by
   *               the Teuchos::ReaderTables::lexer starting from its current
   *               stream position, then repeat from the new stream position.
   */
  virtual void at_shift(any& result, int token, std::string& text);

  /** \brief User-overridable REDUCE (production) method
   *
   * This method will be called when the lexer (ReaderTables::lexer) accepts a token
   * in the text stream.
   * The default implementation does nothing.
   *
   * \param result A reference to the Teuchos::any result object.
   *               The contents of this object are user-defined, and mayb
   *               be passed later to the (rhs) argument of Reader::at_reduce
   *               when this nonterminal is reduced along with other Language symbols
   *               into a nonterminal, of it will be assigned to the result
   *               value of Teuchos::Reader::read_stream itself iff the LHS
   *               of (production) is the root nonterminal.
   * \param token  An integer index denoting which production was reduced.
   *               This is the index of the production in the Language::productions vector.
   * \param text   A vector containing the result objects associated with each
   *               RHS symbol in the production.
   *               These objects are the result of either previous calls to Reader::at_shift
   *               (if their associated symbol is a token) 
   *               or Reader::at_reduce (if the associated symbol is a nonterminal).
   *               This is a non-const reference so that users may do things
   *               like move or swap it (or its elements) directly into the result object.
   */
  virtual void at_reduce(any& result, int production, std::vector<any>& rhs);

 protected: // variables for normal Context Free Grammar parsing and error reporting

  ReaderTablesPtr tables;
  Parser const& parser;
  FiniteAutomaton const& lexer;
  GrammarPtr grammar;
  std::size_t line;
  std::size_t column;
  int lexer_state;
  std::string lexer_text;
  std::string line_text;
  int lexer_token;
  std::size_t last_lexer_accept;
  std::size_t last_lexer_accept_line;
  std::size_t last_lexer_accept_column;
  std::string last_lexer_accept_line_text;
  int parser_state;
  std::vector<int> parser_stack;
  std::vector<any> value_stack;
  std::vector<any> reduction_rhs;
  std::string stream_name;
  bool did_accept;

 protected: // variables for indentation-sensitive language parsing

  bool sensing_indent;
  std::string indent_text;
  struct IndentStackEntry {
    std::size_t line;
    std::size_t start_length;
    std::size_t end_length;
    IndentStackEntry(std::size_t l, std::size_t s, std::size_t e);
  };
  // this is the stack that shows, for the current leading indentation
  // characters, which subset of them came from each nested increase
  // in indentation
  std::vector<IndentStackEntry> indent_stack;
  // this stack notes, for each symbol in the pushdown automaton
  // stack, how many characters indent the line that that symbol
  // starts on
  std::vector<std::size_t> symbol_indentation_stack;

 private: // helper methods

  void at_token(std::istream& stream);
  void indent_mismatch();
  void at_token_indent(std::istream& stream);
  void at_lexer_end(std::istream& stream);
  void backtrack_to_last_accept(std::istream& stream);
  void reset_lexer_state();
  void update_position(char c);
  void error_print_line(std::istream& is, std::ostream& os);
};

class DebugReader : public Reader {
 public:
  DebugReader(ReaderTablesPtr tables_in, std::ostream& os_in);
  virtual ~DebugReader() {}
 protected:
  virtual void at_shift(any& result, int token, std::string& text);
  virtual void at_reduce(any& result, int token, std::vector<any>& rhs);
 private:
  std::ostream& os;
};

}

#endif
