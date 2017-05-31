#ifndef TEUCHOS_READER_HPP
#define TEUCHOS_READER_HPP

#include <iosfwd>
#include <functional>

#include <Teuchos_ReaderTables.hpp>
#include <Teuchos_any.hpp>

namespace Teuchos {

class Reader {
 public:
  virtual ~Reader() {}
  Reader(ReaderTablesPtr tables_in);
  void read_stream(any& result, std::istream& stream, std::string const& stream_name_in);
  void read_string(any& result, std::string const& string, std::string const& string_name);
  void read_file(any& result, std::string const& file_name);
 protected:
  /* these if there can be parse errors beyond what can
   * be specified by the ReaderTables. return true iff the operation succeeds. */
  virtual void at_shift(any& result, int token, std::string& text);
  virtual void at_reduce(any& result, int production, std::vector<any>& rhs);
 protected:
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
 protected:
  bool sensing_indent;
  std::string indent_text;
  struct IndentStackEntry {
    std::size_t line;
    std::size_t start_length;
    std::size_t end_length;
    IndentStackEntry(std::size_t l, std::size_t s, std::size_t e);
  };
  std::vector<IndentStackEntry> indent_stack;
 private:
  void at_token();
  void indent_mismatch();
  void at_token_indent();
  void at_lexer_end(std::istream& stream);
  void backtrack_to_last_accept(std::istream& stream);
  void reset_lexer_state();
  void update_position(char c);
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
