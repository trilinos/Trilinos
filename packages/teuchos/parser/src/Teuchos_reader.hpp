#ifndef TEUCHOS_READER_HPP
#define TEUCHOS_READER_HPP

#include <iosfwd>
#include <functional>

#include <Teuchos_reader_tables.hpp>
#include <any/Teuchos_any.hpp>

namespace Teuchos {

class Reader {
 public:
  Reader() = default;
  Reader(Reader const& other) = default;
  virtual ~Reader() = default;
  Reader(ReaderTablesPtr tables_in);
  bool read_stream(std::string const& stream_name_in, std::istream& stream);
  bool read_string(std::string const& string_name, std::string const& string);
  bool read_file(std::string const& file_name);
  template <typename T>
  T&& move_result() {
    assert(size(value_stack) == 1);
    return move_value<T>(value_stack.back());
  }
 protected:
  /* override these if there can be parse errors beyond what can
   * be specified by the ReaderTables. return true iff the operation succeeds. */
  virtual bool at_shift(int token, std::string& text, any& result);
  virtual bool at_reduce(int token, std::vector<any>& rhs, any& result);
  /* override these if parsing cannot fail as long as the text matches
   * the specification in the ReaderTables. (simple languages) */
  virtual any at_shift(int token, std::string& text);
  virtual any at_reduce(int token, std::vector<any>& rhs);
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
  };
  std::vector<IndentStackEntry> indent_stack;
 private:
  bool at_token();
  void indent_mismatch();
  bool at_token_indent();
  bool at_lexer_end(std::istream& stream);
  void backtrack_to_last_accept(std::istream& stream);
  void reset_lexer_state();
  void update_position(char c);
};

class DebugReader : public Reader {
 public:
  DebugReader(ReaderTablesPtr tables_in);
  DebugReader(DebugReader const& other) = default;
  virtual ~DebugReader() override = default;
 protected:
  virtual any at_shift(int token, std::string& text) override;
  virtual any at_reduce(int token, std::vector<any>& rhs) override;
};

}

#endif
