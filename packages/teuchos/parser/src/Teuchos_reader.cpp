#include "Teuchos_reader.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <ios>
#include <cstdlib>
#include <set>

#include "Teuchos_string.hpp"

namespace Teuchos {

namespace {

void print_indicator(std::ostream& os, std::string const& above, std::size_t pos) {
  for (std::size_t i = 0; i < pos; ++i) {
    if (above.at(i) == '\t') os << '\t';
    else os << ' ';
  }
  os << "^\n";
}

void print_underline(std::ostream& os, std::string const& above, std::size_t start, std::size_t end) {
  for (std::size_t i = 0; i < start; ++i) {
    if (above.at(i) == '\t') os << '\t';
    else os << ' ';
  }
  for (auto i = start; i < end; ++i) os << '~';
  os << '\n';
}

} // end anonymous namespace

bool Reader::at_token() {
  bool done = false;
  /* this can loop arbitrarily as reductions are made,
     because they don't consume the token */
  while (!done) {
    auto parser_action = get_action(parser, parser_state, lexer_token);
    if (parser_action.kind == ACTION_NONE) {
      std::cerr << "error: Parser failure at line " << line;
      std::cerr << " column " << column << " of " << stream_name << '\n';
      std::cerr << line_text << '\n';
      print_indicator(std::cerr, line_text, line_text.size());
      std::set<std::string> expect_names;
      for (int expect_token = 0;
           expect_token < grammar->nterminals; ++expect_token) {
        auto expect_action = get_action(parser, parser_state, expect_token);
        if (expect_action.kind != ACTION_NONE) {
          expect_names.insert(at(grammar->symbol_names, expect_token));
        }
      }
      std::cerr << "Expected one of {";
      for (auto it = expect_names.begin(); it != expect_names.end(); ++it) {
        if (it != expect_names.begin()) std::cerr << ", ";
        if (*it == ",") std::cerr << "','";
        else std::cerr << *it;
      }
      std::cerr << "}\n";
      std::cerr << "Got: " << at(grammar->symbol_names, lexer_token) << '\n';
      std::cerr << "Parser was in state " << parser_state << '\n';
      return false;
    } else if (parser_action.kind == ACTION_SHIFT) {
      Teuchos::any shift_result;
      bool shift_ok = this->at_shift(lexer_token, lexer_text, shift_result);
      if (!shift_ok) return false;
      value_stack.emplace_back(std::move(shift_result));
      done = true;
    } else if (parser_action.kind == ACTION_REDUCE) {
      if (parser_action.production == get_accept_production(*grammar)) {
        did_accept = true;
        return true;
      }
      auto& prod = at(grammar->productions, parser_action.production);
      reduction_rhs.clear();
      for (int i = 0; i < size(prod.rhs); ++i) {
        reduction_rhs.emplace_back(std::move(
              at(value_stack, size(value_stack) - size(prod.rhs) + i)));
      }
      resize(value_stack, size(value_stack) - size(prod.rhs));
      Teuchos::any reduce_result;
      bool reduce_ok = this->at_reduce(parser_action.production, reduction_rhs, reduce_result);
      if (!reduce_ok) return false;
      value_stack.emplace_back(std::move(reduce_result));
    } else {
      std::cerr << "SERIOUS BUG: Action::kind enum value not in range\n";
      abort();
    }
    parser_state = execute_action(parser, parser_stack, parser_action);
  }
  return true;
}

void Reader::indent_mismatch() {
  assert(!indent_stack.empty());
  auto top = indent_stack.back();
  std::cerr << "error: Indentation characters beginning line " << line << " of " << stream_name
    << " don't match those beginning line " << top.line << '\n';
  std::cerr << "It is strongly recommended not to mix tabs and spaces in indentation-sensitive formats\n";
}

bool Reader::at_token_indent() {
  if (!sensing_indent || lexer_token != tables->indent_info.nodent_token) {
    return at_token();
  }
  assert(at(lexer_text, 0) == '\n');
  auto lexer_indent = lexer_text.substr(1, std::string::npos);
  auto minlen = std::min(lexer_indent.length(), indent_text.length());
  if (lexer_indent.length() > indent_text.length()) {
    if (0 != lexer_indent.compare(0, indent_text.length(), indent_text)) {
      indent_mismatch();
      return false;
    }
    indent_stack.push_back({line, indent_text.length(), lexer_indent.length()});
    indent_text = lexer_indent;
    lexer_token = tables->indent_info.indent_token;
    return at_token();
  } else if (lexer_indent.length() < indent_text.length()) {
    if (0 != indent_text.compare(0, lexer_indent.length(), lexer_indent)) {
      indent_mismatch();
      return false;
    }
    bool first = true;
    while (!indent_stack.empty()) {
      auto top = indent_stack.back();
      if (top.end_length <= minlen) break;
      indent_stack.pop_back();
      lexer_token = tables->indent_info.dedent_token;
      if (!at_token()) return false;
      if (first) {
        lexer_text.clear();
        first = false;
      }
    }
    if (first) lexer_text.clear();
    indent_text = lexer_indent;
    return true;
  } else {
    if (0 != lexer_indent.compare(indent_text)) {
      indent_mismatch();
      return false;
    }
    lexer_token = tables->indent_info.eqdent_token;
    return at_token();
  }
}

void Reader::backtrack_to_last_accept(std::istream& stream) {
  /* all the last_accept and backtracking is driven by
    the "accept the longest match" rule */
  line = last_lexer_accept_line;
  column = last_lexer_accept_column;
  line_text = last_lexer_accept_line_text;
  while (lexer_text.size() > last_lexer_accept) {
    bool ok = !stream.unget().fail();
    assert(ok);
    lexer_text.pop_back();
  }
}

void Reader::reset_lexer_state() {
  lexer_state = 0;
  lexer_text.clear();
  lexer_token = -1;
}

bool Reader::at_lexer_end(std::istream& stream) {
  if (lexer_token == -1) {
    if (lexer_text.find('\n') == std::string::npos) {
      std::cerr << "error: Could not tokenize this (line " <<  line;
      std::cerr << " column " << column << " of " << stream_name << "):\n";
      std::cerr << line_text << '\n';
      assert(line_text.size() >= lexer_text.size());
      print_underline(std::cerr, line_text, line_text.size() - lexer_text.size(), line_text.size());
    } else {
      std::cerr << "error: Could not tokenize this (ends at line " << line;
      std::cerr << " column " << column << " of " << stream_name << "):\n";
      std::cerr << lexer_text << '\n';
    }
    return false;
  }
  backtrack_to_last_accept(stream);
  if (!at_token_indent()) { return false; }
  reset_lexer_state();
  return true;
}

Reader::Reader(ReaderTablesPtr tables_in):
  tables(tables_in),
  parser(tables->parser),
  lexer(tables->lexer),
  grammar(get_grammar(parser))
{
  assert(get_determinism(lexer));
}

void Reader::update_position(char c) {
  if (c == '\n') {
    ++line;
    column = 1;
    line_text.clear();
  } else {
    ++column;
  }
}

bool Reader::read_stream(std::string const& stream_name_in, std::istream& stream) {
  line = 1;
  column = 1;
  lexer_state = 0;
  lexer_text.clear();
  line_text.clear();
  lexer_token = -1;
  parser_state = 0;
  parser_stack.clear();
  parser_stack.push_back(parser_state);
  value_stack.clear();
  did_accept = false;
  stream_name = stream_name_in;
  if (tables->indent_info.is_sensitive) {
    sensing_indent = true;
    indent_text.clear();
    indent_stack.clear();
    /* pretend the stream starts with a newline so we can
       detect an INDENT on the first line. don't update the
       line/column pointers though. */
    char c = '\n';
    lexer_text.push_back(c);
    auto lexer_symbol = get_symbol(c);
    lexer_state = step(lexer, lexer_state, lexer_symbol);
    assert(lexer_state != -1);
    lexer_token = accepts(lexer, lexer_state);
    assert(lexer_token == tables->indent_info.nodent_token);
    last_lexer_accept = lexer_text.size();
    last_lexer_accept_line = 1;
    last_lexer_accept_column = 1;
    last_lexer_accept_line_text = line_text;
  } else {
    sensing_indent = false;
  }
  char c;
  while (stream.get(c)) {
    if (!is_symbol(c)) {
      std::cerr << "error: Unexpected character code " << int(c);
      std::cerr << " at line " << line << " column " << column;
      std::cerr << " of " << stream_name << '\n';
      if (!line_text.empty()) {
        std::cerr << line_text << '\n';
        print_indicator(std::cerr, line_text, line_text.size());
      }
      return false;
    }
    line_text.push_back(c);
    lexer_text.push_back(c);
    auto lexer_symbol = get_symbol(c);
    lexer_state = step(lexer, lexer_state, lexer_symbol);
    if (lexer_state == -1) {
      if (!at_lexer_end(stream)) return false;
    } else {
      auto token = accepts(lexer, lexer_state);
      update_position(c);
      if (token != -1) {
        lexer_token = token;
        last_lexer_accept = lexer_text.size();
        last_lexer_accept_line = line;
        last_lexer_accept_column = column;
        last_lexer_accept_line_text = line_text;
      }
    }
  }
  if (last_lexer_accept < lexer_text.size()) {
    auto bad_str = lexer_text.substr(last_lexer_accept, std::string::npos);
    std::cerr << "error: Could not tokenize \"" << bad_str;
    std::cerr << "\" at end of " << stream_name << '\n';
    return false;
  }
  if (!at_lexer_end(stream)) return false;
  lexer_token = get_end_terminal(*grammar);
  if (!at_token()) return false;
  return did_accept;
}

bool Reader::read_string(std::string const& string_name, std::string const& string) {
  std::istringstream stream(string);
  return read_stream(string_name, stream);
}

bool Reader::read_file(std::string const& file_name) {
  std::ifstream stream(file_name.c_str());
  return read_stream(file_name, stream);
}

bool Reader::at_shift(int token, std::string& text, any& result) {
  result = this->at_shift(token, text);
  return true;
}

bool Reader::at_reduce(int prod, std::vector<any>& rhs, any& result) {
  result = this->at_reduce(prod, rhs);
  return true;
}

any Reader::at_shift(int, std::string&) {
  return any();
}

any Reader::at_reduce(int, std::vector<any>&) {
  return any();
}

DebugReader::DebugReader(ReaderTablesPtr tables_in):Reader(tables_in) {}

any DebugReader::at_shift(int token, std::string& text) {
  std::string text_escaped;
  for (auto c : text) {
    switch (c) {
      case '\n': text_escaped.append("\\n"); break;
      case '\t': text_escaped.append("\\t"); break;
      case '\r': text_escaped.append("\\r"); break;
      default: text_escaped.push_back(c);
    }
  }
  std::cerr << "SHIFT (" << at(grammar->symbol_names, token) << ")[" << text_escaped << "]\n";
  return text_escaped;
}

any DebugReader::at_reduce(int prod_i, std::vector<any>& rhs) {
  std::cerr << "REDUCE";
  std::string lhs_text;
  auto& prod = at(grammar->productions, prod_i);
  for (int i = 0; i < size(prod.rhs); ++i) {
    auto& rhs_name = at(grammar->symbol_names, at(prod.rhs, i));
    auto rhs_text = move_value<std::string>(at(rhs, i));
    std::cerr << " (" << rhs_name << ")[" << rhs_text << "]";
    lhs_text.append(rhs_text);
  }
  auto& lhs_name = at(grammar->symbol_names, prod.lhs);
  std::cerr << " -> (" << lhs_name << ")[" << lhs_text << "]\n";
  return lhs_text;
}

}  // end namespace Teuchos
