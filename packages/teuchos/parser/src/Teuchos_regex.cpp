#include <Teuchos_regex.hpp>

#include <cassert>
#include <iostream>
#include <sstream>

#include <Teuchos_vector.hpp>
#include <Teuchos_string.hpp>
#include <Teuchos_chartab.hpp>
#include <Teuchos_build_parser.hpp>
#include <Teuchos_set.hpp>
#include <Teuchos_reader.hpp>
#include <Teuchos_chartab.hpp>

namespace Teuchos {
namespace regex {

Language build_language() {
  /* The top produtions were from the "grep.y" YACC grammar in the source
     code for Plan 9's grep utility, see here:
https://github.com/wangeguo/plan9/blob/master/sys/src/cmd/grep/grep.y
     The "set" related productions
     are from a grammar intended to be used by ProLog to parse Perl's regular
     expressions, see here:
http://www.cs.sfu.ca/~cameron/Teaching/384/99-3/regexp-plg.html */
  Language out;
  auto& prods = out.productions;
  prods.resize(NPRODS);
  prods[PROD_REGEX] = {"regex", {"union"}};
  prods[PROD_UNION_DECAY] = {"union", {"concat"}};
  prods[PROD_UNION] = {"union", {"union", "|", "concat"}}; // union
  prods[PROD_CONCAT_DECAY] = {"concat", {"qualified"}};
  prods[PROD_CONCAT] = {"concat", {"concat", "qualified"}}; // concatenation
  prods[PROD_QUAL_DECAY] = {"qualified", {"single"}};
  prods[PROD_STAR] = {"qualified", {"qualified", "*"}};
  prods[PROD_PLUS] = {"qualified", {"qualified", "+"}};
  prods[PROD_MAYBE] = {"qualified", {"qualified", "?"}};
  prods[PROD_SINGLE_CHAR] = {"single", {"char"}};
  prods[PROD_ANY] = {"single", {"."}}; // any
  prods[PROD_SINGLE_SET] = {"single", {"set"}};
  prods[PROD_PARENS_UNION] = {"single", {"(", "union", ")"}};
  prods[PROD_SET_POSITIVE] = {"set", {"positive-set"}};
  prods[PROD_SET_NEGATIVE] = {"set", {"negative-set"}};
  prods[PROD_POSITIVE_SET] = {"positive-set", {"[", "set-items", "]"}};
  prods[PROD_NEGATIVE_SET] = {"negative-set", {"[", "^", "set-items", "]"}};
  prods[PROD_SET_ITEMS_DECAY] = {"set-items", {"set-item"}};
  prods[PROD_SET_ITEMS_ADD] = {"set-items", {"set-items", "set-item"}};
  prods[PROD_SET_ITEM_CHAR] = {"set-item", {"char"}};
  prods[PROD_SET_ITEM_RANGE] = {"set-item", {"range"}};
  prods[PROD_RANGE] = {"range", {"char", "-", "char"}};
  out.tokens.resize(NTOKS);
  /* either one of the non-meta characters, or anything preceded by the escape slash */
  out.tokens[TOK_CHAR] = {"char", "[^\\\\\\.\\[\\]\\(\\)\\|\\-\\^\\*\\+\\?]|\\\\."};
  out.tokens[TOK_DOT] = {".", "\\."};
  out.tokens[TOK_LRANGE] = {"[", "\\]"};
  out.tokens[TOK_RRANGE] = {"]", "\\]"};
  out.tokens[TOK_LPAREN] = {"(", "\\("};
  out.tokens[TOK_RPAREN] = {")", "\\)"};
  out.tokens[TOK_UNION] = {"|", "\\|"};
  out.tokens[TOK_RANGE] = {"-", "\\-"};
  out.tokens[TOK_NEGATE] = {"^", "\\^"};
  out.tokens[TOK_STAR] = {"*", "\\*"};
  out.tokens[TOK_PLUS] = {"+", "\\+"};
  out.tokens[TOK_MAYBE] = {"?", "\\?"};
  return out;
}

/* bootstrap ! This lexer is used to build the ReaderTables that read
   regular expressions themselves, so it can't depend on that reader ! */
FiniteAutomaton build_lexer() {
  std::string meta_chars_str = ".[]()|-^*+?";
  std::set<int> all_chars;
  for (int i = 0; i < NCHARS; ++i) all_chars.insert(i);
  auto nonmeta_chars = all_chars;
  for (auto meta_char : meta_chars_str) {
    auto it = nonmeta_chars.find(get_symbol(meta_char));
    nonmeta_chars.erase(it);
  }
  auto lex_nonmeta = make_set_nfa(NCHARS, nonmeta_chars, TOK_CHAR);
  auto lex_slash = make_char_single_nfa('\\');
  auto lex_any = make_set_nfa(NCHARS, all_chars);
  auto lex_escaped = concat(lex_slash, lex_any, TOK_CHAR);
  auto lex_char = unite(lex_nonmeta, lex_escaped);
  FiniteAutomaton lex_metachars;
  for (int i = 0; i < size(meta_chars_str); ++i) {
    int token = TOK_CHAR + i + 1;
    auto lex_metachar = make_char_single_nfa(at(meta_chars_str, i), token);
    if (i) lex_metachars = unite(lex_metachars, lex_metachar);
    else lex_metachars = lex_metachar;
  }
  auto out = unite(lex_char, lex_metachars);
  return simplify(make_deterministic(out));
}

ReaderTablesPtr ask_reader_tables() {
  static ReaderTablesPtr ptr;
  if (ptr.use_count() == 0) {
    auto lang = regex::ask_language();
    auto grammar = build_grammar(*lang);
    auto parser = accept_parser(build_lalr1_parser(grammar));
    auto lexer = regex::build_lexer();
    IndentInfo indent_info;
    indent_info.is_sensitive = false;
    indent_info.indent_token = -1;
    indent_info.dedent_token = -1;
    ptr.reset(new ReaderTables{parser, lexer, indent_info});
  }
  return ptr;
}

LanguagePtr ask_language() {
  static LanguagePtr ptr;
  if (ptr.use_count() == 0) {
    ptr.reset(new Language(build_language()));
  }
  return ptr;
}

FiniteAutomaton build_dfa(std::string const& name, std::string const& regex, int token) {
  /* special "regex"es for indentation sensitivity support */
  if (regex == "]INDENT[" || regex == "]DEDENT[" || regex == "]EQDENT[" || regex == "]NODENT[") {
    return build_dfa(name, "\r?\n[ \t]*", token);
  }
  auto reader = regex::Reader(token);
  bool ok = reader.read_string(name, regex);
  assert(ok);
  return reader.move_result<FiniteAutomaton>();
}

regex::Reader::Reader(int result_token_in):
  Teuchos::Reader(regex::ask_reader_tables()),
  result_token(result_token_in) {
}

any regex::Reader::at_shift(int token, std::string& text) {
  if (token != TOK_CHAR) { return any(); }
  if (size(text) == 1) { return any(text[0]); }
  else if (size(text) == 2) {
    assert(text[0] == '\\');
    return any(text[1]);
  } else {
    std::cerr << "BUG: regex char text is \"" << text << "\"\n";
    abort();
  }
}

any regex::Reader::at_reduce(int production, std::vector<any>& rhs) {
  switch (production) {
    case PROD_REGEX:
      return simplify(make_deterministic(
            move_value<FiniteAutomaton>(at(rhs, 0))));
    case PROD_UNION_DECAY:
    case PROD_CONCAT_DECAY:
    case PROD_QUAL_DECAY:
    case PROD_SET_ITEMS_DECAY:
    case PROD_SET_ITEM_RANGE:
      return at(rhs, 0);
    case PROD_UNION:
      return unite(
          move_value<FiniteAutomaton>(at(rhs, 0)),
          move_value<FiniteAutomaton>(at(rhs, 2)));
    case PROD_CONCAT:
      {
        auto& a_any = at(rhs, 0);
        auto& b_any = at(rhs, 1);
        auto a = move_value<FiniteAutomaton>(a_any);
        auto b = move_value<FiniteAutomaton>(b_any);
        return concat(a, b, result_token);
      }
    case PROD_STAR:
      return star(
          move_value<FiniteAutomaton>(at(rhs, 0)), result_token);
    case PROD_PLUS:
      return plus(
          move_value<FiniteAutomaton>(at(rhs, 0)), result_token);
    case PROD_MAYBE:
      return maybe(
          move_value<FiniteAutomaton>(at(rhs, 0)), result_token);
    case PROD_SINGLE_CHAR:
      return make_char_single_nfa(
          move_value<char>(at(rhs, 0)), result_token);
    case PROD_ANY:
      return make_range_nfa(NCHARS, 0, NCHARS - 1, result_token);
    case PROD_SINGLE_SET:
      return make_char_set_nfa(
          move_value<std::set<char>>(at(rhs, 0)), result_token);
    case PROD_PARENS_UNION:
      return at(rhs, 1);
    case PROD_SET_POSITIVE:
      return at(rhs, 0);
    case PROD_SET_NEGATIVE:
      return negate_set(move_value<std::set<char>>(at(rhs, 0)));
    case PROD_POSITIVE_SET:
      return at(rhs, 1);
    case PROD_NEGATIVE_SET:
      return at(rhs, 2);
    case PROD_SET_ITEMS_ADD:
      return unite(
          move_value<std::set<char>>(at(rhs, 0)),
          move_value<std::set<char>>(at(rhs, 1)));
    case PROD_SET_ITEM_CHAR:
      return std::set<char>({move_value<char>(at(rhs, 0))});
    case PROD_RANGE:
      {
        std::set<char> set;
        for (char c = move_value<char>(at(rhs, 0));
             c <= move_value<char>(at(rhs, 2)); ++c) {
          set.insert(c);
        }
        return set;
      }
  }
  std::cerr << "BUG: unexpected production " << production << '\n';
  abort();
}

}  // end namespace regex
}  // end namespace Teuchos
