// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_set.hpp"

#include "Teuchos_regex.hpp"

#include <iostream>
#include <sstream>

#include "Teuchos_Assert.hpp"
#include "Teuchos_Parser.hpp"
#include "Teuchos_vector.hpp"
#include "Teuchos_string.hpp"
#include "Teuchos_chartab.hpp"
#include "Teuchos_Reader.hpp"
#include "Teuchos_chartab.hpp"

namespace Teuchos {
namespace regex {

Language make_language() {
  /* The top produtions were from the "grep.y" YACC grammar in the source
     code for Plan 9's grep utility, see here:
https://github.com/wangeguo/plan9/blob/master/sys/src/cmd/grep/grep.y
     The "set" related productions
     are from a grammar intended to be used by ProLog to parse Perl's regular
     expressions, see here:
http://www.cs.sfu.ca/~cameron/Teaching/384/99-3/regexp-plg.html */
  Language out;
  Language::Productions& prods = out.productions;
  prods.resize(NPRODS);
  prods[PROD_REGEX]("regex") >> "union";
  prods[PROD_UNION_DECAY]("union") >> "concat";
  prods[PROD_UNION]("union") >> "union", "|", "concat"; // union
  prods[PROD_CONCAT_DECAY]("concat") >> "qualified";
  prods[PROD_CONCAT]("concat") >> "concat", "qualified"; // concatenation
  prods[PROD_QUAL_DECAY]("qualified") >> "single";
  prods[PROD_STAR]("qualified") >> "qualified", "*";
  prods[PROD_PLUS]("qualified") >> "qualified", "+";
  prods[PROD_MAYBE]("qualified") >> "qualified", "?";
  prods[PROD_SINGLE_CHAR]("single") >> "char";
  prods[PROD_ANY]("single") >> "."; // any
  prods[PROD_SINGLE_SET]("single") >> "set";
  prods[PROD_PARENS_UNION]("single") >> "(", "union", ")";
  prods[PROD_SET_POSITIVE]("set") >> "positive-set";
  prods[PROD_SET_NEGATIVE]("set") >> "negative-set";
  prods[PROD_POSITIVE_SET]("positive-set") >> "[", "set-items", "]";
  prods[PROD_NEGATIVE_SET]("negative-set") >> "[", "^", "set-items", "]";
  prods[PROD_SET_ITEMS_DECAY]("set-items") >> "set-item";
  prods[PROD_SET_ITEMS_ADD]("set-items") >> "set-items", "set-item";
  prods[PROD_SET_ITEM_CHAR]("set-item") >> "char";
  prods[PROD_SET_ITEM_RANGE]("set-item") >> "range";
  prods[PROD_RANGE]("range") >> "char", "-", "char";
  out.tokens.resize(NTOKS);
  /* either one of the non-meta characters, or anything preceded by the escape slash */
  out.tokens[TOK_CHAR]("char", "[^\\\\\\.\\[\\]\\(\\)\\|\\-\\^\\*\\+\\?]|\\\\.");
  out.tokens[TOK_DOT](".", "\\.");
  out.tokens[TOK_LRANGE]("[", "\\]");
  out.tokens[TOK_RRANGE]("]", "\\]");
  out.tokens[TOK_LPAREN]("(", "\\(");
  out.tokens[TOK_RPAREN](")", "\\)");
  out.tokens[TOK_UNION]("|", "\\|");
  out.tokens[TOK_RANGE]("-", "\\-");
  out.tokens[TOK_NEGATE]("^", "\\^");
  out.tokens[TOK_STAR]("*", "\\*");
  out.tokens[TOK_PLUS]("+", "\\+");
  out.tokens[TOK_MAYBE]("?", "\\?");
  return out;
}

/* bootstrap ! This lexer is used to build the ReaderTables that read
   regular expressions themselves, so it can't depend on that reader ! */
void make_lexer(FiniteAutomaton& result) {
  std::string meta_chars_str = ".[]()|-^*+?";
  std::set<int> all_chars;
  for (int i = 0; i < NCHARS; ++i) all_chars.insert(i);
  std::set<int> nonmeta_chars = all_chars;
  for (int i = 0; i < Teuchos::size(meta_chars_str); ++i) {
    int meta_char = at(meta_chars_str, i);
    std::set<int>::iterator it = nonmeta_chars.find(get_symbol(meta_char));
    nonmeta_chars.erase(it);
  }
  FiniteAutomaton lex_nonmeta;
  make_set_nfa(lex_nonmeta, NCHARS, nonmeta_chars, TOK_CHAR);
  FiniteAutomaton lex_slash;
  make_char_single_nfa(lex_slash, '\\');
  FiniteAutomaton lex_any;
  make_set_nfa(lex_any, NCHARS, all_chars);
  FiniteAutomaton lex_escaped;
  concat(lex_escaped, lex_slash, lex_any, TOK_CHAR);
  FiniteAutomaton lex_char;
  unite(lex_char, lex_nonmeta, lex_escaped);
  FiniteAutomaton lex_metachars;
  for (int i = 0; i < Teuchos::size(meta_chars_str); ++i) {
    int token = TOK_CHAR + i + 1;
    if (i) {
      FiniteAutomaton lex_metachar;
      make_char_single_nfa(lex_metachar, at(meta_chars_str, i), token);
      unite(lex_metachars, lex_metachars, lex_metachar);
    } else {
      make_char_single_nfa(lex_metachars, at(meta_chars_str, i), token);
    }
  }
  unite(result, lex_metachars, lex_char);
  make_deterministic(result, result);
  simplify(result, result);
}

ReaderTablesPtr ask_reader_tables() {
  static ReaderTablesPtr ptr;
  if (ptr.strong_count() == 0) {
    RCP<ReaderTables> newptr(new ReaderTables());
    LanguagePtr lang = regex::ask_language();
    GrammarPtr grammar = make_grammar(*lang);
    newptr->parser = make_lalr1_parser(grammar);
    regex::make_lexer(newptr->lexer);
    newptr->indent_info.is_sensitive = false;
    newptr->indent_info.indent_token = -1;
    newptr->indent_info.dedent_token = -1;
    ptr = newptr;
  }
  return ptr;
}

LanguagePtr ask_language() {
  static LanguagePtr ptr;
  if (ptr.strong_count() == 0) {
    ptr.reset(new Language(make_language()));
  }
  return ptr;
}

void make_dfa(FiniteAutomaton& result, std::string const& name, std::string const& regex, int token) {
  using std::swap;
  regex::Reader reader(token);
  any result_any;
  try {
    reader.read_string(result_any, regex, name);
  } catch (const Teuchos::ParserFail& e) {
    std::stringstream ss;
    ss << e.what() << '\n';
    ss << "error: couldn't build DFA for token \"" << name << "\" regex \"" << regex << "\"\n";
    ss << "repeating with DebugReader:\n";
    DebugReader debug_reader(regex::ask_reader_tables(), ss);
    debug_reader.read_string(result_any, regex, name);
    throw ParserFail(ss.str());
  }
  swap(any_ref_cast<FiniteAutomaton>(result_any), result);
}

regex::Reader::Reader(int result_token_in):
  Teuchos::Reader(regex::ask_reader_tables()),
  result_token(result_token_in) {
}

void regex::Reader::at_shift(any& result, int token, std::string& text) {
  if (token != TOK_CHAR) return;
  if (Teuchos::size(text) == 1) {
    result = text[0];
  } else if (Teuchos::size(text) == 2) {
    TEUCHOS_ASSERT(text[0] == '\\');
    result = text[1];
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "BUG: regex char text is \"" << text << "\"\n");
  }
}

void regex::Reader::at_reduce(any& result_any, int production, std::vector<any>& rhs) {
  using std::swap;
  switch (production) {
    case PROD_REGEX: {
      swap(result_any, at(rhs, 0));
      FiniteAutomaton& result = any_ref_cast<FiniteAutomaton>(result_any);
      make_deterministic(result, result);
      simplify(result, result);
      return;
    }
    case PROD_UNION_DECAY:
    case PROD_CONCAT_DECAY:
    case PROD_QUAL_DECAY:
    case PROD_SET_ITEMS_DECAY:
    case PROD_SET_ITEM_RANGE: {
      swap(result_any, at(rhs, 0));
      return;
    }
    case PROD_UNION: {
      FiniteAutomaton& result = make_any_ref<FiniteAutomaton>(result_any);
      FiniteAutomaton& a = any_ref_cast<FiniteAutomaton>(at(rhs, 0));
      FiniteAutomaton& b = any_ref_cast<FiniteAutomaton>(at(rhs, 2));
      unite(result, a, b);
      return;
    }
    case PROD_CONCAT: {
      FiniteAutomaton& result = make_any_ref<FiniteAutomaton>(result_any);
      FiniteAutomaton& a = any_ref_cast<FiniteAutomaton>(at(rhs, 0));
      FiniteAutomaton& b = any_ref_cast<FiniteAutomaton>(at(rhs, 1));
      concat(result, a, b, result_token);
      return;
    }
    case PROD_STAR: {
      FiniteAutomaton& result = make_any_ref<FiniteAutomaton>(result_any);
      FiniteAutomaton& a = any_ref_cast<FiniteAutomaton>(at(rhs, 0));
      star(result, a, result_token);
      return;
    }
    case PROD_PLUS: {
      FiniteAutomaton& result = make_any_ref<FiniteAutomaton>(result_any);
      FiniteAutomaton& a = any_ref_cast<FiniteAutomaton>(at(rhs, 0));
      plus(result, a, result_token);
      return;
    }
    case PROD_MAYBE: {
      FiniteAutomaton& result = make_any_ref<FiniteAutomaton>(result_any);
      FiniteAutomaton& a = any_ref_cast<FiniteAutomaton>(at(rhs, 0));
      maybe(result, a, result_token);
      return;
    }
    case PROD_SINGLE_CHAR: {
      FiniteAutomaton& result = make_any_ref<FiniteAutomaton>(result_any);
      char c = any_cast<char>(at(rhs, 0));
      make_char_single_nfa(result, c, result_token);
      return;
    }
    case PROD_ANY: {
      FiniteAutomaton& result = make_any_ref<FiniteAutomaton>(result_any);
      make_range_nfa(result, NCHARS, 0, NCHARS - 1, result_token);
      return;
    }
    case PROD_SINGLE_SET: {
      FiniteAutomaton& result = make_any_ref<FiniteAutomaton>(result_any);
      std::set<char>& charset = any_ref_cast<std::set<char> >(at(rhs, 0));
      make_char_set_nfa(result, charset, result_token);
      return;
    }
    case PROD_PARENS_UNION: {
      swap(result_any, at(rhs, 1));
      return;
    }
    case PROD_SET_POSITIVE: {
      swap(result_any, at(rhs, 0));
      return;
    }
    case PROD_SET_NEGATIVE: {
      std::set<char>& result = make_any_ref<std::set<char> >(result_any);
      std::set<char> const& charset = any_ref_cast<std::set<char> >(at(rhs, 0));
      negate_set(result, charset);
      return;
    }
    case PROD_POSITIVE_SET: {
      swap(result_any, at(rhs, 1));
      return;
    }
    case PROD_NEGATIVE_SET: {
      swap(result_any, at(rhs, 2));
      return;
    }
    case PROD_SET_ITEMS_ADD: {
      std::set<char>& result = make_any_ref<std::set<char> >(result_any);
      std::set<char>& a = any_ref_cast<std::set<char> >(at(rhs, 0));
      std::set<char> const& b = any_ref_cast<std::set<char> >(at(rhs, 1));
      swap(result, a);
      unite_with(result, b);
      return;
    }
    case PROD_SET_ITEM_CHAR: {
      std::set<char>& result = make_any_ref<std::set<char> >(result_any);
      char c = any_cast<char>(at(rhs, 0));
      result.insert(c);
      return;
    }
    case PROD_RANGE: {
      std::set<char>& result = make_any_ref<std::set<char> >(result_any);
      char a = any_cast<char>(at(rhs, 0));
      char b = any_cast<char>(at(rhs, 2));
      for (char c = a; c <= b; ++c) {
        result.insert(c);
      }
      return;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "BUG: unexpected production " << production << '\n');
}

}  // namespace regex
}  // namespace Teuchos
