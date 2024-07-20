// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Language.hpp"

#include <set>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstdarg>

#include "Teuchos_vector.hpp"
#include "Teuchos_regex.hpp"
#include "Teuchos_Parser.hpp"

namespace Teuchos {

void Language::Token::operator()(std::string const& name_in, std::string const& regex_in) {
  name = name_in;
  regex = regex_in;
}

Language::RHSBuilder::RHSBuilder(Production& prod_in):
  prod(prod_in) {
}

Language::RHSBuilder& Language::RHSBuilder::operator,(std::string const& rhs_item) {
  prod.rhs.push_back(rhs_item);
  return *this;
}

Language::RHSBuilder& Language::RHSBuilder::operator>>(std::string const& rhs_item) {
  prod.rhs.push_back(rhs_item);
  return *this;
}

Language::RHSBuilder Language::Production::operator()(std::string const& lhs_in) {
  lhs = lhs_in;
  return Language::RHSBuilder(*this);
}

GrammarPtr make_grammar(Language const& language) {
  std::map<std::string, int> symbol_map;
  int nterminals = 0;
  for (Language::Tokens::const_iterator it = language.tokens.begin();
       it != language.tokens.end(); ++it) {
    const Language::Token& token = *it;
    TEUCHOS_TEST_FOR_EXCEPTION(token.name.empty(), ParserFail,
        "ERROR: token " << it - language.tokens.begin() << " has an empty name\n");
    symbol_map[token.name] = nterminals++;
  }
  int nsymbols = nterminals;
  for (Language::Productions::const_iterator it = language.productions.begin();
       it != language.productions.end(); ++it) {
    const Language::Production& production = *it;
    TEUCHOS_TEST_FOR_EXCEPTION(production.lhs.empty(), ParserFail,
        "ERROR: production " << it - language.productions.begin() << " has an empty LHS name\n");
    if (symbol_map.count(production.lhs)) continue;
    symbol_map[production.lhs] = nsymbols++;
  }
  RCP<Grammar> out(new Grammar());
  out->nsymbols = nsymbols;
  out->nterminals = nterminals;
  for (Language::Productions::const_iterator it = language.productions.begin();
       it != language.productions.end(); ++it) {
    const Language::Production& lang_prod = *it;
    out->productions.push_back(Grammar::Production());
    Grammar::Production& gprod = out->productions.back();
    TEUCHOS_ASSERT(symbol_map.count(lang_prod.lhs));
    gprod.lhs = symbol_map[lang_prod.lhs];
    for (Language::RHS::const_iterator it2 = lang_prod.rhs.begin();
         it2 != lang_prod.rhs.end(); ++it2) {
      const std::string& lang_symb = *it2;
      TEUCHOS_TEST_FOR_EXCEPTION(!symbol_map.count(lang_symb), ParserFail,
           "RHS entry \"" << lang_symb <<
           "\" is neither a nonterminal (LHS of a production) nor a token!\n");
      gprod.rhs.push_back(symbol_map[lang_symb]);
    }
  }
  out->symbol_names = make_vector<std::string>(nsymbols);
  for (std::map<std::string, int>::const_iterator it = symbol_map.begin();
       it != symbol_map.end(); ++it) {
    const std::pair<std::string, int>& pair = *it;
    at(out->symbol_names, pair.second) = pair.first;
  }
  add_end_terminal(*out);
  add_accept_production(*out);
  return out;
}

std::ostream& operator<<(std::ostream& os, Language const& lang) {
  for (Language::Tokens::const_iterator it = lang.tokens.begin();
       it != lang.tokens.end(); ++it) {
    const Language::Token& token = *it;
    os << "token " << token.name << " regex \'" << token.regex << "\'\n";
  }
  std::set<std::string> nonterminal_set;
  std::vector<std::string> nonterminal_list;
  for (Language::Productions::const_iterator it = lang.productions.begin();
       it != lang.productions.end(); ++it) {
    const Language::Production& prod = *it;
    if (!nonterminal_set.count(prod.lhs)) {
      nonterminal_set.insert(prod.lhs);
      nonterminal_list.push_back(prod.lhs);
    }
  }
  for (std::vector<std::string>::const_iterator it = nonterminal_list.begin();
       it != nonterminal_list.end(); ++it) {
    const std::string& nonterminal = *it;
    std::stringstream ss;
    ss << nonterminal << " ::=";
    std::string lead = ss.str();
    os << lead;
    for (std::string::iterator it2 = lead.begin(); it2 != lead.end(); ++it2) {
      *it2 = ' ';
    }
    bool first = true;
    for (Language::Productions::const_iterator it2 = lang.productions.begin();
         it2 != lang.productions.end(); ++it2) {
      const Language::Production& prod = *it2;
      if (prod.lhs != nonterminal) continue;
      if (first) first = false;
      else os << " |\n" << lead;
      for (Language::RHS::const_iterator it3 = prod.rhs.begin();
           it3 != prod.rhs.end(); ++it3) {
        const std::string& symb = *it3;
        if (symb == "|") os << " '|'";
        else os << " " << symb;
      }
    }
    os << "\n";
  }
  os << "\n";
  return os;
}

void make_lexer(FiniteAutomaton& result, Language const& language) {
  using std::swap;
  for (int i = 0; i < Teuchos::size(language.tokens); ++i) {
    const std::string& name = at(language.tokens, i).name;
    const std::string& regex = at(language.tokens, i).regex;
    if (i == 0) {
      regex::make_dfa(result, name, regex, i);
    } else {
      FiniteAutomaton b;
      regex::make_dfa(b, name, regex, i);
      unite(result, result, b);
    }
  }
  make_deterministic(result, result);
  simplify(result, result);
}

static void make_indent_info(IndentInfo& out, Language const& language) {
  out.is_sensitive = false;
  out.indent_token = -1;
  out.dedent_token = -1;
  out.newline_token = -1;
  for (int tok_i = 0; tok_i < Teuchos::size(language.tokens); ++tok_i) {
    const Language::Token& token = at(language.tokens, tok_i);
    if (token.name == "INDENT") {
      TEUCHOS_TEST_FOR_EXCEPTION(out.indent_token != -1, ParserFail,
          "error: Language has two or more INDENT tokens\n");
      out.indent_token = tok_i;
      out.is_sensitive = true;
    } else if (token.name == "DEDENT") {
      TEUCHOS_TEST_FOR_EXCEPTION(out.dedent_token != -1, ParserFail,
          "error: Language has two or more DEDENT tokens\n");
      out.dedent_token = tok_i;
    } else if (token.name == "NEWLINE") {
      TEUCHOS_TEST_FOR_EXCEPTION(out.newline_token != -1, ParserFail,
          "error: Language has two or more NEWLINE tokens\n");
      out.newline_token = tok_i;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(out.is_sensitive && out.indent_token == -1,
      ParserFail,
      "error: Indentation-sensitive language has no INDENT token\n");
  TEUCHOS_TEST_FOR_EXCEPTION(out.is_sensitive && out.dedent_token == -1,
      ParserFail,
      "error: Indentation-sensitive language has no DEDENT token\n");
  TEUCHOS_TEST_FOR_EXCEPTION(out.is_sensitive && out.newline_token == -1,
      ParserFail,
      "error: Indentation-sensitive language has no NEWLINE token\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
      (out.indent_token < out.newline_token ||
       out.dedent_token < out.newline_token),
      ParserFail,
      "error: NEWLINE needs to come before all other indent tokens\n");
}

ReaderTablesPtr make_reader_tables(Language const& language) {
  RCP<ReaderTables> out(new ReaderTables());
  make_lexer(out->lexer, language);
  make_indent_info(out->indent_info, language);
  GrammarPtr grammar = make_grammar(language);
  out->parser = make_lalr1_parser(grammar);
  return out;
}

}
