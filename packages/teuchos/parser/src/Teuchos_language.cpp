#include "Teuchos_language.hpp"

#include <set>
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "Teuchos_vector.hpp"
#include "Teuchos_regex.hpp"
#include "Teuchos_build_parser.hpp"

namespace Teuchos {

GrammarPtr build_grammar(Language const& language) {
  std::map<std::string, int> symbol_map;
  int nterminals = 0;
  for (auto& token : language.tokens) {
    symbol_map[token.name] = nterminals++;
  }
  int nsymbols = nterminals;
  for (auto& production : language.productions) {
    if (symbol_map.count(production.lhs)) continue;
    symbol_map[production.lhs] = nsymbols++;
  }
  Grammar out;
  out.nsymbols = nsymbols;
  out.nterminals = nterminals;
  for (auto& lang_prod : language.productions) {
    Grammar::Production gprod;
    assert(symbol_map.count(lang_prod.lhs));
    gprod.lhs = symbol_map[lang_prod.lhs];
    for (auto& lang_symb : lang_prod.rhs) {
      if (!symbol_map.count(lang_symb)) {
        std::cerr << "RHS entry \"" << lang_symb << "\" is neither a nonterminal (LHS of a production) nor a token!\n";
        abort();
      }
      gprod.rhs.push_back(symbol_map[lang_symb]);
    }
    out.productions.emplace_back(std::move(gprod));
  }
  out.symbol_names = make_vector<std::string>(nsymbols);
  for (auto& pair : symbol_map) {
    at(out.symbol_names, pair.second) = pair.first;
  }
  add_end_terminal(out);
  add_accept_production(out);
  return std::make_shared<Grammar>(std::move(out));
}

std::ostream& operator<<(std::ostream& os, Language const& lang) {
  for (auto& token : lang.tokens) {
    os << "token " << token.name << " regex \'" << token.regex << "\'\n";
  }
  std::set<std::string> nonterminal_set;
  std::vector<std::string> nonterminal_list;
  for (auto& prod : lang.productions) {
    if (!nonterminal_set.count(prod.lhs)) {
      nonterminal_set.insert(prod.lhs);
      nonterminal_list.push_back(prod.lhs);
    }
  }
  for (auto& nonterminal : nonterminal_list) {
    std::stringstream ss;
    ss << nonterminal << " ::=";
    auto lead = ss.str();
    os << lead;
    for (auto& c : lead) c = ' ';
    bool first = true;
    for (auto& prod : lang.productions) {
      if (prod.lhs != nonterminal) continue;
      if (first) first = false;
      else os << " |\n" << lead;
      for (auto& symb : prod.rhs) {
        if (symb == "|") os << " '|'";
        else os << " " << symb;
      }
    }
    os << "\n";
  }
  os << "\n";
  return os;
}

FiniteAutomaton build_lexer(Language const& language) {
  FiniteAutomaton lexer;
  for (int i = 0; i < size(language.tokens); ++i) {
    auto& name = at(language.tokens, i).name;
    auto& regex = at(language.tokens, i).regex;
    if (i == 0) {
      lexer = regex::build_dfa(name, regex, i);
    } else {
      lexer = unite(lexer, regex::build_dfa(name, regex, i));
    }
  }
  lexer = simplify(make_deterministic(lexer));
  return lexer;
}

static IndentInfo build_indent_info(Language const& language) {
  IndentInfo out;
  out.is_sensitive = false;
  out.indent_token = -1;
  out.dedent_token = -1;
  out.eqdent_token = -1;
  out.nodent_token = -1;
  for (int tok_i = 0; tok_i < size(language.tokens); ++tok_i) {
    auto& token = at(language.tokens, tok_i);
    if (token.regex == "]INDENT[") {
      if (out.indent_token != -1) {
        std::cerr << "ERROR: Language has two or more ]INDENT[ tokens\n";
        abort();
      }
      out.indent_token = tok_i;
      out.is_sensitive = true;
    } else if (token.regex == "]DEDENT[") {
      if (out.dedent_token != -1) {
        std::cerr << "ERROR: Language has two or more ]INDENT[ tokens\n";
        abort();
      }
      out.dedent_token = tok_i;
    } else if (token.regex == "]EQDENT[") {
      if (out.eqdent_token != -1) {
        std::cerr << "ERROR: Language has two or more ]EQDENT[ tokens\n";
        abort();
      }
      out.eqdent_token = tok_i;
    } else if (token.regex == "]NODENT[") {
      if (out.nodent_token != -1) {
        std::cerr << "ERROR: Language has two or more ]EQDENT[ tokens\n";
        abort();
      }
      out.nodent_token = tok_i;
    }
  }
  if (out.is_sensitive && out.indent_token == -1) {
    std::cerr << "ERROR: Indentation-sensitive language has no ]INDENT[ token\n";
    abort();
  }
  if (out.is_sensitive && out.dedent_token == -1) {
    std::cerr << "ERROR: Indentation-sensitive language has no ]DEDENT[ token\n";
    abort();
  }
  if (out.is_sensitive && out.eqdent_token == -1) {
    std::cerr << "ERROR: Indentation-sensitive language has no ]EQDENT[ token\n";
    abort();
  }
  if (out.is_sensitive && out.nodent_token == -1) {
    std::cerr << "ERROR: Indentation-sensitive language has no ]NODENT[ token\n";
    abort();
  }
  if (out.indent_token < out.nodent_token ||
      out.dedent_token < out.nodent_token ||
      out.eqdent_token < out.nodent_token) {
    std::cerr << "ERROR: ]NODENT[ needs to come before all other indent tokens\n";
    abort();
  }
  return out;
}

ReaderTablesPtr build_reader_tables(Language const& language) {
  auto lexer = build_lexer(language);
  auto indent_info = build_indent_info(language);
  auto grammar = build_grammar(language);
  auto parser = accept_parser(build_lalr1_parser(grammar));
  return ReaderTablesPtr(new ReaderTables({parser, lexer, indent_info}));
}

}
