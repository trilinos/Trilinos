// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Grammar.hpp"

#include <Teuchos_Assert.hpp>
#include <set>
#include <iostream>

#include "Teuchos_vector.hpp"
#include "Teuchos_Parser.hpp"

namespace Teuchos {

int get_nnonterminals(Grammar const& g) {
  return g.nsymbols - g.nterminals;
}

bool is_terminal(Grammar const& g, int symbol) {
  TEUCHOS_DEBUG_ASSERT(0 <= symbol);
  TEUCHOS_DEBUG_ASSERT(symbol <= g.nsymbols);
  return symbol < g.nterminals;
}

bool is_nonterminal(Grammar const& g, int symbol) {
  return !is_terminal(g, symbol);
}

int as_nonterminal(Grammar const& g, int symbol) {
  return symbol - g.nterminals;
}

int find_goal_symbol(Grammar const& g) {
  std::set<int> nonterminals_in_rhss;
  for (int i = 0; i < Teuchos::size(g.productions); ++i) {
    const Grammar::Production& p = at(g.productions, i);
    for (int j = 0; j < Teuchos::size(p.rhs); ++j) {
      const int s = at(p.rhs, j);
      TEUCHOS_DEBUG_ASSERT(0 <= s);
      if (is_nonterminal(g, s)) nonterminals_in_rhss.insert(s);
    }
  }
  int result = -1;
  for (int s = g.nterminals; s < g.nsymbols; ++s) {
    if (!nonterminals_in_rhss.count(s)) {
      TEUCHOS_TEST_FOR_EXCEPTION(result != -1, ParserFail,
          "ERROR: there is more than one root nonterminal ("
          << at(g.symbol_names, result) << " [" << result << "] and "
          << at(g.symbol_names, s) << " [" << s << "]) in this grammar\n");
      result = s;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(result == -1, ParserFail,
      "ERROR: the root nonterminal is unclear for this grammar\n"
      "usually this means all nonterminals appear in the RHS of a production\n"
      "and can be fixed by adding a new nonterminal root2, root2 -> root\n");
  return result;
}

void add_end_terminal(Grammar& g) {
  for (int i = 0; i < Teuchos::size(g.productions); ++i) {
    Grammar::Production& prod = at(g.productions, i);
    if (is_nonterminal(g, prod.lhs)) prod.lhs++;
    for (int j = 0; j < Teuchos::size(prod.rhs); ++j) {
      int& rhs_symb = at(prod.rhs, j);
      if (is_nonterminal(g, rhs_symb)) rhs_symb++;
    }
  }
  g.symbol_names.insert(g.symbol_names.begin() + g.nterminals, "EOF");
  g.nterminals++;
  g.nsymbols++;
}

int get_end_terminal(Grammar const& g) {
  return g.nterminals - 1;
}

void add_accept_production(Grammar& g) {
  int goal_symbol = find_goal_symbol(g);
  Grammar::Production p;
  p.lhs = g.nsymbols;
  p.rhs.push_back(goal_symbol);
  g.productions.push_back(p);
  g.symbol_names.push_back("ACCEPT");
  g.nsymbols++;
}

int get_accept_production(Grammar const& g) {
  return Teuchos::size(g.productions) - 1;
}

int get_accept_nonterminal(Grammar const& g) {
  return g.nsymbols - 1;
}

std::ostream& operator<<(std::ostream& os, Grammar const& g) {
  os << "symbols:\n";
  for (int i = 0; i < Teuchos::size(g.symbol_names); ++i) {
    os << i << ": " << at(g.symbol_names, i) << "\n";
  }
  os << "productions:\n";
  for (int i = 0; i < Teuchos::size(g.productions); ++i) {
    const Grammar::Production& prod = at(g.productions, i);
    os << i << ": " << prod.lhs << " ::=";
    for (int j = 0; j < Teuchos::size(prod.rhs); ++j) {
      int symb = at(prod.rhs, j);
      os << ' ' << symb;
    }
    os << '\n';
  }
  os << '\n';
  return os;
}

}
