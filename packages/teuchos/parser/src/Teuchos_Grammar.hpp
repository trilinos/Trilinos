// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_GRAMMAR_HPP
#define TEUCHOS_GRAMMAR_HPP

#include <string>
#include <vector>

#include <Teuchos_RCP.hpp>

namespace Teuchos {

/* convention: symbols are numbered with all
   terminal symbols first, all non-terminal symbols after */

struct Grammar {
  typedef std::vector<int> RHS;
  struct Production {
    int lhs;
    RHS rhs;
  };
  typedef std::vector<Production> Productions;
  int nsymbols;
  int nterminals;
  Productions productions;
  std::vector<std::string> symbol_names;
};

typedef RCP<const Grammar> GrammarPtr;

int get_nnonterminals(Grammar const& g);
bool is_terminal(Grammar const& g, int symbol);
bool is_nonterminal(Grammar const& g, int symbol);
int as_nonterminal(Grammar const& g, int symbol);
int find_goal_symbol(Grammar const& g);
void add_end_terminal(Grammar& g);
int get_end_terminal(Grammar const& g);
void add_accept_production(Grammar& g);
int get_accept_production(Grammar const& g);
int get_accept_nonterminal(Grammar const& g);

std::ostream& operator<<(std::ostream& os, Grammar const& g);

}

#endif
