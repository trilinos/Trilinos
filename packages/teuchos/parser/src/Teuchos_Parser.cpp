// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Parser.hpp"

#include "Teuchos_Table.hpp"
#include "Teuchos_make_lalr1_parser.hpp"

namespace Teuchos {

template struct Table<Action>;

Parser::Parser(GrammarPtr g, int nstates_reserve):
  grammar(g),
  terminal_table(g->nterminals, nstates_reserve),
  nonterminal_table(get_nnonterminals(*g), nstates_reserve) {
}

int get_nstates(Parser const& p) {
  return get_nrows(p.terminal_table);
}

int add_state(Parser& p) {
  int state = get_nstates(p);
  resize(p.terminal_table, state + 1, get_ncols(p.terminal_table));
  resize(p.nonterminal_table, state + 1, get_ncols(p.nonterminal_table));
  for (int t = 0; t < p.grammar->nterminals; ++t) {
    Action action;
    action.kind = ACTION_NONE;
    at(p.terminal_table, state, t) = action;
  }
  for (int nt = 0; nt < get_nnonterminals(*(p.grammar)); ++nt) {
    at(p.nonterminal_table, state, nt) = -1;
  }
  return state;
}

void add_terminal_action(Parser& p, int state, int terminal, Action action) {
  TEUCHOS_ASSERT(at(p.terminal_table, state, terminal).kind == ACTION_NONE);
  TEUCHOS_ASSERT(action.kind != ACTION_NONE);
  if (action.kind == ACTION_SHIFT) {
    TEUCHOS_ASSERT(0 <= action.next_state);
    TEUCHOS_ASSERT(action.next_state < get_nstates(p));
  } else {
    TEUCHOS_ASSERT(0 <= action.production);
    TEUCHOS_ASSERT(action.production < Teuchos::size(p.grammar->productions));
  }
  at(p.terminal_table, state, terminal) = action;
}

void add_nonterminal_action(Parser& p, int state, int nonterminal, int next_state) {
  TEUCHOS_ASSERT(0 <= next_state);
  TEUCHOS_ASSERT(next_state < get_nstates(p));
  TEUCHOS_ASSERT(at(p.nonterminal_table, state, nonterminal) == -1);
  at(p.nonterminal_table, state, nonterminal) = next_state;
}

Action const& get_action(Parser const& p, int state, int terminal) {
  return at(p.terminal_table, state, terminal);
}

int execute_action(Parser const& p, std::vector<int>& stack, Action const& action) {
  TEUCHOS_ASSERT(action.kind != ACTION_NONE);
  if (action.kind == ACTION_SHIFT) {
    stack.push_back(action.next_state);
  } else {
    const Grammar::Production& prod = at(p.grammar->productions, action.production);
    for (int i = 0; i < Teuchos::size(prod.rhs); ++i) stack.pop_back();
    TEUCHOS_ASSERT(p.grammar.get());
    const Grammar& grammar = *(p.grammar);
    int nt = as_nonterminal(grammar, prod.lhs);
    TEUCHOS_ASSERT(!stack.empty());
    int next_state = at(p.nonterminal_table, stack.back(), nt);
    stack.push_back(next_state);
  }
  return stack.back();
}

GrammarPtr const& get_grammar(Parser const& p) { return p.grammar; }

ParserFail::ParserFail(const std::string& msg):
  std::invalid_argument(msg) {
}

} // end namespace Teuchos
