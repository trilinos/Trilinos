#include "Teuchos_parser.hpp"

namespace Teuchos {

Parser::Parser(GrammarPtr g, int nstates_reserve):
  grammar(g),
  terminal_table(g->nterminals, nstates_reserve),
  nonterminal_table(get_nnonterminals(*g), nstates_reserve) {
}

int get_nstates(Parser const& p) {
  return get_nrows(p.terminal_table);
}

int add_state(Parser& p) {
  auto state = get_nstates(p);
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
  assert(at(p.terminal_table, state, terminal).kind == ACTION_NONE);
  assert(action.kind != ACTION_NONE);
  if (action.kind == ACTION_SHIFT) {
    assert(0 <= action.next_state);
    assert(action.next_state < get_nstates(p));
  } else {
    assert(0 <= action.production);
    assert(action.production < size(p.grammar->productions));
  }
  at(p.terminal_table, state, terminal) = action;
}

void add_nonterminal_action(Parser& p, int state, int nonterminal, int next_state) {
  assert(0 <= next_state);
  assert(next_state < get_nstates(p));
  assert(at(p.nonterminal_table, state, nonterminal) == -1);
  at(p.nonterminal_table, state, nonterminal) = next_state;
}

Action const& get_action(Parser const& p, int state, int terminal) {
  return at(p.terminal_table, state, terminal);
}

int execute_action(Parser const& p, std::vector<int>& stack, Action const& action) {
  assert(action.kind != ACTION_NONE);
  if (action.kind == ACTION_SHIFT) {
    stack.push_back(action.next_state);
  } else {
    auto& prod = at(p.grammar->productions, action.production);
    for (int i = 0; i < size(prod.rhs); ++i) stack.pop_back();
    assert(p.grammar.get());
    auto& grammar = *(p.grammar);
    auto nt = as_nonterminal(grammar, prod.lhs);
    assert(!stack.empty());
    auto next_state = at(p.nonterminal_table, stack.back(), nt);
    stack.push_back(next_state);
  }
  return stack.back();
}

GrammarPtr const& get_grammar(Parser const& p) { return p.grammar; }

} // end namespace Teuchos
