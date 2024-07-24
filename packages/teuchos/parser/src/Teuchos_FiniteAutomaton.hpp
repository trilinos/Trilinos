// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_FINITE_AUTOMATON_HPP
#define TEUCHOS_FINITE_AUTOMATON_HPP

#include <Teuchos_TableDecl.hpp>
#include <iosfwd>
#include <set>
#include <stdexcept>

namespace Teuchos {

#ifdef HAVE_TEUCHOSCORE_CXX11
extern template struct Table<int>;
#endif

/* This is basically a weird mix between a DFA and
   an NFA-epsilon. It is really a DFA that can have two extra
   epsilon symbols that it accepts transitions with.
   We can simulate epsilon-transitions to multiple new
   states by making trees of nodes connected by
   epsilon-transitions.

   by convention, the start state is state 0
 */
struct FiniteAutomaton {
  Table<int> table;
  std::vector<int> accepted_tokens;
  bool is_deterministic;
  FiniteAutomaton() {}
  FiniteAutomaton(int nsymbols_init, bool is_deterministic_init, int nstates_reserve);
  void swap(FiniteAutomaton& other);
};

/* NOTE: this is only needed by Teuchos::any to support its non-standard operator== */
inline bool operator==(FiniteAutomaton const&, FiniteAutomaton const&) {
  return false;
}

inline void swap(FiniteAutomaton& a, FiniteAutomaton& b) { a.swap(b); }

int get_nstates(FiniteAutomaton const& fa);
int get_nsymbols(FiniteAutomaton const& fa);
bool get_determinism(FiniteAutomaton const& fa);
int get_epsilon0(FiniteAutomaton const& fa);
int get_epsilon1(FiniteAutomaton const& fa);
int add_state(FiniteAutomaton& fa);
void add_transition(FiniteAutomaton& fa, int from_state, int at_symbol, int to_state);
void add_accept(FiniteAutomaton& fa, int state, int token);
void remove_accept(FiniteAutomaton& fa, int state);
int step(FiniteAutomaton const& fa, int state, int symbol);
int accepts(FiniteAutomaton const& fa, int state);
int get_nsymbols_eps(FiniteAutomaton const& fa);
void append_states(FiniteAutomaton& fa, FiniteAutomaton const& other);

void make_single_nfa(FiniteAutomaton& result, int nsymbols, int symbol, int token = 0);
void make_set_nfa(FiniteAutomaton& result, int nsymbols, std::set<int> const& accepted, int token = 0);
void make_range_nfa(FiniteAutomaton& result, int nsymbols, int range_start, int range_end, int token = 0);
void unite(FiniteAutomaton& result, FiniteAutomaton const& a, FiniteAutomaton const& b);
void concat(FiniteAutomaton& result, FiniteAutomaton const& a, FiniteAutomaton const& b, int token = 0);
void plus(FiniteAutomaton& result, FiniteAutomaton const& a, int token = 0);
void maybe(FiniteAutomaton& result, FiniteAutomaton const& a, int token = 0);
void star(FiniteAutomaton& result, FiniteAutomaton const& a, int token = 0);
void make_deterministic(FiniteAutomaton& result, FiniteAutomaton& nfa);
void simplify_once(FiniteAutomaton& result, FiniteAutomaton const& fa);
void simplify(FiniteAutomaton& result, FiniteAutomaton const& fa);

void add_char_transition(FiniteAutomaton& fa, int from_state, char at_char, int to_state);
bool is_symbol(char c);
int get_symbol(char c);
char get_char(int symbol);
void make_char_set_nfa(FiniteAutomaton& result, std::set<char> const& accepted, int token = 0);
void make_char_range_nfa(FiniteAutomaton& result, char range_start, char range_end, int token = 0);
void make_char_single_nfa(FiniteAutomaton& result, char symbol_char, int token = 0);
void negate_set(std::set<char>& result, std::set<char> const& s);

std::ostream& operator<<(std::ostream& os, FiniteAutomaton const& fa);

}

#endif
