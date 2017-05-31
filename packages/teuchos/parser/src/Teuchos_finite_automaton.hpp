#ifndef TEUCHOS_FINITE_AUTOMATON_HPP
#define TEUCHOS_FINITE_AUTOMATON_HPP

#include <Teuchos_table.hpp>
#include <iosfwd>
#include <set>

namespace Teuchos {

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
  FiniteAutomaton() = default;
  FiniteAutomaton(int nsymbols_init, bool is_deterministic_init, int nstates_reserve);
};

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

FiniteAutomaton make_single_nfa(int nsymbols, int symbol, int token = 0);
FiniteAutomaton make_set_nfa(int nsymbols, std::set<int> const& accepted, int token = 0);
FiniteAutomaton make_range_nfa(int nsymbols, int range_start, int range_end, int token = 0);
FiniteAutomaton unite(FiniteAutomaton const& a, FiniteAutomaton const& b);
FiniteAutomaton concat(FiniteAutomaton const& a, FiniteAutomaton const& b, int token = 0);
FiniteAutomaton plus(FiniteAutomaton const& a, int token = 0);
FiniteAutomaton maybe(FiniteAutomaton const& a, int token = 0);
FiniteAutomaton star(FiniteAutomaton const& a, int token = 0);
FiniteAutomaton make_deterministic(FiniteAutomaton const& nfa);
FiniteAutomaton simplify_once(FiniteAutomaton const& fa);
FiniteAutomaton simplify(FiniteAutomaton const& fa);

FiniteAutomaton make_char_nfa(bool is_deterministic_init, int nstates_reserve);
void add_char_transition(FiniteAutomaton& fa, int from_state, char at_char, int to_state);
bool is_symbol(char c);
int get_symbol(char c);
char get_char(int symbol);
FiniteAutomaton make_char_set_nfa(std::set<char> const& accepted, int token = 0);
FiniteAutomaton make_char_range_nfa(char range_start, char range_end, int token = 0);
FiniteAutomaton make_char_single_nfa(char symbol_char, int token = 0);
std::set<char> negate_set(std::set<char> const& s);

std::ostream& operator<<(std::ostream& os, FiniteAutomaton const& fa);

}

#endif
