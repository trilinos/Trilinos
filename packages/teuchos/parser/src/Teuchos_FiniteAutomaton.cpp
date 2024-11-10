// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_FiniteAutomaton.hpp"

#include <set>
#include <map>
#include <queue>
#include <utility>
#include <memory>
#include <limits>

#include "Teuchos_chartab.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Table.hpp"

namespace Teuchos {

template struct Table<int>;

FiniteAutomaton::FiniteAutomaton(int nsymbols_init, bool is_deterministic_init,
    int nstates_reserve):
  table(nsymbols_init + (is_deterministic_init ? 0 : 2), nstates_reserve),
  is_deterministic(is_deterministic_init)
{
  reserve(accepted_tokens, nstates_reserve);
}

void FiniteAutomaton::swap(FiniteAutomaton& other) {
  using std::swap;
  swap(table, other.table);
  swap(accepted_tokens, other.accepted_tokens);
  swap(is_deterministic, other.is_deterministic);
}

int get_nstates(FiniteAutomaton const& fa) {
  return get_nrows(fa.table);
}

int get_nsymbols(FiniteAutomaton const& fa) {
  return get_ncols(fa.table) - (fa.is_deterministic ? 0 : 2);
}

bool get_determinism(FiniteAutomaton const& fa) {
  return fa.is_deterministic;
}

int get_epsilon0(FiniteAutomaton const& fa) {
  TEUCHOS_ASSERT(!fa.is_deterministic);
  return get_ncols(fa.table) - 2;
}

int get_epsilon1(FiniteAutomaton const& fa) {
  TEUCHOS_ASSERT(!fa.is_deterministic);
  return get_ncols(fa.table) - 1;
}

int add_state(FiniteAutomaton& fa) {
  int state = get_nstates(fa);
  resize(fa.table, state + 1, get_ncols(fa.table));
  for (int j = 0; j < get_ncols(fa.table); ++j) {
    at(fa.table, state, j) = -1;
  }
  fa.accepted_tokens.push_back(-1);
  return state;
}

void add_transition(FiniteAutomaton& fa, int from_state, int at_symbol, int to_state) {
  TEUCHOS_ASSERT(0 <= to_state);
  TEUCHOS_ASSERT(to_state < get_nstates(fa));
  TEUCHOS_ASSERT(0 <= at_symbol);
  TEUCHOS_ASSERT(at_symbol < get_ncols(fa.table)); // allow setting epsilon transitions
  TEUCHOS_ASSERT(at(fa.table, from_state, at_symbol) == -1);
  at(fa.table, from_state, at_symbol) = to_state;
}

void add_accept(FiniteAutomaton& fa, int state, int token) {
  TEUCHOS_ASSERT(0 <= token);
  at(fa.accepted_tokens, state) = token;
}

void remove_accept(FiniteAutomaton& fa, int state) {
  at(fa.accepted_tokens, state) = -1;
}

int step(FiniteAutomaton const& fa, int state, int symbol) {
  TEUCHOS_ASSERT(0 <= state);
  TEUCHOS_ASSERT(state < get_nstates(fa));
  TEUCHOS_ASSERT(0 <= symbol);
  TEUCHOS_ASSERT(symbol < get_ncols(fa.table)); // allow getting epsilon transitions
  return at(fa.table, state, symbol);
}

int accepts(FiniteAutomaton const& fa, int state) {
  return at(fa.accepted_tokens, state);
}

int get_nsymbols_eps(FiniteAutomaton const& fa) {
  return get_ncols(fa.table);
}

void append_states(FiniteAutomaton& fa, FiniteAutomaton const& other) {
  TEUCHOS_ASSERT(get_nsymbols(other) == get_nsymbols(fa));
  bool other_determ = get_determinism(other);
  if (!other_determ) TEUCHOS_ASSERT(!fa.is_deterministic);
  int offset = get_nstates(fa);
  for (int other_state = 0; other_state < get_nstates(other); ++other_state) {
    int my_state = add_state(fa);
    int token = accepts(other, other_state);
    if (0 <= token) add_accept(fa, my_state, token);
  }
  for (int other_state = 0; other_state < get_nstates(other); ++other_state) {
    int my_state = other_state + offset;
    for (int symbol = 0; symbol < get_nsymbols_eps(other); ++symbol) {
      int other_next = step(other, other_state, symbol);
      if (other_next < 0) continue;
      int my_next = other_next + offset;
      add_transition(fa, my_state, symbol, my_next);
    }
  }
}

void make_single_nfa(FiniteAutomaton& result, int nsymbols, int symbol, int token) {
  make_range_nfa(result, nsymbols, symbol, symbol, token);
}

void make_set_nfa(FiniteAutomaton& result, int nsymbols, std::set<int> const& accepted, int token) {
  using std::swap;
  FiniteAutomaton out(nsymbols, true, 2);
  int start_state = add_state(out);
  int accept_state = add_state(out);
  for (std::set<int>::const_iterator it = accepted.begin(); it != accepted.end(); ++it) {
    int i = *it;
    add_transition(out, start_state, i, accept_state);
  }
  add_accept(out, accept_state, token);
  swap(result, out);
}

void make_range_nfa(FiniteAutomaton& result, int nsymbols, int range_start, int range_end, int token) {
  using std::swap;
  TEUCHOS_ASSERT(0 <= range_start);
  TEUCHOS_ASSERT(range_start <= range_end);
  TEUCHOS_ASSERT(range_end <= nsymbols);
  FiniteAutomaton out(nsymbols, true, 2);
  int start_state = add_state(out);
  int accept_state = add_state(out);
  for (int i = range_start; i <= range_end; ++i) {
    add_transition(out, start_state, i, accept_state);
  }
  add_accept(out, accept_state, token);
  swap(result, out);
}

void unite(FiniteAutomaton& result, FiniteAutomaton const& a, FiniteAutomaton const& b) {
  using std::swap;
  int nsymbols = get_nsymbols(a);
  FiniteAutomaton out(nsymbols, false, 1 + get_nstates(a) + get_nstates(b));
  int start_state = add_state(out);
  int a_offset = get_nstates(out);
  append_states(out, a);
  int b_offset = get_nstates(out);
  append_states(out, b);
  int epsilon0 = get_epsilon0(out);
  int epsilon1 = get_epsilon1(out);
  add_transition(out, start_state, epsilon0, a_offset);
  add_transition(out, start_state, epsilon1, b_offset);
  using std::swap;
  swap(out, result);
}

void concat(FiniteAutomaton& result, FiniteAutomaton const& a, FiniteAutomaton const& b, int token) {
  int nsymbols = get_nsymbols(a);
  FiniteAutomaton out(nsymbols, false, get_nstates(a) + get_nstates(b));
  append_states(out, a);
  int b_offset = get_nstates(out);
  append_states(out, b);
  int epsilon0 = get_epsilon0(out);
  for (int i = 0; i < get_nstates(a); ++i) {
    if (accepts(a, i) != -1) {
      add_transition(out, i, epsilon0, b_offset);
      remove_accept(out, i);
    }
  }
  for (int i = 0; i < get_nstates(b); ++i) {
    if (accepts(b, i) != -1) {
      add_accept(out, i + b_offset, token);
    }
  }
  swap(result, out);
}

void plus(FiniteAutomaton& result, FiniteAutomaton const& a, int token) {
  using std::swap;
  FiniteAutomaton out(get_nsymbols(a), false, get_nstates(a) + 1);
  append_states(out, a);
  int new_accept_state = add_state(out);
  add_accept(out, new_accept_state, token);
  int epsilon0 = get_epsilon0(out);
  int epsilon1 = get_epsilon1(out);
  for (int i = 0; i < get_nstates(a); ++i) {
    if (accepts(a, i) != -1) {
      add_transition(out, i, epsilon0, new_accept_state);
      /* we follow a convention that accepting
         states should not have epsilon transitions */
      add_transition(out, i, epsilon1, 0);
      remove_accept(out, i);
    }
  }
  swap(result, out);
}

void maybe(FiniteAutomaton& result, FiniteAutomaton const& a, int token) {
  using std::swap;
  FiniteAutomaton out(get_nsymbols(a), false, get_nstates(a) + 2);
  int new_start_state = add_state(out);
  int offset = get_nstates(out);
  append_states(out, a);
  int new_accept_state = add_state(out);
  int epsilon0 = get_epsilon0(out);
  int epsilon1 = get_epsilon1(out);
  add_transition(out, new_start_state, epsilon1, offset);
  /* form an epsilon0 linked list of new start state,
     all old accepting states, and new accepting state */
  int last = new_start_state;
  for (int i = 0; i < get_nstates(a); ++i) {
    if (accepts(a, i) != -1) {
      add_transition(out, last, epsilon0, i + offset);
      remove_accept(out, i + offset);
      last = i + offset;
    }
  }
  add_transition(out, last, epsilon0, new_accept_state);
  add_accept(out, new_accept_state, token);
  swap(result, out);
}

void star(FiniteAutomaton& result, FiniteAutomaton const& a, int token) {
  using std::swap;
  plus(result, a, token);
  maybe(result, result, token);
}

typedef std::set<int> StateSet;

static StateSet step(StateSet const& ss, int symbol, FiniteAutomaton const& fa) {
  StateSet next_ss;
  for (StateSet::const_iterator it = ss.begin(); it != ss.end(); ++it) {
    int state = *it;
    int next_state = step(fa, state, symbol);
    if (next_state != -1) next_ss.insert(next_state);
  }
  return next_ss;
}

typedef std::queue<int> StateQueue;

static StateSet get_epsilon_closure(StateSet ss, FiniteAutomaton const& fa) {
  StateQueue q;
  for (StateSet::const_iterator it = ss.begin(); it != ss.end(); ++it) {
    int state = *it;
    q.push(state);
  }
  int epsilon0 = get_epsilon0(fa);
  int epsilon1 = get_epsilon1(fa);
  while (!q.empty()) {
    int state = q.front(); q.pop();
    for (int epsilon = epsilon0; epsilon <= epsilon1; ++epsilon) {
      int next_state = step(fa, state, epsilon);
      if (next_state == -1) continue;
      if (!ss.count(next_state)) {
        ss.insert(next_state);
        q.push(next_state);
      }
    }
  }
  return ss;
}

struct StateSetPtrLess {
  bool operator()(StateSet* a, StateSet* b) const {
    return *a < *b;
  }
};

typedef std::map<StateSet*,int,StateSetPtrLess> StateSetPtr2State;
typedef RCP<StateSet> StateSetPtr;
typedef std::vector<StateSetPtr> StateSetUniqPtrVector;

static void add_back(StateSetUniqPtrVector& ssupv, StateSet& ss) {
  using std::swap;
  StateSetPtr ptr(new StateSet());
  swap(*ptr, ss);
  ssupv.push_back(ptr);
}

/* powerset construction, NFA -> DFA */
void make_deterministic(FiniteAutomaton& result, FiniteAutomaton& nfa) {
  using std::swap;
  if (get_determinism(nfa)) {
    swap(result, nfa);
    return;
  }
  StateSetPtr2State ssp2s;
  StateSetUniqPtrVector ssupv;
  FiniteAutomaton out(get_nsymbols(nfa), true, 0);
  StateSet start_ss;
  start_ss.insert(0);
  start_ss = get_epsilon_closure(start_ss, nfa);
  add_back(ssupv, start_ss);
  ssp2s[ssupv.back().get()] = add_state(out);
  int front = 0;
  while (front < int(ssupv.size())) {
    int state = front;
    StateSet& ss = *at(ssupv, front);
    ++front;
    for (int symbol = 0; symbol < get_nsymbols(nfa); ++symbol) {
      StateSet next_ss = Teuchos::step(ss, symbol, nfa);
      if (next_ss.empty()) continue;
      next_ss = get_epsilon_closure(next_ss, nfa);
      int next_state;
      StateSetPtr2State::iterator it = ssp2s.find(&next_ss);
      if (it == ssp2s.end()) {
        next_state = add_state(out);
        add_back(ssupv, next_ss);
        ssp2s[ssupv.back().get()] = next_state;
      } else {
        next_state = it->second;
      }
      add_transition(out, state, symbol, next_state);
    }
    int min_accepted = -1;
    for (StateSet::const_iterator it = ss.begin(); it != ss.end(); ++it) {
      int nfa_state = *it;
      int nfa_token = accepts(nfa, nfa_state);
      if (nfa_token == -1) continue;
      if (min_accepted == -1 || nfa_token < min_accepted) {
        min_accepted = nfa_token;
      }
    }
    if (min_accepted != -1) add_accept(out, state, min_accepted);
  }
  swap(result, out);
}

struct StateRowLess {
  Table<int> const& table;
  std::vector<int> const& accepted;
  bool operator()(int const& a, int const& b) const {
    int aa = at(accepted, a);
    int ab = at(accepted, b);
    if (aa != ab) return aa < ab;
    for (int symbol = 0, ncols = get_ncols(table); symbol < ncols; ++symbol) {
      int na = at(table, a, symbol);
      int nb = at(table, b, symbol);
      if (na != nb) return na < nb;
    }
    return false;
  }
  StateRowLess(Table<int> const& t, std::vector<int> const& a):
    table(t),accepted(a) {
  }
};

typedef std::map<int, int, StateRowLess> StateRow2SimpleState;

void simplify_once(FiniteAutomaton& result, FiniteAutomaton const& fa) {
  using std::swap;
  StateRow2SimpleState sr2ss(StateRowLess(fa.table, fa.accepted_tokens));
  int nsimple = 0;
  for (int state = 0; state < get_nstates(fa); ++state) {
    std::pair<StateRow2SimpleState::iterator, bool> res =
      sr2ss.insert(std::make_pair(state, nsimple));
    if (res.second) {
      ++nsimple;
    }
  }
  FiniteAutomaton out(get_nsymbols(fa), get_determinism(fa), nsimple);
  for (int simple = 0; simple < nsimple; ++simple) {
    add_state(out);
  }
  std::vector<bool> did_simple(size_t(nsimple), false);
  for (int state = 0; state < get_nstates(fa); ++state) {
    TEUCHOS_ASSERT(sr2ss.count(state));
    int simple = sr2ss[state];
    if (at(did_simple, simple)) continue;
    for (int symbol = 0; symbol < get_nsymbols_eps(fa); ++symbol) {
      int next_state = step(fa, state, symbol);
      if (next_state == -1) continue;
      TEUCHOS_ASSERT(sr2ss.count(next_state));
      int next_simple = sr2ss[next_state];
      add_transition(out, simple, symbol, next_simple);
    }
    int token = accepts(fa, state);
    if (token != -1) {
      add_accept(out, simple, token);
    }
    at(did_simple, simple) = true;
  }
  swap(result, out);
}

void simplify(FiniteAutomaton& result, FiniteAutomaton const& fa) {
  using std::swap;
  FiniteAutomaton prev;
  FiniteAutomaton next = fa;
  int nstates_next = get_nstates(next);
  int nstates;
  int i = 0;
  do {
    swap(prev, next);
    nstates = nstates_next;
    simplify_once(next, prev);
    ++i;
    nstates_next = get_nstates(next);
  } while (nstates_next < nstates);
  swap(result, next);
}

void add_char_transition(FiniteAutomaton& fa, int from_state, char at_char, int to_state) {
  add_transition(fa, from_state, get_symbol(at_char), to_state);
}

template <typename T, bool is_signed = std::numeric_limits<T>::is_signed>
struct IsSymbol;

template <typename T>
struct IsSymbol<T, true> {
  static bool eval(T c) {
    if (c < 0) return false;
    return 0 <= Teuchos::chartab[int(c)];
  }
};

template <typename T>
struct IsSymbol<T, false> {
  static bool eval(T c) {
    if (c >= TEUCHOS_CHARTAB_SIZE) return false;
    return 0 <= Teuchos::chartab[int(c)];
  }
};

bool is_symbol(char c) {
  return IsSymbol<char>::eval(c);
}

template <typename T, bool is_signed = std::numeric_limits<T>::is_signed>
struct GetSymbol;

template <typename T>
struct GetSymbol<T, true> {
  static int eval(T c) {
    TEUCHOS_ASSERT(0 <= c);
    int symbol = Teuchos::chartab[int(c)];
    TEUCHOS_ASSERT(0 <= symbol);
    return symbol;
  }
};

template <typename T>
struct GetSymbol<T, false> {
  static int eval(T c) {
    int symbol = Teuchos::chartab[int(c)];
    TEUCHOS_ASSERT(0 <= symbol);
    return symbol;
  }
};

int get_symbol(char c) {
  return GetSymbol<char>::eval(c);
}

char get_char(int symbol) {
  TEUCHOS_ASSERT(0 <= symbol);
  TEUCHOS_ASSERT(symbol < Teuchos::NCHARS);
  return inv_chartab[symbol];
}

void make_char_set_nfa(FiniteAutomaton& result, std::set<char> const& accepted, int token) {
  std::set<int> symbol_set;
  for (std::set<char>::const_iterator it = accepted.begin(); it != accepted.end(); ++it) {
    char c = *it;
    symbol_set.insert(get_symbol(c));
  }
  return make_set_nfa(result, Teuchos::NCHARS, symbol_set, token);
}

void make_char_range_nfa(FiniteAutomaton& result, char range_start, char range_end, int token) {
  return make_range_nfa(result, Teuchos::NCHARS, get_symbol(range_start), get_symbol(range_end), token);
}

void make_char_single_nfa(FiniteAutomaton& result, char symbol_char, int token) {
  return make_range_nfa(result, Teuchos::NCHARS, get_symbol(symbol_char), get_symbol(symbol_char), token);
}

void negate_set(std::set<char>& result, std::set<char> const& s) {
  using std::swap;
  std::set<char> out;
  for (int symbol = 0; symbol < NCHARS; ++symbol) {
    char c = inv_chartab[symbol];
    if (!s.count(c)) out.insert(c);
  }
  swap(result, out);
}

std::ostream& operator<<(std::ostream& os, FiniteAutomaton const& fa) {
  if (get_determinism(fa)) os << "dfa ";
  else os << "nfa ";
  os << get_nstates(fa) << " states " << get_nsymbols(fa) << " symbols\n";
  for (int state = 0; state < get_nstates(fa); ++state) {
    for (int symbol = 0; symbol < get_nsymbols(fa); ++symbol) {
      int next_state = step(fa, state, symbol);
      if (next_state != -1) os << "(" << state << ", " << symbol << ") -> " << next_state << '\n';
    }
    if (!get_determinism(fa)) {
      for (int symbol = get_epsilon0(fa); symbol <= get_epsilon1(fa); ++symbol) {
        int next_state = step(fa, state, symbol);
        if (next_state != -1) os << "(" << state << ", eps" << (symbol - get_epsilon0(fa)) << ") -> " << next_state << '\n';
      }
    }
    int token = accepts(fa, state);
    if (token != -1) os << state << " accepts " << token << '\n';
  }
  return os;
}

}  // end namespace Teuchos
