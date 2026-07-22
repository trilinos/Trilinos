// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_make_lalr1_parser.hpp"

#include <map>
#include <iostream>
#include <cstdlib>
#include <queue>
#include <algorithm>
#include <fstream>

#include "Teuchos_vector.hpp"
#include "Teuchos_Graph.hpp"
#include "Teuchos_stack.hpp"
#include "Teuchos_set.hpp"

namespace Teuchos {

/* The LALR(1) parser construction implemented here is based on
   David Pager's work:

  Pager, David.
  "The lane-tracing algorithm for constructing LR (k) parsers
   and ways of enhancing its efficiency."
  Information Sciences 12.1 (1977): 19-42.

  The identifiers used in this code are consistent with the terminology
  in that paper, except where we bring in FIRST set terminology, which
  Pager doesn't go into detail about. */

Config::Config(int p, int d):
  production(p),
  dot(d)
{
}

StateConfig::StateConfig(int s, int cis):
  state(s),
  config_in_state(cis)
{
}

void swap(StateInProgress& a, StateInProgress& b) {
  using std::swap;
  swap(a.configs, b.configs);
  swap(a.actions, b.actions);
}

// expand the grammar productions into marked productions
static Configs make_configs(Grammar const& g) {
  Configs configs;
  for (int i = 0; i < Teuchos::size(g.productions); ++i) {
    const Grammar::Production& production = at(g.productions, i);
    for (int j = 0; j <= Teuchos::size(production.rhs); ++j) {
      configs.push_back(Config(i,j));
    }
  }
  return configs;
}

static Graph get_left_hand_sides_to_start_configs(
    Configs const& cs, Grammar const& grammar) {
  Graph lhs2sc = make_graph_with_nnodes(grammar.nsymbols);
  for (int c_i = 0; c_i < Teuchos::size(cs); ++c_i) {
    const Config& c = at(cs, c_i);
    if (c.dot > 0) continue;
    int p_i = c.production;
    const Grammar::Production& p = at(grammar.productions, p_i);
    add_edge(lhs2sc, p.lhs, c_i);
  }
  return lhs2sc;
}

struct StateCompare {
  typedef StateInProgress const* Value;
  bool operator()(Value const& a, Value const& b) const {
    return a->configs < b->configs;
  }
};

typedef std::map<StateInProgress const*, int, StateCompare> StatePtr2StateIndex;

static void close(StateInProgress& state,
    Configs const& cs, Grammar const& grammar,
    Graph const& lhs2sc) {
  std::queue<int> config_q;
  std::set<int> config_set;
  for (std::vector<int>::const_iterator it = state.configs.begin();
       it != state.configs.end(); ++it) {
    int config_i = *it;
    config_q.push(config_i);
    TEUCHOS_ASSERT(!config_set.count(config_i));
    config_set.insert(config_i);
  }
  while (!config_q.empty()) {
    int config_i = config_q.front(); config_q.pop();
    const Config& config = at(cs, config_i);
    int prod_i = config.production;
    const Grammar::Production& prod = at(grammar.productions, prod_i);
    if (config.dot == Teuchos::size(prod.rhs)) continue;
    int symbol_after_dot = at(prod.rhs, config.dot);
    if (is_terminal(grammar, symbol_after_dot)) continue;
    const NodeEdges& edges = get_edges(lhs2sc, symbol_after_dot);
    for (NodeEdges::const_iterator it = edges.begin();
         it != edges.end(); ++it) {
      int sc = *it;
      if (!config_set.count(sc)) {
        config_set.insert(sc);
        config_q.push(sc);
      }
    }
  }
  state.configs.assign(config_set.begin(), config_set.end());
}

static void add_back(StatesInProgress& sips, StateInProgress& sip) {
  using std::swap;
  StateInProgressPtr ptr(new StateInProgress());
  swap(*ptr, sip);
  sips.push_back(ptr);
}

static void add_reduction_actions(StatesInProgress& states,
    Configs const& cs, Grammar const& grammar) {
  for (StatesInProgress::iterator it = states.begin();
       it != states.end(); ++it) {
    StateInProgressPtr& state_uptr = *it;
    StateInProgress& state = *state_uptr;
    for (std::vector<int>::const_iterator it2 = state.configs.begin();
         it2 != state.configs.end(); ++it2) {
      int config_i = *it2;
      const Config& config = at(cs, config_i);
      int prod_i = config.production;
      const Grammar::Production& prod = at(grammar.productions, prod_i);
      if (config.dot != Teuchos::size(prod.rhs)) continue;
      ActionInProgress reduction;
      reduction.action.kind = ACTION_REDUCE;
      reduction.action.production = config.production;
      state.actions.push_back(reduction);
    }
  }
}

static void set_lr0_contexts(
    StatesInProgress& states,
    Grammar const& grammar) {
  for (StatesInProgress::iterator it = states.begin();
       it != states.end(); ++it) {
    StateInProgressPtr& state_uptr = *it;
    StateInProgress& state = *state_uptr;
    for (StateInProgress::Actions::iterator it2 = state.actions.begin();
         it2 != state.actions.end(); ++it2) {
      ActionInProgress& action = *it2;
      if (action.action.kind != ACTION_REDUCE) continue;
      if (action.action.production == get_accept_production(grammar)) {
        action.context.insert(get_end_terminal(grammar));
      } else {
        for (int terminal = 0; terminal < grammar.nterminals; ++terminal) {
          action.context.insert(terminal);
        }
      }
    }
  }
}

static StatesInProgress make_lr0_parser(Configs const& cs, Grammar const& grammar,
    Graph const& lhs2sc) {
  StatesInProgress states;
  StatePtr2StateIndex state_ptrs2idxs;
  std::queue<int> state_q;
  { /* start state */
    StateInProgress start_state;
    int accept_nt = get_accept_nonterminal(grammar);
    /* there should only be one start configuration for the accept symbol */
    int start_accept_config = get_edges(lhs2sc, accept_nt).front();
    start_state.configs.push_back(start_accept_config);
    close(start_state, cs, grammar, lhs2sc);
    int start_state_i = Teuchos::size(states);
    state_q.push(start_state_i);
    add_back(states, start_state);
    state_ptrs2idxs[states.back().get()] = start_state_i;
  }
  while (!state_q.empty()) {
    int state_i = state_q.front(); state_q.pop();
    StateInProgress& state = *at(states, state_i);
    std::set<int> transition_symbols;
    for (std::vector<int>::const_iterator it = state.configs.begin();
         it != state.configs.end(); ++it) {
      int config_i = *it;
      const Config& config = at(cs, config_i);
      int prod_i = config.production;
      const Grammar::Production& prod = at(grammar.productions, prod_i);
      if (config.dot == Teuchos::size(prod.rhs)) continue;
      int symbol_after_dot = at(prod.rhs, config.dot);
      transition_symbols.insert(symbol_after_dot);
    }
    for (std::set<int>::const_iterator it = transition_symbols.begin();
         it != transition_symbols.end(); ++it) {
      int transition_symbol = *it;
      StateInProgress next_state;
      for (std::vector<int>::const_iterator it2 = state.configs.begin();
           it2 != state.configs.end(); ++it2) {
        int config_i = *it2;
        const Config& config = at(cs, config_i);
        int prod_i = config.production;
        const Grammar::Production& prod = at(grammar.productions, prod_i);
        if (config.dot == Teuchos::size(prod.rhs)) continue;
        int symbol_after_dot = at(prod.rhs, config.dot);
        if (symbol_after_dot != transition_symbol) continue;
        /* transition successor should just be the next index */
        int next_config_i = config_i + 1;
        next_state.configs.push_back(next_config_i);
      }
      close(next_state, cs, grammar, lhs2sc);
      StatePtr2StateIndex::iterator it2 = state_ptrs2idxs.find(&next_state);
      int next_state_i;
      if (it2 == state_ptrs2idxs.end()) {
        next_state_i = Teuchos::size(states);
        state_q.push(next_state_i);
        add_back(states, next_state);
        state_ptrs2idxs[states.back().get()] = next_state_i;
      } else {
        next_state_i = it2->second;
      }
      ActionInProgress transition;
      transition.action.kind = ACTION_SHIFT;
      transition.action.next_state = next_state_i;
      transition.context.insert(transition_symbol);
      state.actions.push_back(transition);
    }
  }
  add_reduction_actions(states, cs, grammar);
  set_lr0_contexts(states, grammar);
  return states;
}

static Graph get_productions_by_lhs(Grammar const& grammar) {
  int nsymbols = grammar.nsymbols;
  Graph lhs2prods = make_graph_with_nnodes(nsymbols);
  for (int prod_i = 0; prod_i < Teuchos::size(grammar.productions); ++prod_i) {
    const Grammar::Production& prod = at(grammar.productions, prod_i);
    add_edge(lhs2prods, prod.lhs, prod_i);
  }
  return lhs2prods;
}

/* compute a graph where symbols are graph nodes, and
   there exists an edge (A, B) if B appears in the RHS of
   any production in which A is the LHS */
static Graph get_symbol_graph(Grammar const& grammar, Graph const& lhs2prods) {
  int nsymbols = grammar.nsymbols;
  Graph out = make_graph_with_nnodes(nsymbols);
  for (int lhs = 0; lhs < nsymbols; ++lhs) {
    std::set<int> dependees;
    for (int i = 0; i < count_edges(lhs2prods, lhs); ++i) {
      int prod_i = at(lhs2prods, lhs, i);
      const Grammar::Production& prod = at(grammar.productions, prod_i);
      for (int j = 0; j < Teuchos::size(prod.rhs); ++j) {
        int rhs_symb = at(prod.rhs, j);
        dependees.insert(rhs_symb);
      }
    }
    at(out, lhs).assign(dependees.begin(), dependees.end());
  }
  return out;
}

/* the "FIRST" set, i.e. the set of 1-heads of non-null terminal descendants of
   some string.
   As suggested by Westley Weimer here:
   https://www.cs.virginia.edu/~weimer/2008-415/reading/FirstFollowLL.pdf
   we will also use the FIRST set for determining whether the string has
   a null terminal descendant, indicated by the prescence of a special
   FIRST_NULL symbol in the FIRST set */
enum { FIRST_NULL = -425 };
typedef std::set<int> FirstSet;

static void print_set(std::set<int> const& set, Grammar const& grammar) {
  std::cerr << "{";
  for (std::set<int>::const_iterator it = set.begin(); it != set.end(); ++it) {
    if (it != set.begin()) std::cerr << ", ";
    int symb = *it;
    if (symb == FIRST_NULL) std::cerr << "null";
    else {
      const std::string& symb_name = at(grammar.symbol_names, symb);
      if (symb_name == ",") std::cerr << "','";
      else std::cerr << symb_name;
    }
  }
  std::cerr << "}";
}

static FirstSet get_first_set_of_string(std::vector<int> const& string,
    std::vector<FirstSet> const& first_sets) {
  FirstSet out;
  /* walk the string, stop when any symbol is found that doesn't
     have a null terminal descendant */
  int i;
  for (i = 0; i < Teuchos::size(string); ++i) {
    int symbol = at(string, i);
    bool has_null = false;
    const FirstSet& first_set = at(first_sets, symbol);
    for (FirstSet::const_iterator it = first_set.begin();
         it != first_set.end(); ++it) {
      int first_symbol = *it;
      if (first_symbol == FIRST_NULL) has_null = true;
      else out.insert(first_symbol);
    }
    if (!has_null) break;
  }
  if (i == Teuchos::size(string)) out.insert(FIRST_NULL);
  return out;
}

struct Event {
  int added_symbol;
  int dependee;
  Event(int as, int d):
    added_symbol(as),
    dependee(d)
  {}
};

/* figure out the FIRST sets for each non-terminal in the grammar.
   I couldn't find a super-efficient way to do this, so here is a
   free-for-all event-driven implementation. */
static std::vector<FirstSet> compute_first_sets(Grammar const& grammar,
    bool verbose) {
  if (verbose) std::cerr << "computing FIRST sets...\n";
  std::queue<Event> event_q;
  int nsymbols = grammar.nsymbols;
  std::vector<FirstSet> first_sets = make_vector<FirstSet>(nsymbols);
  Graph lhs2prods = get_productions_by_lhs(grammar);
  for (int symbol = 0; symbol < nsymbols; ++symbol) {
    if (is_terminal(grammar, symbol)) {
      event_q.push(Event(symbol, symbol));
    } else {
      for (int i = 0; i < count_edges(lhs2prods, symbol); ++i) {
        int prod_i = at(lhs2prods, symbol, i);
        const Grammar::Production& prod = at(grammar.productions, prod_i);
        if (prod.rhs.empty()) {
          event_q.push(Event(FIRST_NULL, symbol));
          break;
        }
      }
    }
  }
  Graph dependers2dependees = get_symbol_graph(grammar, lhs2prods);
  Graph dependees2dependers = make_transpose(dependers2dependees);
  while (!event_q.empty()) {
    Event event = event_q.front(); event_q.pop();
    int added_symb = event.added_symbol;
    int dependee = event.dependee;
    FirstSet& dependee_firsts = at(first_sets, dependee);
    /* hopefully we don't get too many duplicate events piled up... */
    if (dependee_firsts.count(added_symb)) continue;
    dependee_firsts.insert(added_symb);
    for (int i = 0; i < count_edges(dependees2dependers, dependee); ++i) {
      int depender = at(dependees2dependers, dependee, i);
      TEUCHOS_ASSERT(is_nonterminal(grammar, depender));
      const FirstSet& depender_firsts = at(first_sets, depender);
      for (int j = 0; j < count_edges(lhs2prods, depender); ++j) {
        int prod_i = at(lhs2prods, depender, j);
        const Grammar::Production& prod = at(grammar.productions, prod_i);
        FirstSet rhs_first_set = get_first_set_of_string(prod.rhs, first_sets);
        for (FirstSet::iterator it = rhs_first_set.begin();
             it != rhs_first_set.end(); ++it) {
          int rhs_first_symb = *it;
          if (!depender_firsts.count(rhs_first_symb)) {
            event_q.push(Event(rhs_first_symb, depender));
          }
        }
      }
    }
  }
  if (verbose) {
    for (int symb = 0; symb < nsymbols; ++symb) {
      const std::string& symb_name = at(grammar.symbol_names, symb);
      std::cerr << "FIRST(" << symb_name << ") = {";
      const FirstSet& c = at(first_sets, symb);
      for (FirstSet::const_iterator it = c.begin(); it != c.end(); ++it) {
        if (it != c.begin()) std::cerr << ", ";
        int first_symb = *it;
        if (first_symb == FIRST_NULL) {
          std::cerr << "null";
        } else {
          const std::string& first_name = at(grammar.symbol_names, first_symb);
          std::cerr << first_name;
        }
      }
      std::cerr << "}\n";
    }
    std::cerr << '\n';
  }
  return first_sets;
}

StateConfigs form_state_configs(StatesInProgress const& states) {
  StateConfigs out;
  for (int i = 0; i < Teuchos::size(states); ++i) {
    StateInProgress& state = *at(states, i);
    for (int j = 0; j < Teuchos::size(state.configs); ++j) {
      out.push_back(StateConfig(i, j));
    }
  }
  return out;
}

Graph form_states_to_state_configs(StateConfigs const& scs,
    StatesInProgress const& states) {
  Graph out = make_graph_with_nnodes(Teuchos::size(states));
  for (int i = 0; i < Teuchos::size(scs); ++i) {
    const StateConfig& sc = at(scs, i);
    at(out, sc.state).push_back(i);
  }
  return out;
}

static std::string escape_dot(std::string const& s) {
  std::string out;
  for (std::size_t i = 0; i < s.size(); ++i) {
    char c = s[i];
    if (c == '\\' || c == '|' || c == '\"' || c == '<' || c == '>' || c == '{' || c == '}') {
      out.push_back('\\');
      out.push_back(c);
    } else if (c == '.') {
      out = "'";
      out += s;
      out += "'";
      return out;
    } else {
      out.push_back(c);
    }
  }
  return out;
}

void print_graphviz(
    std::string const& filepath,
    ParserInProgress const& pip,
    bool /* verbose */,
    std::ostream& os
    ) {
  const StatesInProgress& sips = pip.states;
  const Configs& cs = pip.configs;
  GrammarPtr grammar = pip.grammar;
  const Graph& states2scs = pip.states2state_configs;
  os << "writing GraphViz file \"" << filepath << "\"\n";
  os << "process with:\n";
  os << "  dot -Tpdf -o \"" << filepath << ".pdf\" \"" << filepath << "\"\n";
  std::ofstream file(filepath.c_str());
  TEUCHOS_ASSERT(file.is_open());
  file << "digraph {\n";
  file << "graph [\n";
  file << "rankdir = \"LR\"\n";
  file << "]\n";
  for (int s_i = 0; s_i < Teuchos::size(sips); ++s_i) {
    const StateInProgress& state = *at(sips, s_i);
    file << s_i << " [\n";
    file << "label = \"";
    file << "State " << s_i << "\\l";
    for (int cis_i = 0; cis_i < Teuchos::size(state.configs); ++cis_i) {
      int c_i = at(state.configs, cis_i);
      const Config& config = at(cs, c_i);
      const Grammar::Production& prod = at(grammar->productions, config.production);
      int sc_i = at(states2scs, s_i, cis_i);
      file << sc_i << ": ";
      const std::string& lhs_name = at(grammar->symbol_names, prod.lhs);
      file << escape_dot(lhs_name) << " ::= ";
      for (int rhs_i = 0; rhs_i <= Teuchos::size(prod.rhs); ++rhs_i) {
        if (rhs_i == config.dot) file << " .";
        if (rhs_i < Teuchos::size(prod.rhs)) {
          int rhs_symb = at(prod.rhs, rhs_i);
          const std::string& rhs_symb_name = at(grammar->symbol_names, rhs_symb);
          file << " " << escape_dot(rhs_symb_name);
        }
      }
      if (config.dot == Teuchos::size(prod.rhs)) {
        file << ", \\{";
        bool found = false;
        for (int a_i = 0; a_i < Teuchos::size(state.actions); ++a_i) {
          const ActionInProgress& action = at(state.actions, a_i);
          if (action.action.kind == ACTION_REDUCE &&
              action.action.production == config.production) {
            found = true;
            const Context& ac = action.context;
            for (Context::const_iterator it = ac.begin(); it != ac.end(); ++it) {
              if (it != ac.begin()) file << ", ";
              int symb = *it;
              const std::string& symb_name = at(grammar->symbol_names, symb);
              file << escape_dot(symb_name);
            }
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error,
            "BUG: missing reduce action in state " << s_i << " !!!\n");
        file << "\\}";
      }
      file << "\\l";
    }
    file << "\"\n";
    file << "shape = \"record\"\n";
    file << "]\n";
    for (int a_i = 0; a_i < Teuchos::size(state.actions); ++a_i) {
      const ActionInProgress& action = at(state.actions, a_i);
      if (action.action.kind == ACTION_SHIFT) {
        int symb = *(action.context.begin());
        const std::string& symb_name = at(grammar->symbol_names, symb);
        int next = action.action.next_state;
        file << s_i << " -> " << next << " [\n";
        file << "label = \"" << escape_dot(symb_name) << "\"\n";
        file << "]\n";
      }
    }
  }
  file << "}\n";
}

static Graph make_immediate_predecessor_graph(
    StateConfigs const& scs,
    StatesInProgress const& states,
    Graph const& states2scs,
    Configs const& cs,
    GrammarPtr grammar) {
  Graph out = make_graph_with_nnodes(Teuchos::size(scs));
  for (int s_i = 0; s_i < Teuchos::size(states); ++s_i) {
    const StateInProgress& state = *at(states, s_i);
    for (int cis_i = 0; cis_i < Teuchos::size(state.configs); ++cis_i) {
      int config_i = at(state.configs, cis_i);
      const Config& config = at(cs, config_i);
      const Grammar::Production& prod = at(grammar->productions, config.production);
      int dot = config.dot;
      if (dot == Teuchos::size(prod.rhs)) continue;
      int s = at(prod.rhs, dot);
      if (is_terminal(*grammar, s)) continue;
      for (int cis_j = 0; cis_j < Teuchos::size(state.configs); ++cis_j) {
        int config_j = at(state.configs, cis_j);
        const Config& config2 = at(cs, config_j);
        const Grammar::Production& prod2 = at(grammar->productions, config2.production);
        if (prod2.lhs == s) {
          int sc_i = at(states2scs, s_i, cis_i);
          int sc_j = at(states2scs, s_i, cis_j);
          add_edge(out, sc_j, sc_i);
        }
      }
    }
  }
  return out;
}

static Graph find_transition_predecessors(
    StateConfigs const& scs,
    StatesInProgress const& states,
    Graph const& states2scs,
    Configs const& cs,
    GrammarPtr grammar) {
  Graph out = make_graph_with_nnodes(Teuchos::size(scs));
  for (int state_i = 0; state_i < Teuchos::size(states); ++state_i) {
    const StateInProgress& state = *at(states, state_i);
    for (int a_i = 0; a_i < Teuchos::size(state.actions); ++a_i) {
      const ActionInProgress& action = at(state.actions, a_i);
      if (action.action.kind != ACTION_SHIFT) continue;
      TEUCHOS_ASSERT(action.context.size() == 1);
      int symbol = *(action.context.begin());
      int state_j = action.action.next_state;
      const StateInProgress& state2 = *at(states, state_j);
      for (int cis_i = 0; cis_i < Teuchos::size(state.configs); ++cis_i) {
        int config_i = at(state.configs, cis_i);
        const Config& config = at(cs, config_i);
        for (int cis_j = 0; cis_j < Teuchos::size(state2.configs); ++cis_j) {
          int config_j = at(state2.configs, cis_j);
          const Config& config2 = at(cs, config_j);
          if (config.production == config2.production &&
              config.dot + 1 == config2.dot) {
            const Grammar::Production& prod = at(grammar->productions, config.production);
            int rhs_symbol = at(prod.rhs, config.dot);
            if (rhs_symbol == symbol) {
              int sc_i = at(states2scs, state_i, cis_i);
              int sc_j = at(states2scs, state_j, cis_j);
              add_edge(out, sc_j, sc_i);
            }
          }
        }
      }
    }
  }
  return out;
}

static Graph make_originator_graph(
    StateConfigs const& scs,
    StatesInProgress const& states,
    Graph const& states2scs,
    Configs const& cs,
    GrammarPtr grammar) {
  Graph out = make_graph_with_nnodes(Teuchos::size(scs));
  Graph ipg = make_immediate_predecessor_graph(
      scs, states, states2scs, cs, grammar);
  Graph tpg = find_transition_predecessors(
      scs, states, states2scs, cs, grammar);
  for (int sc_i = 0; sc_i < Teuchos::size(scs); ++sc_i) {
    std::set<int> originators;
    /* breadth-first search through the Transition
       Precessor Graph (tpg), followed by a single hop
       along the Immediate Predecessor Graph (ipg) */
    std::queue<int> tpq;
    std::set<int> tps;
    tpq.push(sc_i);
    tps.insert(sc_i);
    while (!tpq.empty()) {
      int tpp = tpq.front(); tpq.pop();
      for (int i = 0; i < count_edges(tpg, tpp); ++i) {
        int tpc = at(tpg, tpp, i);
        if (tps.count(tpc)) continue;
        tpq.push(tpc);
        tps.insert(tpc);
      }
      for (int i = 0; i < count_edges(ipg, tpp); ++i) {
        int ip_i = at(ipg, tpp, i);
        originators.insert(ip_i);
      }
    }
    at(out, sc_i).assign(originators.begin(), originators.end());
  }
  return out;
}

static std::vector<int> get_follow_string(
    int sc_addr,
    StateConfigs const& scs,
    StatesInProgress const& states,
    Configs const& cs,
    GrammarPtr grammar) {
  const StateConfig& sc = at(scs, sc_addr);
  const StateInProgress& state = *at(states, sc.state);
  int config_i = at(state.configs, sc.config_in_state);
  const Config& config = at(cs, config_i);
  const Grammar::Production& prod = at(grammar->productions, config.production);
  int out_size = Teuchos::size(prod.rhs) - (config.dot + 1);
  std::vector<int> out;
  /* out_size can be negative */
  if (out_size < 1) return out;
  reserve(out, out_size);
  for (int i = config.dot + 1; i < Teuchos::size(prod.rhs); ++i) {
    out.push_back(at(prod.rhs, i));
  }
  return out;
}

static void print_string(std::vector<int> const& str, GrammarPtr grammar) {
  std::cerr << "\"";
  for (int i = 0; i < Teuchos::size(str); ++i) {
    int symb = at(str, i);
    const std::string& symb_name = at(grammar->symbol_names, symb);
    std::cerr << symb_name;
  }
  std::cerr << "\"";
}

static bool has_non_null_terminal_descendant(FirstSet const& first_set) {
  if (first_set.empty()) return false;
  if (first_set.size() > 1) return true;
  return *(first_set.begin()) != FIRST_NULL;
}

static Context get_contexts(FirstSet first_set) {
  FirstSet::iterator it = first_set.find(FIRST_NULL);
  if (it != first_set.end()) first_set.erase(it);
  return first_set;
}

enum { MARKER = -433 };
enum { ZERO = -100 }; // actual zero is a valid index for us

static void print_stack(std::vector<int> const& stack) {
  for (int i = 0; i < Teuchos::size(stack); ++i) {
    int symb = at(stack, i);
    if (symb == MARKER) std::cerr << " M";
    else if (symb == ZERO) std::cerr << " Z";
    else std::cerr << " " << symb;
  }
  std::cerr << '\n';
}

static void move_markers(
    std::vector<int>& lane,
    std::vector<int>& in_lane,
    int zeta_prime_addr,
    int zeta_pointer,
    bool tests_failed
    ) {
  int loc_of_zeta_prime = at(in_lane, zeta_prime_addr);
  TEUCHOS_ASSERT(loc_of_zeta_prime != -1);
  int r = 0;
  for (int i = loc_of_zeta_prime + 1; i < zeta_pointer; ++i) {
    if (at(lane, i) == MARKER) {
      ++r;
      at(lane, i) = ZERO;
    }
  }
  int top_addr = -1;
  if (tests_failed) {
    top_addr = lane.back();
    lane.resize(lane.size() - 1); // pop
  }
  for (int i = 0; i < r; ++i) lane.push_back(MARKER);
  if (tests_failed) {
    lane.push_back(top_addr);
    at(in_lane, top_addr) = Teuchos::size(lane) - 1;
  }
}

typedef std::vector<Context> Contexts;

static void context_adding_routine(
    std::vector<int> const& lane,
    int zeta_pointer,
    Context& contexts_generated,
    Contexts& contexts,
    bool verbose,
    GrammarPtr grammar
    ) {
  if (verbose) {
    std::cerr << "  CONTEXT ADDING ROUTINE\n";
    std::cerr << "  LANE:";
    print_stack(lane);
    std::cerr << "  $\\zeta$-POINTER = " << zeta_pointer << '\n';
  }
  for (int r = zeta_pointer; r >= 0 && (!contexts_generated.empty()); --r) {
    int v_r = at(lane, r);
    if (verbose) std::cerr << "    r = " << r << ", $v_r$ = ";
    if (v_r < 0) {
      if (verbose) {
        if (v_r == MARKER) std::cerr << "marker\n";
        else if (v_r == ZERO) std::cerr << "zero\n";
      }
      continue;
    }
    int tau_r_addr = v_r;
    if (verbose) {
      std::cerr << "$\\tau_r$ = " << tau_r_addr << '\n';
      std::cerr << "    CONTEXTS_GENERATED = ";
      print_set(contexts_generated, *grammar);
      std::cerr << "\n    CONTEXTS_$\\tau_r$ = ";
      print_set(at(contexts, tau_r_addr), *grammar);
      std::cerr << "\n    CONTEXTS_GENERATED <- CONTEXTS_GENERATED - CONTEXTS_$\\tau_r$";
    }
    subtract_from(contexts_generated, at(contexts, tau_r_addr));
    if (verbose) {
      std::cerr << "\n    CONTEXTS_GENERATED = ";
      print_set(contexts_generated, *grammar);
      std::cerr << "\n    CONTEXTS_$\\tau_r$ <- CONTEXTS_$\\tau_r$ U CONTEXTS_GENERATED";
    }
    unite_with(at(contexts, tau_r_addr), contexts_generated);
    if (verbose) {
      std::cerr << "\n    CONTEXTS_$\\tau_r$ = ";
      print_set(at(contexts, tau_r_addr), *grammar);
      std::cerr << "\n";
    }
  }
}

static void deal_with_tests_failed(
    int& num_originators_failed,
    int& first_originator_failed,
    int zeta_prime_addr,
    bool& tests_failed,
    std::vector<int>& lane,
    std::vector<int>& in_lane,
    int zeta_addr,
    std::vector<int>& stack,
    bool verbose
    ) {
  if (verbose) std::cerr << "  Dealing with test failures\n";
  if (num_originators_failed == 0) {
    if (verbose) std::cerr << "    " << zeta_prime_addr << " is the first originator of " << zeta_addr << " to fail the tests\n";
    first_originator_failed = zeta_prime_addr;
    if (verbose) std::cerr << "    pushing " << zeta_prime_addr << " onto LANE:\n    ";
    lane.push_back(zeta_prime_addr);
    if (verbose) print_stack(lane);
    at(in_lane, zeta_prime_addr) = Teuchos::size(lane) - 1;
    if (verbose) std::cerr << "    IN_LANE(" << zeta_prime_addr << ") <- ON\n";
    tests_failed = true;
    if (verbose) std::cerr << "    TESTS_FAILED <- ON\n";
  } else if (num_originators_failed == 1) {
    if (verbose) std::cerr << "    " << zeta_prime_addr << " is the second originator of " << zeta_addr << " to fail the tests\n";
    int zeta_double_prime_addr = first_originator_failed;
    if (verbose) std::cerr << "    the first was " << zeta_double_prime_addr << '\n';
    TEUCHOS_ASSERT(at(lane, Teuchos::size(lane) - 1) == zeta_double_prime_addr);
    TEUCHOS_ASSERT(at(lane, Teuchos::size(lane) - 2) == zeta_addr);
    if (verbose) std::cerr << "    pop LANE, push {marker, " << zeta_double_prime_addr << "} onto it:\n    ";
    lane.resize(lane.size() - 1);
    lane.push_back(MARKER);
    lane.push_back(zeta_double_prime_addr);
    if (verbose) print_stack(lane);
    if (verbose) std::cerr << "    push {marker, " << zeta_prime_addr << "} onto STACK:\n    ";
    stack.push_back(MARKER);
    stack.push_back(zeta_prime_addr);
    if (verbose) print_stack(stack);
  } else {
    if (verbose) std::cerr << "    " << zeta_prime_addr << " is the third or later originator of " << zeta_addr << " to fail the tests\n";
    if (verbose) std::cerr << "    pushing " << zeta_prime_addr << " onto STACK:\n    ";
    stack.push_back(zeta_prime_addr);
    if (verbose) print_stack(stack);
  }
  ++num_originators_failed;
}

static void heuristic_propagation_of_context_sets(
    int tau_addr,
    Contexts& contexts,
    std::vector<bool>& complete,
    StateConfigs const& scs,
    StatesInProgress const& states,
    Graph const& states2scs,
    Configs const& cs,
    GrammarPtr grammar
    ) {
  const StateConfig& tau = at(scs, tau_addr);
  const StateInProgress& state = *at(states, tau.state);
  int config_i = at(state.configs, tau.config_in_state);
  const Config& config = at(cs, config_i);
  if (config.dot != 0) return;
  const Grammar::Production& prod = at(grammar->productions, config.production);
  for (int cis_j = 0; cis_j < Teuchos::size(state.configs); ++cis_j) {
    int config_j = at(state.configs, cis_j);
    if (config_j == config_i) continue;
    const Config& config2 = at(cs, config_j);
    if (config2.dot != 0) continue;
    const Grammar::Production& prod2 = at(grammar->productions, config2.production);
    if (prod.lhs != prod2.lhs) continue;
    int tau_prime_addr = at(states2scs, tau.state, cis_j);
    at(contexts, tau_prime_addr) = at(contexts, tau_addr);
    at(complete, tau_prime_addr) = true;
  }
}

/* Here it is! The magical algorithm described by a flowchart in
   Figure 7 of David Pager's paper. */
static void compute_context_set(
    int zeta_j_addr,
    Contexts& contexts,
    std::vector<bool>& complete,
    StateConfigs const& scs,
    Graph const& originator_graph,
    StatesInProgress const& states,
    Graph const& states2scs,
    Configs const& cs,
    std::vector<FirstSet> const& first_sets,
    GrammarPtr grammar,
    bool verbose,
    ParserInProgress const& pip
    ) {
  if (verbose) std::cerr << "Computing context set for $\\zeta_j$ = " << zeta_j_addr << "...\n";
  if (verbose) std::cerr << "BEGIN PROGRAM\n";
  if (at(complete, zeta_j_addr)) {
    if (verbose) std::cerr << zeta_j_addr << " was already complete!\nEND PROGRAM\n\n";
    return;
  }
  std::vector<int> stack;
  // need random access, inner insert, which std::stack doesn't provide
  std::vector<int> lane;
  std::vector<int> in_lane = make_vector<int>(Teuchos::size(scs), -1);
  lane.push_back(zeta_j_addr);
  at(in_lane, zeta_j_addr) = Teuchos::size(lane) - 1;
  bool tests_failed = false;
  Context contexts_generated;
  if (verbose) {
    std::cerr << "Initial LANE:";
    print_stack(lane);
  }
  while (true) {
    TEUCHOS_ASSERT(!lane.empty());
    int zeta_addr = lane.back();
    if (verbose) {
      std::cerr << "Top of LANE is $\\zeta$ = " << zeta_addr << '\n';
    }
    int zeta_pointer = Teuchos::size(lane) - 1;
    if (verbose) std::cerr << "$\\zeta$-POINTER <- " << zeta_pointer << '\n';
    int num_originators_failed = 0;
    int first_originator_failed = -1;
    if (verbose) std::cerr << "DO_LOOP:\n";
    /* DO_LOOP */
    for (int i = 0; i < count_edges(originator_graph, zeta_addr); ++i) {
      int zeta_prime_addr = at(originator_graph, zeta_addr, i);
      if (verbose) {
        std::cerr << "Next originator of $\\zeta$ = " << zeta_addr << " is $\\zeta'$ = " << zeta_prime_addr << '\n';
      }
      std::vector<int> gamma = get_follow_string(zeta_prime_addr, scs, states, cs, grammar);
      if (verbose) {
        std::cerr << "  FOLLOW string of $\\zeta'$ = " << zeta_prime_addr << " is ";
        print_string(gamma, grammar);
        std::cerr << '\n';
      }
      FirstSet gamma_first = get_first_set_of_string(gamma, first_sets);
      if (verbose) {
        std::cerr << "  FIRST set of ";
        print_string(gamma, grammar);
        std::cerr << " is ";
        print_set(gamma_first, *grammar);
        std::cerr << "\n";
      }
      if (has_non_null_terminal_descendant(gamma_first)) { // test A
        if (verbose) {
          std::cerr << "  ";
          print_string(gamma, grammar);
          std::cerr << " has a non-null terminal descendant\n";
        }
        contexts_generated = get_contexts(gamma_first);
        if (verbose) {
          std::cerr << "  CONTEXTS_GENERATED = ";
          print_set(contexts_generated, *grammar);
          std::cerr << " = 1-heads of non-null descendants of ";
          print_string(gamma, grammar);
          std::cerr << '\n';
        }
        if (gamma_first.count(FIRST_NULL)) {
          if (verbose) {
            std::cerr << "  ";
            print_string(gamma, grammar);
            std::cerr << " has a null terminal descendant\n";
          }
          if (at(complete, zeta_prime_addr)) {
            unite_with(contexts_generated, at(contexts, zeta_prime_addr));
            context_adding_routine(lane, zeta_pointer, contexts_generated, contexts,
                verbose, grammar);
          } else if (-1 == at(in_lane, zeta_prime_addr)) {
            context_adding_routine(lane, zeta_pointer, contexts_generated, contexts,
                verbose, grammar);
            /* TRACE_FURTHER */
            deal_with_tests_failed(num_originators_failed, first_originator_failed,
                zeta_prime_addr, tests_failed, lane, in_lane, zeta_addr, stack,
                verbose);
          } else {
            print_graphviz("ambiguous.dot", pip, true, std::cerr);
            std::stringstream ss;
            ss << "error: grammar is ambiguous.\n";
            ss << "zeta_j is " << zeta_j_addr << ", zeta is " << zeta_addr << ", and zeta prime is " << zeta_prime_addr << '\n';
            throw ParserBuildFail(ss.str());
          }
        } else {
          context_adding_routine(lane, zeta_pointer, contexts_generated, contexts,
              verbose, grammar);
        }
      } else {
        if (verbose) {
          std::cerr << "  ";
          print_string(gamma, grammar);
          std::cerr << " does not have a non-null terminal descendant\n";
        }
        if (at(complete, zeta_prime_addr)) { // test B
          if (verbose) std::cerr << "  COMPLETE(" << zeta_prime_addr << ") is ON\n";
          contexts_generated = at(contexts, zeta_prime_addr);
          context_adding_routine(lane, zeta_pointer, contexts_generated, contexts,
              verbose, grammar);
        } else {
          if (verbose) std::cerr << "  COMPLETE(" << zeta_prime_addr << ") is OFF\n";
          if (-1 != at(in_lane, zeta_prime_addr)) { // test C
            if (verbose) std::cerr << "  IN_LANE(" << zeta_prime_addr << ") is ON\n";
            move_markers(lane, in_lane, zeta_prime_addr, zeta_pointer, tests_failed);
            contexts_generated = at(contexts, zeta_prime_addr);
            context_adding_routine(lane, zeta_pointer, contexts_generated, contexts,
                verbose, grammar);
          } else {
            if (verbose) std::cerr << "  IN_LANE(" << zeta_prime_addr << ") is OFF\n";
            deal_with_tests_failed(num_originators_failed, first_originator_failed,
                zeta_prime_addr, tests_failed, lane, in_lane, zeta_addr, stack,
                verbose);
          }
        }
      }
    } /* END DO_LOOP */
    if (verbose) std::cerr << "END DO_LOOP\n";
    if (tests_failed) {
      if (verbose) {
        std::cerr << "  TESTS_FAILED was on, turning it off and going to next configuration\n";
      }
      tests_failed = false;
      continue;
    }
    bool keep_lane_popping = true;
    if (verbose) std::cerr << "  Start LANE popping\n";
    while (keep_lane_popping) { // LANE popping loop
      TEUCHOS_ASSERT(!lane.empty());
      if (verbose) {
        std::cerr << "  LANE:";
        print_stack(lane);
      }
      if (at(lane, Teuchos::size(lane) - 1) == MARKER) {
        if (verbose) std::cerr << "  Top of LANE is a marker\n";
        if (verbose) std::cerr << "  Start STACK popping\n";
        while (true) { // STACK popping loop
          TEUCHOS_ASSERT(!stack.empty());
          if (verbose) {
            std::cerr << "    STACK:";
            print_stack(stack);
            std::cerr << "    LANE:";
            print_stack(lane);
          }
          if (stack.back() == MARKER) {
            if (verbose) std::cerr << "    Top of STACK is a marker, pop STACK and LANE\n";
            resize(stack, Teuchos::size(stack) - 1);
            resize(lane, Teuchos::size(lane) - 1);
            break; // out of STACK popping, back into LANE popping
          } else if (at(complete, stack.back())) {
            if (verbose) std::cerr << "    Top of STACK is has COMPLETE flag, pop STACK\n";
            resize(stack, Teuchos::size(stack) - 1);
            // back into STACK popping
          } else {
            int addr = stack.back();
            if (verbose) std::cerr << "    Top of STACK is " << addr << ", pop STACK\n";
            resize(stack, Teuchos::size(stack) - 1);
            if (verbose) std::cerr << "    Push " << addr << " onto LANE\n";
            lane.push_back(addr);
            if (verbose) std::cerr << "    IN_LANE(" << addr << ") <- ON\n";
            at(in_lane, addr) = Teuchos::size(lane) - 1;
            keep_lane_popping = false;
            break; // out of STACK and LANE popping, into top-level loop
          } // end STACK top checks
        } // end STACK popping loop
      } else if (at(lane, Teuchos::size(lane) - 1) == ZERO) {
        if (verbose) std::cerr << "  Top of LANE is a zero\n";
        if (verbose) std::cerr << "  Pop LANE\n";
        resize(lane, Teuchos::size(lane) - 1); // pop LANE
        // back to top of LANE popping loop
      } else { // top of LANE neither marker nor zero
        int tau_addr = lane.back();
        if (verbose) std::cerr << "  Top of LANE is $\\tau$ = " << tau_addr << "\n";
        at(in_lane, tau_addr) = -1;
        if (verbose) std::cerr << "  IN_LANE(" << tau_addr << ") <- OFF\n";
        at(complete, tau_addr) = true;
        if (verbose) std::cerr << "  COMPLETE(" << tau_addr << ") <- ON\n";
        if (verbose) std::cerr << "  HEURISTIC PROPAGATION OF CONTEXT SETS\n";
        heuristic_propagation_of_context_sets(
            tau_addr, contexts, complete,
            scs, states, states2scs, cs, grammar);
        if (Teuchos::size(lane) == 1 && at(lane, 0) == zeta_j_addr) {
          if (verbose) std::cerr << "END PROGRAM\n\n";
          return;
        }
        if (verbose) std::cerr << "  Pop LANE\n";
        resize(lane, Teuchos::size(lane) - 1); // pop LANE
        // back to top of LANE popping loop
      } // end top of LANE checks
    } // end LANE popping loop
  } // end top-level while(1) loop
}

static std::vector<bool> determine_adequate_states(
    StatesInProgress const& states,
    GrammarPtr grammar,
    bool verbose,
    std::ostream& os) {
  std::vector<bool> out = make_vector<bool>(Teuchos::size(states));
  for (int s_i = 0; s_i < Teuchos::size(states); ++s_i) {
    const StateInProgress& state = *at(states, s_i);
    bool state_is_adequate = true;
    for (int a_i = 0; a_i < Teuchos::size(state.actions); ++a_i) {
      const ActionInProgress& action = at(state.actions, a_i);
      if (action.action.kind == ACTION_SHIFT &&
          is_nonterminal(*grammar, *(action.context.begin()))) {
        continue;
      }
      for (int a_j = a_i + 1; a_j < Teuchos::size(state.actions); ++a_j) {
        const ActionInProgress& action2 = at(state.actions, a_j);
        if (action2.action.kind == ACTION_SHIFT &&
            is_nonterminal(*grammar, *(action2.context.begin()))) {
          continue;
        }
        if (intersects(action2.context, action.context)) {
          if (verbose) {
            const ActionInProgress* ap1 = &action;
            const ActionInProgress* ap2 = &action2;
            if (ap1->action.kind == ACTION_SHIFT) {
              std::swap(ap1, ap2);
            }
            TEUCHOS_ASSERT(ap1->action.kind == ACTION_REDUCE);
            os << "shift-reduce conflict in state " << s_i << ":\n";
            os << "reduce ";
            const Grammar::Production& prod = at(grammar->productions, ap1->action.production);
            const std::string& lhs_name = at(grammar->symbol_names, prod.lhs);
            os << lhs_name << " ::=";
            for (int rhs_i = 0; rhs_i < Teuchos::size(prod.rhs); ++rhs_i) {
              int rhs_symb = at(prod.rhs, rhs_i);
              const std::string& rhs_symb_name = at(grammar->symbol_names, rhs_symb);
              os << " " << rhs_symb_name;
            }
            int shift_symb = *(ap2->context.begin());
            const std::string& shift_name = at(grammar->symbol_names, shift_symb);
            os << "\nshift " << shift_name << '\n';
          }
          state_is_adequate = false;
          break;
        }
      }
      if (!state_is_adequate) break;
    }
    at(out, s_i) = state_is_adequate;
  }
  if (verbose) os << '\n';
  return out;
}

ParserInProgress draft_lalr1_parser(GrammarPtr grammar, bool verbose) {
  ParserInProgress out;
  Configs& cs = out.configs;
  StatesInProgress& states = out.states;
  StateConfigs& scs = out.state_configs;
  Graph& states2scs = out.states2state_configs;
  out.grammar = grammar;
  cs = make_configs(*grammar);
  Graph lhs2cs = get_left_hand_sides_to_start_configs(cs, *grammar);
  if (verbose) std::cerr << "Building LR(0) parser\n";
  states = make_lr0_parser(cs, *grammar, lhs2cs);
  scs = form_state_configs(states);
  states2scs = form_states_to_state_configs(scs, states);
  if (verbose) print_graphviz("lr0.dot", out, true, std::cerr);
  if (verbose) std::cerr << "Checking adequacy of LR(0) machine\n";
  std::vector<bool> adequate = determine_adequate_states(states, grammar, verbose,
      std::cerr);
  if (*(std::min_element(adequate.begin(), adequate.end()))) {
    if (verbose) std::cerr << "The grammar is LR(0)!\n";
    return out;
  }
  std::vector<bool> complete = make_vector<bool>(Teuchos::size(scs), false);
  std::vector<Context> contexts = make_vector<Context>(Teuchos::size(scs));
  int accept_prod_i = get_accept_production(*grammar);
  /* initialize the accepting state-configs as described in
     footnote 8 at the bottom of page 37 */
  for (int sc_i = 0; sc_i < Teuchos::size(scs); ++sc_i) {
    StateConfig& sc = at(scs, sc_i);
    StateInProgress& state = *at(states, sc.state);
    int config_i = at(state.configs, sc.config_in_state);
    Config& config = at(cs, config_i);
    if (config.production == accept_prod_i) {
      at(complete, sc_i) = true;
      at(contexts, sc_i).insert(get_end_terminal(*grammar));
    }
  }
  Graph og = make_originator_graph(scs, states, states2scs, cs, grammar);
  if (verbose) std::cerr << "Originator Graph:\n";
  if (verbose) std::cerr << og << '\n';
  std::vector<FirstSet> first_sets = compute_first_sets(*grammar, verbose);
  /* compute context sets for all state-configs associated with reduction
     actions that are part of an inadequate state */
  for (int s_i = 0; s_i < Teuchos::size(states); ++s_i) {
    if (at(adequate, s_i)) continue;
    StateInProgress& state = *at(states, s_i);
    for (int cis_i = 0; cis_i < Teuchos::size(state.configs); ++cis_i) {
      int config_i = at(state.configs, cis_i);
      const Config& config = at(cs, config_i);
      const Grammar::Production& prod = at(grammar->productions, config.production);
      if (config.dot != Teuchos::size(prod.rhs)) continue;
      int zeta_j_addr = at(states2scs, s_i, cis_i);
      compute_context_set(zeta_j_addr, contexts, complete, scs,
          og, states, states2scs, cs, first_sets, grammar, verbose, out);
    }
  }
  /* update the context sets for all reduction state-configs
     which are marked complete, even if they aren't in inadequate states */
  for (int s_i = 0; s_i < Teuchos::size(states); ++s_i) {
    StateInProgress& state = *at(states, s_i);
    for (int cis_i = 0; cis_i < Teuchos::size(state.configs); ++cis_i) {
      int sc_i = at(states2scs, s_i, cis_i);
      if (!at(complete, sc_i)) continue;
      int config_i = at(state.configs, cis_i);
      Config& config = at(cs, config_i);
      const Grammar::Production& prod = at(grammar->productions, config.production);
      if (config.dot != Teuchos::size(prod.rhs)) continue;
      for (int a_i = 0; a_i < Teuchos::size(state.actions); ++a_i) {
        ActionInProgress& action = at(state.actions, a_i);
        if (action.action.kind == ACTION_REDUCE &&
            action.action.production == config.production) {
          action.context = at(contexts, sc_i);
        }
      }
    }
  }
  if (verbose) std::cerr << "Checking adequacy of LALR(1) machine\n";
  adequate = determine_adequate_states(states, grammar, verbose, std::cerr);
  if (!(*(std::min_element(adequate.begin(), adequate.end())))) {
    std::stringstream ss;
    ss << "error: The grammar is not LALR(1).\n";
    determine_adequate_states(states, grammar, true, ss);
    print_graphviz("error.dot", out, true, ss);
    std::string s = ss.str();
    throw ParserBuildFail(s);
  }
  if (verbose) std::cerr << "The grammar is LALR(1)!\n";
  if (verbose) print_graphviz("lalr1.dot", out, true, std::cerr);
  return out;
}

Parser accept_parser(ParserInProgress const& pip) {
  const StatesInProgress& sips = pip.states;
  GrammarPtr grammar = pip.grammar;
  Parser out(grammar, Teuchos::size(sips));
  for (int s_i = 0; s_i < Teuchos::size(sips); ++s_i) {
    add_state(out);
  }
  for (int s_i = 0; s_i < Teuchos::size(sips); ++s_i) {
    const StateInProgress& sip = *at(sips, s_i);
    for (int a_i = 0; a_i < Teuchos::size(sip.actions); ++a_i) {
      const ActionInProgress& action = at(sip.actions, a_i);
      if (action.action.kind == ACTION_SHIFT &&
          is_nonterminal(*grammar, *(action.context.begin()))) {
        int nt = as_nonterminal(*grammar, *(action.context.begin()));
        add_nonterminal_action(out, s_i, nt, action.action.next_state);
      } else {
        for (Context::const_iterator it = action.context.begin();
             it != action.context.end(); ++it) {
          int terminal = *it;
          TEUCHOS_ASSERT(is_terminal(*grammar, terminal));
          add_terminal_action(out, s_i, terminal, action.action);
        }
      }
    }
  }
  return out;
}

ParserBuildFail::ParserBuildFail(const std::string& msg):
  std::invalid_argument(msg) {
}

Parser make_lalr1_parser(GrammarPtr grammar, bool verbose) {
  ParserInProgress pip = draft_lalr1_parser(grammar, verbose);
  return accept_parser(pip);
}

}
