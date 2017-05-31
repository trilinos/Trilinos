#include "Teuchos_build_parser.hpp"

#include <map>
#include <iostream>
#include <cstdlib>
#include <queue>
#include <algorithm>
#include <fstream>

#include "Teuchos_vector.hpp"
#include "Teuchos_graph.hpp"
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

// expand the grammar productions into marked productions
static Configs make_configs(Grammar const& g) {
  Configs configs;
  for (int i = 0; i < size(g.productions); ++i) {
    auto& production = at(g.productions, i);
    for (int j = 0; j <= size(production.rhs); ++j) {
      configs.push_back({i,j});
    }
  }
  return configs;
}

static Graph get_left_hand_sides_to_start_configs(
    Configs const& cs, Grammar const& grammar) {
  auto lhs2sc = make_graph_with_nnodes(grammar.nsymbols);
  for (int c_i = 0; c_i < size(cs); ++c_i) {
    auto& c = at(cs, c_i);
    if (c.dot > 0) continue;
    auto p_i = c.production;
    auto& p = at(grammar.productions, p_i);
    add_edge(lhs2sc, p.lhs, c_i);
  }
  return lhs2sc;
}

struct StateCompare {
  using Value = StateInProgress const*;
  bool operator()(Value const& a, Value const& b) const {
    return a->configs < b->configs;
  }
};

using StatePtr2StateIndex = std::map<StateInProgress const*, int, StateCompare>;

static void close(StateInProgress& state,
    Configs const& cs, Grammar const& grammar,
    Graph const& lhs2sc) {
  std::queue<int> config_q;
  std::set<int> config_set;
  for (auto config_i : state.configs) {
    config_q.push(config_i);
    assert(!config_set.count(config_i));
    config_set.insert(config_i);
  }
  while (!config_q.empty()) {
    auto config_i = config_q.front(); config_q.pop();
    auto& config = at(cs, config_i);
    auto prod_i = config.production;
    auto& prod = at(grammar.productions, prod_i);
    if (config.dot == size(prod.rhs)) continue;
    auto symbol_after_dot = at(prod.rhs, config.dot);
    if (is_terminal(grammar, symbol_after_dot)) continue;
    for (auto sc : get_edges(lhs2sc, symbol_after_dot)) {
      if (!config_set.count(sc)) {
        config_set.insert(sc);
        config_q.push(sc);
      }
    }
  }
  state.configs.assign(config_set.begin(), config_set.end());
}

static void emplace_back(StatesInProgress& sips, StateInProgress& sip) {
  sips.push_back(std::unique_ptr<StateInProgress>(new StateInProgress(std::move(sip))));
}

static void add_reduction_actions(StatesInProgress& states,
    Configs const& cs, Grammar const& grammar) {
  for (auto& state_uptr : states) {
    auto& state = *state_uptr;
    for (auto config_i : state.configs) {
      auto& config = at(cs, config_i);
      auto prod_i = config.production;
      auto& prod = at(grammar.productions, prod_i);
      if (config.dot != size(prod.rhs)) continue;
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
  for (auto& state_uptr : states) {
    auto& state = *state_uptr;
    for (auto& action : state.actions) {
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

static StatesInProgress build_lr0_parser(Configs const& cs, Grammar const& grammar,
    Graph const& lhs2sc) {
  StatesInProgress states;
  StatePtr2StateIndex state_ptrs2idxs;
  std::queue<int> state_q;
  { /* start state */
    StateInProgress start_state;
    auto accept_nt = get_accept_nonterminal(grammar);
    /* there should only be one start configuration for the accept symbol */
    auto start_accept_config = get_edges(lhs2sc, accept_nt).front();
    start_state.configs.push_back(start_accept_config);
    close(start_state, cs, grammar, lhs2sc);
    auto start_state_i = size(states);
    state_q.push(start_state_i);
    emplace_back(states, start_state);
    state_ptrs2idxs[states.back().get()] = start_state_i;
  }
  while (!state_q.empty()) {
    auto state_i = state_q.front(); state_q.pop();
    auto& state = *at(states, state_i);
    std::set<int> transition_symbols;
    for (auto config_i : state.configs) {
      auto& config = at(cs, config_i);
      auto prod_i = config.production;
      auto& prod = at(grammar.productions, prod_i);
      if (config.dot == size(prod.rhs)) continue;
      auto symbol_after_dot = at(prod.rhs, config.dot);
      transition_symbols.insert(symbol_after_dot);
    }
    for (auto transition_symbol : transition_symbols) {
      StateInProgress next_state;
      for (auto config_i : state.configs) {
        auto& config = at(cs, config_i);
        auto prod_i = config.production;
        auto& prod = at(grammar.productions, prod_i);
        if (config.dot == size(prod.rhs)) continue;
        auto symbol_after_dot = at(prod.rhs, config.dot);
        if (symbol_after_dot != transition_symbol) continue;
        /* transition successor should just be the next index */
        auto next_config_i = config_i + 1;
        next_state.configs.push_back(next_config_i);
      }
      close(next_state, cs, grammar, lhs2sc);
      auto it = state_ptrs2idxs.find(&next_state);
      int next_state_i;
      if (it == state_ptrs2idxs.end()) {
        next_state_i = size(states);
        state_q.push(next_state_i);
        emplace_back(states, next_state);
        state_ptrs2idxs[states.back().get()] = next_state_i;
      } else {
        next_state_i = it->second;
      }
      ActionInProgress transition;
      transition.action.kind = ACTION_SHIFT;
      transition.action.next_state = next_state_i;
      transition.context.insert(transition_symbol);
      state.actions.emplace_back(std::move(transition));
    }
  }
  add_reduction_actions(states, cs, grammar);
  set_lr0_contexts(states, grammar);
  return states;
}

static Graph get_productions_by_lhs(Grammar const& grammar) {
  auto nsymbols = grammar.nsymbols;
  auto lhs2prods = make_graph_with_nnodes(nsymbols);
  for (int prod_i = 0; prod_i < size(grammar.productions); ++prod_i) {
    auto& prod = at(grammar.productions, prod_i);
    add_edge(lhs2prods, prod.lhs, prod_i);
  }
  return lhs2prods;
}

/* compute a graph where symbols are graph nodes, and
   there exists an edge (A, B) if B appears in the RHS of
   any production in which A is the LHS */
static Graph get_symbol_graph(Grammar const& grammar, Graph const& lhs2prods) {
  auto nsymbols = grammar.nsymbols;
  auto out = make_graph_with_nnodes(nsymbols);
  for (int lhs = 0; lhs < nsymbols; ++lhs) {
    std::set<int> dependees;
    auto& lhs_prods = get_edges(lhs2prods, lhs);
    for (auto prod_i : lhs_prods) {
      auto& prod = at(grammar.productions, prod_i);
      for (auto rhs_symb : prod.rhs) {
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
using FirstSet = std::set<int>;

static void print_set(std::set<int> const& set, Grammar const& grammar) {
  std::cerr << "{";
  for (auto it = set.begin(); it != set.end(); ++it) {
    if (it != set.begin()) std::cerr << ", ";
    auto symb = *it;
    if (symb == FIRST_NULL) std::cerr << "null";
    else {
      auto& symb_name = at(grammar.symbol_names, symb);
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
  for (i = 0; i < size(string); ++i) {
    auto symbol = at(string, i);
    bool has_null = false;
    for (auto first_symbol : at(first_sets, symbol)) {
      if (first_symbol == FIRST_NULL) has_null = true;
      else out.insert(first_symbol);
    }
    if (!has_null) break;
  }
  if (i == size(string)) out.insert(FIRST_NULL);
  return out;
}

/* figure out the FIRST sets for each non-terminal in the grammar.
   I couldn't find a super-efficient way to do this, so here is a
   free-for-all event-driven implementation. */
static std::vector<FirstSet> compute_first_sets(Grammar const& grammar,
    bool verbose) {
  if (verbose) std::cerr << "computing FIRST sets...\n";
  struct Event {
    int added_symbol;
    int dependee;
  };
  std::queue<Event> event_q;
  auto nsymbols = grammar.nsymbols;
  auto first_sets = make_vector<FirstSet>(nsymbols);
  auto lhs2prods = get_productions_by_lhs(grammar);
  for (int symbol = 0; symbol < nsymbols; ++symbol) {
    if (is_terminal(grammar, symbol)) {
      event_q.push({symbol, symbol});
    } else {
      auto& lhs_prods = get_edges(lhs2prods, symbol);
      for (auto prod_i : lhs_prods) {
        auto& prod = at(grammar.productions, prod_i);
        if (prod.rhs.empty()) {
          event_q.push({FIRST_NULL, symbol});
          break;
        }
      }
    }
  }
  auto dependers2dependees = get_symbol_graph(grammar, lhs2prods);
  auto dependees2dependers = make_transpose(dependers2dependees);
  while (!event_q.empty()) {
    auto event = event_q.front(); event_q.pop();
    auto added_symb = event.added_symbol;
    auto dependee = event.dependee;
    auto& dependee_firsts = at(first_sets, dependee);
    /* hopefully we don't get too many duplicate events piled up... */
    if (dependee_firsts.count(added_symb)) continue;
    dependee_firsts.insert(added_symb);
    auto& dependers = get_edges(dependees2dependers, dependee);
    for (auto depender : dependers) {
      assert(is_nonterminal(grammar, depender));
      auto& prods = get_edges(lhs2prods, depender);
      auto const& depender_firsts = at(first_sets, depender);
      for (auto prod_i : prods) {
        auto& prod = at(grammar.productions, prod_i);
        auto rhs_first_set = get_first_set_of_string(prod.rhs, first_sets);
        for (auto rhs_first_symb : rhs_first_set) {
          if (!depender_firsts.count(rhs_first_symb)) {
            event_q.push({rhs_first_symb, depender});
          }
        }
      }
    }
  }
  if (verbose) {
    for (int symb = 0; symb < nsymbols; ++symb) {
      auto& symb_name = at(grammar.symbol_names, symb);
      std::cerr << "FIRST(" << symb_name << ") = {";
      auto& c = at(first_sets, symb);
      for (auto it = c.begin(); it != c.end(); ++it) {
        if (it != c.begin()) std::cerr << ", ";
        auto first_symb = *it;
        if (first_symb == FIRST_NULL) {
          std::cerr << "null";
        } else {
          auto& first_name = at(grammar.symbol_names, first_symb);
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
  for (int i = 0; i < size(states); ++i) {
    auto& state = *at(states, i);
    for (int j = 0; j < size(state.configs); ++j) out.push_back({i, j});
  }
  return out;
}

Graph form_states_to_state_configs(StateConfigs const& scs,
    StatesInProgress const& states) {
  auto out = make_graph_with_nnodes(size(states));
  for (int i = 0; i < size(scs); ++i) {
    auto& sc = at(scs, i);
    at(out, sc.state).push_back(i);
  }
  return out;
}

static std::string escape_dot(std::string const& s) {
  std::string out;
  for (auto c : s) {
    if (c == '\\' || c == '|' || c == '\"' || c == '<' || c == '>') {
      out.push_back('\\');
      out.push_back(c);
    } else if (c == '.') {
      out += "\'.\'";
    } else {
      out.push_back(c);
    }
  }
  return out;
}

void print_dot(
    std::string const& filepath,
    ParserInProgress const& pip
    ) {
  auto& sips = pip.states;
  auto& cs = pip.configs;
  auto& grammar = pip.grammar;
  auto& states2scs = pip.states2state_configs;
  std::cerr << "writing " << filepath << "\n\n";
  std::ofstream file(filepath.c_str());
  assert(file.is_open());
  file << "digraph {\n";
  file << "graph [\n";
  file << "rankdir = \"LR\"\n";
  file << "]\n";
  for (int s_i = 0; s_i < size(sips); ++s_i) {
    auto& state = *at(sips, s_i);
    file << s_i << " [\n";
    file << "label = \"";
    file << "State " << s_i << "\\l";
    for (int cis_i = 0; cis_i < size(state.configs); ++cis_i) {
      auto c_i = at(state.configs, cis_i);
      auto& config = at(cs, c_i);
      auto& prod = at(grammar->productions, config.production);
      auto sc_i = at(states2scs, s_i, cis_i);
      file << sc_i << ": ";
      auto lhs_name = at(grammar->symbol_names, prod.lhs);
      file << escape_dot(lhs_name) << " ::= ";
      for (int rhs_i = 0; rhs_i <= size(prod.rhs); ++rhs_i) {
        if (rhs_i == config.dot) file << " .";
        if (rhs_i < size(prod.rhs)) {
          auto rhs_symb = at(prod.rhs, rhs_i);
          auto rhs_symb_name = at(grammar->symbol_names, rhs_symb);
          file << " " << escape_dot(rhs_symb_name);
        }
      }
      if (config.dot == size(prod.rhs)) {
        file << ", \\{";
        bool found = false;
        for (auto& action : state.actions) {
          if (action.action.kind == ACTION_REDUCE &&
              action.action.production == config.production) {
            found = true;
            auto& ac = action.context;
            for (auto it = ac.begin(); it != ac.end(); ++it) {
              if (it != ac.begin()) file << ", ";
              auto symb = *it;
              auto& symb_name = at(grammar->symbol_names, symb);
              file << escape_dot(symb_name);
            }
          }
        }
        if (!found) {
          std::cerr << "BUG: missing reduce action in state " << s_i << " !!!\n";
          abort();
        }
        file << "\\}";
      }
      file << "\\l";
    }
    file << "\"\n";
    file << "shape = \"record\"\n";
    file << "]\n";
    for (auto& action : state.actions) {
      if (action.action.kind == ACTION_SHIFT) {
        auto symb = *(action.context.begin());
        auto symb_name = at(grammar->symbol_names, symb);
        auto next = action.action.next_state;
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
  auto out = make_graph_with_nnodes(size(scs));
  for (int s_i = 0; s_i < size(states); ++s_i) {
    auto& state = *at(states, s_i);
    for (int cis_i = 0; cis_i < size(state.configs); ++cis_i) {
      auto config_i = at(state.configs, cis_i);
      auto& config = at(cs, config_i);
      auto& prod = at(grammar->productions, config.production);
      auto dot = config.dot;
      if (dot == size(prod.rhs)) continue;
      auto s = at(prod.rhs, dot);
      if (is_terminal(*grammar, s)) continue;
      for (int cis_j = 0; cis_j < size(state.configs); ++cis_j) {
        auto config_j = at(state.configs, cis_j);
        auto& config2 = at(cs, config_j);
        auto& prod2 = at(grammar->productions, config2.production);
        if (prod2.lhs == s) {
          auto sc_i = at(states2scs, s_i, cis_i);
          auto sc_j = at(states2scs, s_i, cis_j);
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
  auto out = make_graph_with_nnodes(size(scs));
  for (int state_i = 0; state_i < size(states); ++state_i) {
    auto& state = *at(states, state_i);
    for (auto& action : state.actions) {
      if (action.action.kind != ACTION_SHIFT) continue;
      assert(action.context.size() == 1);
      auto symbol = *(action.context.begin());
      auto state_j = action.action.next_state;
      auto& state2 = *at(states, state_j);
      for (int cis_i = 0; cis_i < size(state.configs); ++cis_i) {
        auto config_i = at(state.configs, cis_i);
        auto& config = at(cs, config_i);
        for (int cis_j = 0; cis_j < size(state2.configs); ++cis_j) {
          auto config_j = at(state2.configs, cis_j);
          auto& config2 = at(cs, config_j);
          if (config.production == config2.production &&
              config.dot + 1 == config2.dot) {
            auto& prod = at(grammar->productions, config.production);
            auto rhs_symbol = at(prod.rhs, config.dot);
            if (rhs_symbol == symbol) {
              auto sc_i = at(states2scs, state_i, cis_i);
              auto sc_j = at(states2scs, state_j, cis_j);
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
  auto out = make_graph_with_nnodes(size(scs));
  auto ipg = make_immediate_predecessor_graph(
      scs, states, states2scs, cs, grammar);
  auto tpg = find_transition_predecessors(
      scs, states, states2scs, cs, grammar);
  for (auto sc_i = 0; sc_i < size(scs); ++sc_i) {
    std::set<int> originators;
    /* breadth-first search through the transition
       precessor graph, followed by a single hop
       along the immediate predecessor graph */
    std::queue<int> tpq;
    std::set<int> tps;
    tpq.push(sc_i);
    tps.insert(sc_i);
    while (!tpq.empty()) {
      auto tpp = tpq.front(); tpq.pop();
      for (auto tpc : get_edges(tpg, tpp)) {
        if (tps.count(tpc)) continue;
        tpq.push(tpc);
        tps.insert(tpc);
      }
      for (auto ip_i : get_edges(ipg, tpp)) {
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
  auto& sc = at(scs, sc_addr);
  auto& state = *at(states, sc.state);
  auto config_i = at(state.configs, sc.config_in_state);
  auto& config = at(cs, config_i);
  auto& prod = at(grammar->productions, config.production);
  auto out_size = size(prod.rhs) - (config.dot + 1);
  std::vector<int> out;
  /* out_size can be negative */
  if (out_size < 1) return out;
  reserve(out, out_size);
  for (auto i = config.dot + 1; i < size(prod.rhs); ++i) {
    out.push_back(at(prod.rhs, i));
  }
  return out;
}

static void print_string(std::vector<int> const& str, GrammarPtr grammar) {
  std::cerr << "\"";
  for (auto symb : str) {
    auto& symb_name = at(grammar->symbol_names, symb);
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
  auto it = first_set.find(FIRST_NULL);
  if (it != first_set.end()) first_set.erase(it);
  return first_set;
}

enum { MARKER = -433 };
enum { ZERO = -100 }; // actual zero is a valid index for us

static void print_stack(std::vector<int> const& stack) {
  for (auto& symb : stack) {
    if (symb == MARKER) std::cerr << " M";
    else if (symb == ZERO) std::cerr << " Z";
    else std::cerr << " " << symb;
  }
  std::cerr << '\n';
}

static void move_markers(
    std::vector<int>& lane,
    int zeta_prime_addr,
    int zeta_pointer,
    bool tests_failed
    ) {
/* TODO: change in_lane to contain the index of that config in the lane,
   not just a boolean. this would save us the search here: */
  auto it = std::find_if(lane.begin(), lane.end(),
      [=](int item) { return item == zeta_prime_addr; });
  assert(it != lane.end());
  auto loc_of_zeta_prime = int(it - lane.begin());
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
  if (tests_failed) lane.push_back(top_addr);
}

using Contexts = std::vector<Context>;

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
    auto v_r = at(lane, r);
    if (verbose) std::cerr << "    r = " << r << ", $v_r$ = ";
    if (v_r < 0) {
      if (verbose) {
        if (v_r == MARKER) std::cerr << "marker\n";
        else if (v_r == ZERO) std::cerr << "zero\n";
      }
      continue;
    }
    auto tau_r_addr = v_r;
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
    std::vector<bool>& in_lane,
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
    at(in_lane, zeta_prime_addr) = true;
    if (verbose) std::cerr << "    IN_LANE(" << zeta_prime_addr << ") <- ON\n";
    tests_failed = true;
    if (verbose) std::cerr << "    TESTS_FAILED <- ON\n";
  } else if (num_originators_failed == 1) {
    if (verbose) std::cerr << "    " << zeta_prime_addr << " is the second originator of " << zeta_addr << " to fail the tests\n";
    auto zeta_double_prime_addr = first_originator_failed;
    if (verbose) std::cerr << "    the first was " << zeta_double_prime_addr << '\n';
    assert(at(lane, size(lane) - 1) == zeta_double_prime_addr);
    assert(at(lane, size(lane) - 2) == zeta_addr);
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
  auto& tau = at(scs, tau_addr);
  auto& state = *at(states, tau.state);
  auto config_i = at(state.configs, tau.config_in_state);
  auto& config = at(cs, config_i);
  if (config.dot != 0) return;
  auto& prod = at(grammar->productions, config.production);
  for (int cis_j = 0; cis_j < size(state.configs); ++cis_j) {
    auto config_j = at(state.configs, cis_j);
    if (config_j == config_i) continue;
    auto& config2 = at(cs, config_j);
    if (config2.dot != 0) continue;
    auto& prod2 = at(grammar->productions, config2.production);
    if (prod.lhs != prod2.lhs) continue;
    auto tau_prime_addr = at(states2scs, tau.state, cis_j);
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
    bool verbose
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
  auto in_lane = make_vector<bool>(size(scs), false);
  lane.push_back(zeta_j_addr);
  at(in_lane, zeta_j_addr) = true;
  bool tests_failed = false;
  Context contexts_generated;
  if (verbose) {
    std::cerr << "Initial LANE:";
    print_stack(lane);
  }
  while (true) {
    assert(!lane.empty());
    auto zeta_addr = lane.back();
    if (verbose) {
      std::cerr << "Top of LANE is $\\zeta$ = " << zeta_addr << '\n';
    }
    auto zeta_pointer = size(lane) - 1;
    if (verbose) std::cerr << "$\\zeta$-POINTER <- " << zeta_pointer << '\n';
    int num_originators_failed = 0;
    int first_originator_failed = -1;
    if (verbose) std::cerr << "DO_LOOP:\n";
    /* DO_LOOP */
    for (auto zeta_prime_addr : get_edges(originator_graph, zeta_addr)) {
      if (verbose) {
        std::cerr << "Next originator of $\\zeta$ = " << zeta_addr << " is $\\zeta'$ = " << zeta_prime_addr << '\n';
      }
      auto gamma = get_follow_string(zeta_prime_addr, scs, states, cs, grammar);
      if (verbose) {
        std::cerr << "  FOLLOW string of $\\zeta'$ = " << zeta_prime_addr << " is ";
        print_string(gamma, grammar);
        std::cerr << '\n';
      }
      auto gamma_first = get_first_set_of_string(gamma, first_sets);
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
          } else if (!at(in_lane, zeta_prime_addr)) {
            context_adding_routine(lane, zeta_pointer, contexts_generated, contexts,
                verbose, grammar);
            /* TRACE_FURTHER */
            deal_with_tests_failed(num_originators_failed, first_originator_failed,
                zeta_prime_addr, tests_failed, lane, in_lane, zeta_addr, stack,
                verbose);
          } else {
            std::cerr << "ERROR: grammar is ambiguous.\n";
            abort();
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
          if (at(in_lane, zeta_prime_addr)) { // test C
            if (verbose) std::cerr << "  IN_LANE(" << zeta_prime_addr << ") is ON\n";
            move_markers(lane, zeta_prime_addr, zeta_pointer, tests_failed);
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
      assert(!lane.empty());
      if (verbose) {
        std::cerr << "  LANE:";
        print_stack(lane);
      }
      if (at(lane, size(lane) - 1) == MARKER) {
        if (verbose) std::cerr << "  Top of LANE is a marker\n";
        if (verbose) std::cerr << "  Start STACK popping\n";
        while (true) { // STACK popping loop
          assert(!stack.empty());
          if (verbose) {
            std::cerr << "    STACK:";
            print_stack(stack);
            std::cerr << "    LANE:";
            print_stack(lane);
          }
          if (stack.back() == MARKER) {
            if (verbose) std::cerr << "    Top of STACK is a marker, pop STACK and LANE\n";
            resize(stack, size(stack) - 1);
            resize(lane, size(lane) - 1);
            break; // out of STACK popping, back into LANE popping
          } else if (at(complete, stack.back())) {
            if (verbose) std::cerr << "    Top of STACK is has COMPLETE flag, pop STACK\n";
            resize(stack, size(stack) - 1);
            // back into STACK popping
          } else {
            auto addr = stack.back();
            if (verbose) std::cerr << "    Top of STACK is " << addr << ", pop STACK\n";
            resize(stack, size(stack) - 1);
            if (verbose) std::cerr << "    Push " << addr << " onto LANE\n";
            lane.push_back(addr);
            if (verbose) std::cerr << "    IN_LANE(" << addr << ") <- ON\n";
            at(in_lane, addr) = true;
            keep_lane_popping = false;
            break; // out of STACK and LANE popping, into top-level loop
          } // end STACK top checks
        } // end STACK popping loop
      } else if (at(lane, size(lane) - 1) == ZERO) {
        if (verbose) std::cerr << "  Top of LANE is a zero\n";
        if (verbose) std::cerr << "  Pop LANE\n";
        resize(lane, size(lane) - 1); // pop LANE
        // back to top of LANE popping loop
      } else { // top of LANE neither marker nor zero
        auto tau_addr = lane.back();
        if (verbose) std::cerr << "  Top of LANE is $\\tau$ = " << tau_addr << "\n";
        at(in_lane, tau_addr) = false;
        if (verbose) std::cerr << "  IN_LANE(" << tau_addr << ") <- OFF\n";
        at(complete, tau_addr) = true;
        if (verbose) std::cerr << "  COMPLETE(" << tau_addr << ") <- ON\n";
        if (verbose) std::cerr << "  HEURISTIC PROPAGATION OF CONTEXT SETS\n";
        heuristic_propagation_of_context_sets(
            tau_addr, contexts, complete,
            scs, states, states2scs, cs, grammar);
        if (size(lane) == 1 && at(lane, 0) == zeta_j_addr) {
          if (verbose) std::cerr << "END PROGRAM\n\n";
          return;
        }
        if (verbose) std::cerr << "  Pop LANE\n";
        resize(lane, size(lane) - 1); // pop LANE
        // back to top of LANE popping loop
      } // end top of LANE checks
    } // end LANE popping loop
  } // end top-level while(1) loop
}

static std::vector<bool> determine_adequate_states(
    StatesInProgress const& states,
    GrammarPtr grammar,
    bool verbose) {
  auto out = make_vector<bool>(size(states));
  for (int s_i = 0; s_i < size(states); ++s_i) {
    auto& state = *at(states, s_i);
    bool state_is_adequate = true;
    for (int a_i = 0; a_i < size(state.actions); ++a_i) {
      auto& action = at(state.actions, a_i);
      if (action.action.kind == ACTION_SHIFT &&
          is_nonterminal(*grammar, *(action.context.begin()))) {
        continue;
      }
      for (int a_j = a_i + 1; a_j < size(state.actions); ++a_j) {
        auto& action2 = at(state.actions, a_j);
        if (action2.action.kind == ACTION_SHIFT &&
            is_nonterminal(*grammar, *(action2.context.begin()))) {
          continue;
        }
        if (intersects(action2.context, action.context)) {
          if (verbose) {
            auto* ap1 = &action;
            auto* ap2 = &action2;
            if (ap1->action.kind == ACTION_SHIFT) {
              std::swap(ap1, ap2);
            }
            assert(ap1->action.kind == ACTION_REDUCE);
            std::cerr << "shift-reduce conflict in state " << s_i << ":\n";
            std::cerr << "reduce ";
            auto& prod = at(grammar->productions, ap1->action.production);
            auto& lhs_name = at(grammar->symbol_names, prod.lhs);
            std::cerr << lhs_name << " ::=";
            for (auto rhs_symb : prod.rhs) {
              auto& rhs_symb_name = at(grammar->symbol_names, rhs_symb);
              std::cerr << " " << rhs_symb_name;
            }
            auto shift_symb = *(ap2->context.begin());
            auto shift_name = at(grammar->symbol_names, shift_symb);
            std::cerr << "\nshift " << shift_name << '\n';
          }
          state_is_adequate = false;
          break;
        }
      }
      if (!state_is_adequate) break;
    }
    at(out, s_i) = state_is_adequate;
  }
  if (verbose) std::cerr << '\n';
  return out;
}

ParserInProgress build_lalr1_parser(GrammarPtr grammar, bool verbose) {
  ParserInProgress out;
  auto& cs = out.configs;
  auto& states = out.states;
  auto& scs = out.state_configs;
  auto& states2scs = out.states2state_configs;
  out.grammar = grammar;
  cs = make_configs(*grammar);
  auto lhs2cs = get_left_hand_sides_to_start_configs(cs, *grammar);
  if (verbose) std::cerr << "Building LR(0) parser\n";
  states = build_lr0_parser(cs, *grammar, lhs2cs);
  scs = form_state_configs(states);
  states2scs = form_states_to_state_configs(scs, states);
  if (verbose) print_dot("lr0.dot", out);
  if (verbose) std::cerr << "Checking adequacy of LR(0) machine\n";
  auto adequate = determine_adequate_states(states, grammar, verbose);
  if (*(std::min_element(adequate.begin(), adequate.end()))) {
    if (verbose) std::cerr << "The grammar is LR(0)!\n";
    return out;
  }
  auto complete = make_vector<bool>(size(scs), false);
  auto contexts = make_vector<Context>(size(scs));
  auto accept_prod_i = get_accept_production(*grammar);
  /* initialize the accepting state-configs as described in
     footnote 8 at the bottom of page 37 */
  for (int sc_i = 0; sc_i < size(scs); ++sc_i) {
    auto& sc = at(scs, sc_i);
    auto& state = *at(states, sc.state);
    auto config_i = at(state.configs, sc.config_in_state);
    auto& config = at(cs, config_i);
    if (config.production == accept_prod_i) {
      at(complete, sc_i) = true;
      at(contexts, sc_i).insert(get_end_terminal(*grammar));
    }
  }
  auto og = make_originator_graph(scs, states, states2scs, cs, grammar);
  if (verbose) std::cerr << "Originator Graph:\n";
  if (verbose) std::cerr << og << '\n';
  auto first_sets = compute_first_sets(*grammar, verbose);
  /* compute context sets for all state-configs associated with reduction
     actions that are part of an inadequate state */
  for (int s_i = 0; s_i < size(states); ++s_i) {
    if (at(adequate, s_i)) continue;
    auto& state = *at(states, s_i);
    for (int cis_i = 0; cis_i < size(state.configs); ++cis_i) {
      auto config_i = at(state.configs, cis_i);
      auto& config = at(cs, config_i);
      auto& prod = at(grammar->productions, config.production);
      if (config.dot != size(prod.rhs)) continue;
      auto zeta_j_addr = at(states2scs, s_i, cis_i);
      compute_context_set(zeta_j_addr, contexts, complete, scs,
          og, states, states2scs, cs, first_sets, grammar, verbose);
    }
  }
  /* update the context sets for all reduction state-configs
     which are marked complete, even if they aren't in inadequate states */
  for (int s_i = 0; s_i < size(states); ++s_i) {
    auto& state = *at(states, s_i);
    for (int cis_i = 0; cis_i < size(state.configs); ++cis_i) {
      auto sc_i = at(states2scs, s_i, cis_i);
      if (!at(complete, sc_i)) continue;
      auto config_i = at(state.configs, cis_i);
      auto& config = at(cs, config_i);
      auto& prod = at(grammar->productions, config.production);
      if (config.dot != size(prod.rhs)) continue;
      for (auto& action : state.actions) {
        if (action.action.kind == ACTION_REDUCE &&
            action.action.production == config.production) {
          action.context = at(contexts, sc_i);
        }
      }
    }
  }
  if (verbose) std::cerr << "Checking adequacy of LALR(1) machine\n";
  adequate = determine_adequate_states(states, grammar, verbose);
  if (!(*(std::min_element(adequate.begin(), adequate.end())))) {
    std::cerr << "ERROR: The grammar is not LALR(1).\n";
    determine_adequate_states(states, grammar, true);
    print_dot("error.dot", out);
    abort();
  }
  if (verbose) std::cerr << "The grammar is LALR(1)!\n";
  if (verbose) print_dot("lalr1.dot", out);
  return out;
}

Parser accept_parser(ParserInProgress const& pip) {
  auto& sips = pip.states;
  auto& grammar = pip.grammar;
  auto out = Parser(grammar, size(sips));
  for (int s_i = 0; s_i < size(sips); ++s_i) {
    add_state(out);
  }
  for (int s_i = 0; s_i < size(sips); ++s_i) {
    auto& sip = *at(sips, s_i);
    for (auto& action : sip.actions) {
      if (action.action.kind == ACTION_SHIFT &&
          is_nonterminal(*grammar, *(action.context.begin()))) {
        auto nt = as_nonterminal(*grammar, *(action.context.begin()));
        add_nonterminal_action(out, s_i, nt, action.action.next_state);
      } else {
        for (auto terminal : action.context) {
          assert(is_terminal(*grammar, terminal));
          add_terminal_action(out, s_i, terminal, action.action);
        }
      }
    }
  }
  return out;
}

}
