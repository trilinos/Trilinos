// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_BUILD_PARSER_HPP
#define TEUCHOS_BUILD_PARSER_HPP

#include <set>
#include <memory>

#include <Teuchos_Parser.hpp>
#include <Teuchos_Graph.hpp>
#include <Teuchos_RCP.hpp>

namespace Teuchos {

struct Config {
  int production;
  int dot;
  Config(int p, int d);
};

typedef std::vector<Config> Configs;

typedef std::set<int> Context;

/* nonterminal transitions will be stored as SHIFT
   actions while in progress */
struct ActionInProgress {
  Action action;
  Context context;
};

struct StateInProgress {
  std::vector<int> configs;
  typedef std::vector<ActionInProgress> Actions;
  Actions actions;
};

void swap(StateInProgress& a, StateInProgress& b);

typedef RCP<StateInProgress> StateInProgressPtr;

typedef std::vector<StateInProgressPtr> StatesInProgress;

struct StateConfig {
  int state;
  int config_in_state;
  StateConfig(int s, int cis);
};

typedef std::vector<StateConfig> StateConfigs;

struct ParserInProgress {
  StatesInProgress states;
  Configs configs;
  StateConfigs state_configs;
  Graph states2state_configs;
  GrammarPtr grammar;
};

StateConfigs form_state_configs(StatesInProgress const& states);
Graph form_states_to_state_configs(StateConfigs const& scs,
    StatesInProgress const& states);

void print_graphviz(
    std::string const& filepath,
    ParserInProgress const& pip,
    bool verbose,
    std::ostream& os
    );

ParserInProgress draft_lalr1_parser(GrammarPtr grammar, bool verbose = false);

Parser accept_parser(ParserInProgress const& pip);

class ParserBuildFail: public std::invalid_argument {
 public:
  ParserBuildFail(const std::string& msg);
};

}

#endif
