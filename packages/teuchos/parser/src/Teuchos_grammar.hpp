#ifndef TEUCHOS_GRAMMAR_HPP
#define TEUCHOS_GRAMMAR_HPP

#include <string>
#include <vector>
#include <memory>

namespace Teuchos {

/* convention: symbols are numbered with all
   terminal symbols first, all non-terminal symbols after */

struct Grammar {
  using RHS = std::vector<int>;
  struct Production {
    int lhs;
    RHS rhs;
  };
  using Productions = std::vector<Production>;
  int nsymbols;
  int nterminals;
  Productions productions;
  std::vector<std::string> symbol_names;
};

using GrammarPtr = std::shared_ptr<Grammar const>;

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
