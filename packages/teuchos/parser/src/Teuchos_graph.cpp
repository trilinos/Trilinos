#include "Teuchos_graph.hpp"

#include <iostream>

#include "Teuchos_vector.hpp"

namespace Teuchos {

Graph make_graph_with_nnodes(int nnodes) {
  return Graph(std::size_t(nnodes));
}

int get_nnodes(Graph const& g) {
  return size(g);
}

void add_edge(Graph& g, int i, int j) {
  at(g, i).push_back(j);
}

NodeEdges const& get_edges(Graph const& g, int i) {
  return at(g, i);
}

NodeEdges& get_edges(Graph& g, int i) {
  return at(g, i);
}

Graph make_transpose(Graph const& g) {
  auto nnodes = get_nnodes(g);
  auto transpose = make_graph_with_nnodes(nnodes);
  for (int i = 0; i < nnodes; ++i) {
    for (auto j : get_edges(g, i)) {
      add_edge(transpose, j, i);
    }
  }
  return transpose;
}

int at(Graph const& g, int i, int j) {
  return at(at(g, i), j);
}

std::ostream& operator<<(std::ostream& os, Graph const& g) {
  for (int i = 0; i < get_nnodes(g); ++i) {
    os << i << ":";
    for (auto j : get_edges(g, i)) os << " " << j;
    os << '\n';
  }
  return os;
}

}
