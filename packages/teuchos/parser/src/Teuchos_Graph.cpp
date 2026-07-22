// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Graph.hpp"

#include <iostream>

#include "Teuchos_vector.hpp"

namespace Teuchos {

Graph make_graph_with_nnodes(int nnodes) {
  return Graph(std::size_t(nnodes));
}

int get_nnodes(Graph const& g) {
  return Teuchos::size(g);
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

int count_edges(const Graph& g, int i) {
  return Teuchos::size(at(g, i));
}

Graph make_transpose(Graph const& g) {
  int nnodes = get_nnodes(g);
  Graph transpose = make_graph_with_nnodes(nnodes);
  for (int i = 0; i < nnodes; ++i) {
    const NodeEdges& edges = get_edges(g, i);
    for (NodeEdges::const_iterator it = edges.begin(); it != edges.end(); ++it) {
      int j = *it;
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
    const NodeEdges& edges = get_edges(g, i);
    for (NodeEdges::const_iterator it = edges.begin(); it != edges.end(); ++it) {
      int j = *it;
      os << " " << j;
    }
    os << '\n';
  }
  return os;
}

}
