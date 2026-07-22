// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_GRAPH_HPP
#define TEUCHOS_GRAPH_HPP

#include <vector>
#include <iosfwd>

namespace Teuchos {

typedef std::vector<int> NodeEdges;
typedef std::vector<NodeEdges> Graph;

Graph make_graph_with_nnodes(int nnodes);
int get_nnodes(Graph const& g);
void add_edge(Graph& g, int i, int j);
NodeEdges const& get_edges(Graph const& g, int i);
NodeEdges& get_edges(Graph& g, int i);
int count_edges(const Graph& g, int i);
Graph make_transpose(Graph const& g);
int at(Graph const& g, int i, int j);

std::ostream& operator<<(std::ostream& os, Graph const& g);

}

#endif
