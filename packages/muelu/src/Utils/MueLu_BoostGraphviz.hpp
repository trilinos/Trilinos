// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_BOOSTGRAPHVIZ_HPP
#define MUELU_BOOSTGRAPHVIZ_HPP

// This header file can be used in place of <boost/graph/graphviz.hpp>. It disable the warnings present in boost.

// Note: pragma warnings available since gcc 4.2
//       pragma push/pop available since gcc 4.6
// We no longer check for gcc version as Trilinos requires a minimum 4.7.2.

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_BOOST) && defined(HAVE_MUELU_BOOST_FOR_REAL)

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif  // __GNUC__

#include <boost/graph/graphviz.hpp>

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif  // __GNUC__

// define boost graph types
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
                              boost::property<boost::vertex_name_t, std::string,
                                              boost::property<boost::vertex_color_t, std::string,
                                                              boost::property<boost::vertex_index_t, std::string> > >,
                              boost::property<boost::edge_name_t, std::string,
                                              boost::property<boost::edge_color_t, std::string> > >
    BoostGraph;
typedef boost::dynamic_properties BoostProperties;
typedef boost::graph_traits<BoostGraph>::vertex_descriptor BoostVertex;
typedef boost::graph_traits<BoostGraph>::edge_descriptor BoostEdge;

#endif  // HAVE_MUELU_BOOST && HAVE_MUELU_BOOST_FOR_REAL

#endif  // MUELU_BOOSTGRAPHVIZ_HPP
