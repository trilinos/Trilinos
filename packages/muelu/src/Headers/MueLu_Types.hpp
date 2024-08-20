// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TYPES_HPP
#define MUELU_TYPES_HPP

#include "MueLu_ConfigDefs.hpp"

namespace MueLu {
enum CycleType {
  VCYCLE,
  WCYCLE
};

enum PreOrPost {
  PRE  = 0x1,
  POST = 0x2,
  BOTH = 0x3
};

// In the algorithm, aggStat[] = READY/NOTSEL/SELECTED indicates whether a node has been aggregated
enum NodeState {
  READY = 1,  // indicates that a node is available to be
              // selected as a root node of an aggregate

  NOTSEL = 2,  // indicates that a node has been rejected as a root node.
               // This could perhaps be because if this node had been
               // selected a small aggregate would have resulted
               // This is Phase 1 specific

  AGGREGATED = 3,  // indicates that a node has been assigned
                   // to an aggregate

  ONEPT = 4,  // indicates that a node shall be preserved over
              // all multigrid levels as 1 point aggregate

  IGNORED = 5,  // indicates that the node is removed from consideration,
                // and is not aggregated

  BOUNDARY = 6,  // node is a Dirichlet node
                 // During aggregation, it is transformed either to AGGREGATED
                 // or to IGNORED
  INTERFACE = 7  // node is chosen as root node on an interface where coordinated
                 // coarsening across the interface is required.
};

// This is use by the structured aggregation index manager to keep track of the underlying mesh
// layout.
enum IndexingType {
  UNCOUPLED = 1,  // indicates that the underlying mesh is treated independently from rank to rank

  LOCALLEXI = 2,  // local lexicographic indexing of the mesh, this is similar to uncoupled but
                  // extra data is used to compute indices accross ranks

  GLOBALLEXI = 3  // global lexicographic indexing of the mesh means that the mesh is ordered
                  // lexicographically accorss and subsequently split among ranks.
};

}  // namespace MueLu

#endif  // ifndef MUELU_TYPES_HPP
