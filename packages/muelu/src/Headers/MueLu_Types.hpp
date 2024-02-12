// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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
