// Copyright(C) 1999-2020, 2022, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"

#include "Partition.h"
#include "PartitionTypes.h"

#include <cstdint> // for int64_t
#include <cstdio>  // for stderr, etc
#include <cstdlib> // for exit
#include <cstring>

#include <cstdarg> // for va_end, va_arg, va_list, etc
#include <cstddef> // for size_t
#include <cstdio>  // for stderr
#include <cstdlib> // for exit, malloc
#include <fmt/format.h>
#include <fmt/ostream.h>

#include "DataUtils.h"
#include "IossMesh.h"

namespace utest_util {

  template void fill_node_element_connectivity(IossMesh                      *mesh,
                                               std::vector<std::vector<int>> &surroundingElements,
                                               int &maxNumSurroundingElements);
  template void
  fill_node_element_connectivity(IossMesh                          *mesh,
                                 std::vector<std::vector<int64_t>> &surroundingElements,
                                 int                               &maxNumSurroundingElements);

  template <typename INT>
  void fill_node_element_connectivity(IossMesh                      *mesh,
                                      std::vector<std::vector<INT>> &surroundingElements,
                                      int                           &maxNumSurroundingElements)
  {
    surroundingElements.clear();
    maxNumSurroundingElements = 0;

    size_t numNodes = mesh->get_num_local_nodes();
    size_t numElems = mesh->get_num_local_elements();

    std::vector<int>    surround_count(numNodes, 0);
    std::vector<size_t> last_element(numNodes, 0);

    size_t sur_elem_total_size = 0;
    /* Find the count of surrounding elements for each node in the mesh */
    // The hope is that this code will speed up the entire routine even
    // though we are iterating the complete connectivity array twice.
    for (size_t ecnt = 0; ecnt < numElems; ecnt++) {
      IossElementData elem = mesh->get_local_element(ecnt);

      int nnodes = elem.topology->number_nodes();

      for (int ncnt = 0; ncnt < nnodes; ncnt++) {
        auto node = elem.localConnectivity[ncnt];
        /*
         * in the case of degenerate elements, where a node can be
         * entered into the connect table twice, need to check to
         * make sure that this element is not already listed as
         * surrounding this node
         */
        if (last_element[node] != ecnt + 1) {
          last_element[node] = ecnt + 1;
          surround_count[node]++;
          sur_elem_total_size++;
        }
      }
    }

    vec_free(last_element);

    // Attempt to reserve an array with this size...
    surroundingElements.resize(numNodes);
    for (size_t ncnt = 0; ncnt < numNodes; ncnt++) {
      if (surround_count[ncnt] == 0) {
        fmt::print(stderr, "WARNING: Node = {} has no elements\n", ncnt + 1);
      }
      else {
        surroundingElements[ncnt].reserve(surround_count[ncnt]);
        maxNumSurroundingElements = std::max(surround_count[ncnt], maxNumSurroundingElements);
      }
    }

    /* Find the surrounding elements for each node in the mesh */
    for (size_t ecnt = 0; ecnt < numElems; ecnt++) {
      IossElementData elem = mesh->get_local_element(ecnt);

      int nnodes = elem.topology->number_nodes();
      for (int ncnt = 0; ncnt < nnodes; ncnt++) {
        auto node = elem.localConnectivity[ncnt];

        /*
         * in the case of degenerate elements, where a node can be
         * entered into the connect table twice, need to check to
         * make sure that this element is not already listed as
         * surrounding this node
         */
        if (surroundingElements[node].empty() || ecnt != (size_t)surroundingElements[node].back()) {
          /* Add the element to the list */
          surroundingElements[node].push_back(ecnt);
        }
      }
    }

#ifndef NDEBUG
    for (size_t ncnt = 0; ncnt < numNodes; ncnt++) {
      assert(surroundingElements[ncnt].size() == (size_t)surround_count[ncnt]);
    }
#endif
  }

  Partition::Partition(IossMesh *mesh, const std::vector<EntityProc> &procAssign, const int nProc)
      : m_mesh(mesh), m_numProcs(nProc), m_entityPartition(procAssign)
  {
  }

  void Partition::verify_input_partition()
  {
    EXPECT_TRUE(m_numProcs > 0) << "Invalid processor decomposition count: " << m_numProcs;

    // Verify processor distribution
    for (const auto &entry : m_entityPartition) {
      EXPECT_TRUE(entry.proc >= 0 && entry.proc < m_numProcs)
          << "Invalid processor assignment for entity " << entry.id;
    }

    // Sort and uniquify
    std::sort(m_entityPartition.begin(), m_entityPartition.end(), EntityProcLess());
    auto iter = std::unique(m_entityPartition.begin(), m_entityPartition.end(), EntityProcEqual());
    m_entityPartition.resize(iter - m_entityPartition.begin());

    EXPECT_EQ(get_num_local_entities(), m_entityPartition.size());

    // Make sure every node has an assignment
    for (size_t i = 0; i < get_num_local_entities(); i++) {
      (void)get_processor(i);
    }
  }
} // namespace utest_util
