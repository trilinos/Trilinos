// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <Ioss_ConcreteVariableType.h>  // for StorageInitializer
#include <Ioss_ElementTopology.h>       // for ElementTopology, NameList
#include <Ioss_Initializer.h>           // for Initializer
#include <stddef.h>                     // for NULL
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_io/IossBridge.hpp>        // for map_ioss_topology_to_stk
#include <stk_topology/topology.hpp>    // for topology, etc
#include <gtest/gtest.h>
#include <string>                       // for operator<<, operator!=, etc
#include "gtest/gtest.h"                // for AssertHelper
#include "stk_topology/topology.hpp"    // for topology::num_edges, etc





namespace {
template<typename T, typename U>
int my_assert(T a, T b, U msg) {
  if (a != b) {
    std::cerr << "\tERROR: '" << msg << "' assertion failed: " << a << " is not equal to " << b << "\n";
    return 1;
  }
  return 0;
}
}

namespace {
int testElement(const std::string &name)
{
  int errors = 0;
  Ioss::ElementTopology *element = Ioss::ElementTopology::factory(name);
  if (element == NULL) {
    std::cerr << "\tERROR: Element type '" << name << "' could not be constructed.";
    // Must return since we have a NULL pointer and can't do further tests...
    return 1;
  }
  if (element->name() != name) {
    // This is an alias.  Don't run it through the rest of the tests
    std::cerr << "Element '" << name << "' is an alias for element '" << element->name() << "'\n";
    return 0;
  }

  std::cerr << "Testing element '" << name << "'\n";
  // Currently not supported:
  if (element->name() == "unknown") {
    std::cerr << "\tERROR (EXPECTED): No support for '" << element->name() << "'\n";
    return 0;
  }

  // Get the corresponding stk::topology:
  stk::topology cell = stk::io::map_ioss_topology_to_stk(element);
  if (cell == stk::topology::INVALID_TOPOLOGY) {
    std::cerr << "\tERROR: Could not find a stk::topology corresponding to the Ioss::ElementTopology element '"
              << name << "'.";
    return 1;
  }

  // See if we get the same element back when converting from
  // stk::topology to Ioss::ElementToplogy
  Ioss::ElementTopology *new_element = Ioss::ElementTopology::factory(cell.name());
  if (element->name() != new_element->name()) {
    std::cerr << "\tERROR: New name = '" << new_element->name()
              << "' doesn't match old name '" << element->name()
              << "'\n";
    errors++;
  }

  // At this point, 'element' is the Ioss element topology and
  //                'cell' is the corresponding stk:topology topology.
  // Make sure that they agree on all subcell details...
  // Exceptions:
  // 1. An Ioss Node has 1 node per element; a stk::topology Node has 0 nodes per element...

  errors += my_assert(static_cast<int>(cell.num_nodes()),
                      element->number_nodes(),
                      "node count");
  errors += my_assert(static_cast<int>(cell.num_vertices()),
                      element->number_corner_nodes(),
                      "vertex count");

  // NOTE: CellTopology and Ioss disagree on parametric dimension.
  int add_to = element->spatial_dimension() != element->parametric_dimension() && element->is_element() ? 1 : 0;
  errors += my_assert(static_cast<int>(cell.dimension()),
                      element->parametric_dimension()+add_to,
                      "parametric dimension");
  errors += my_assert(static_cast<int>(cell.num_edges()),
                      element->number_edges(),
                      "edge count");

  // NOTE: Ioss counts edges and faces as boundaries for shell elements
  int add_boundary = 0;
  if (add_to == 1 && element->spatial_dimension() == 3 && element->parametric_dimension() == 2)
    add_boundary = cell.num_edges();
  if (cell == stk::topology::PARTICLE)
    add_boundary = -1;
  errors += my_assert(static_cast<int>(cell.num_sides()) + add_boundary,
                      element->number_boundaries(),
                      "boundary count");


#if 0
  // Check face topologies for all elements...
  if (element->is_element()) {
    if (cell.dimension() == 3) {
      int face_count = element->number_faces();
      for (int i=0; i < face_count; i++) {
        Ioss::ElementTopology *face = element->face_type(i+1);
        stk::topology cell_face = cell.face_topology(i);
        Ioss::ElementTopology *cell_face_top = Ioss::ElementTopology::factory(cell_face.name());
        errors += my_assert(face->name(),
                            cell_face_top->name()
                            "face type");

        Ioss::IntVector fcon = element->face_connectivity(i+1);
        size_t node_count = fcon.size();
        for (size_t j=0; j < node_count; j++) {
          std::ostringstream msg;
          msg << "face node connectivity for node " << j << " on face " << i;
          errors += my_assert(fcon[j],
                              static_cast<int>(cell.getNodeMap(cell.dimension()-1, i, j)),
                              msg.str());
        }
      }

      int edge_count = element->number_edges();
      for (int i=0; i < edge_count; i++) {

        Ioss::IntVector fcon = element->edge_connectivity(i+1);
        size_t node_count = fcon.size();
        for (size_t j=0; j < node_count; j++) {
          std::ostringstream msg;
          msg << "edge node connectivity for node " << j << " on edge " << i;
          errors += my_assert(fcon[j],
                              static_cast<int>(cell.getNodeMap(cell.dimension()-2, i, j)),
                              msg.str());
        }
      }
    }
    else if (cell.getDimension() == 2) {
      int edge_count = element->number_edges();
      for (int i=0; i < edge_count; i++) {
        Ioss::ElementTopology *edge = element->edge_type(i+1);
        const CellTopologyData *cell_edge = cell.getCellTopologyData(cell.dimension()-1,i);
        errors += my_assert(edge->name(),
                            stk::io::map_stk_topology_to_ioss(cell_edge, edge->spatial_dimension()),
                            "edge type");

        Ioss::IntVector econ = element->edge_connectivity(i+1);
        size_t node_count = econ.size();
        for (size_t j=0; j < node_count; j++) {
          std::ostringstream msg;
          msg << "edge node connectivity for node " << j << " on edge " << i;
          errors += my_assert(econ[j],
                              static_cast<int>(cell.getNodeMap(cell.getDimension()-1, i, j)),
                              msg.str());
        }
      }

    }
  }
#endif
  return errors;
}
}

TEST(UnitTestTopology, testUnit)
{
  Ioss::StorageInitializer initialize_storage;
  Ioss::Initializer        initialize_topologies;

  Ioss::NameList elements;
  int element_count = Ioss::ElementTopology::describe(&elements);

  int errors = 0;
  for (int i=0; i < element_count; i++) {
    // FIXME: Need to totally skip tetra7 for now
    if (elements[i] == "tetra7") {
      continue;
    }

    int current_error = testElement(elements[i]);
    if (elements[i] != "node" &&
        elements[i] != "rod3d2" &&
        elements[i] != "rod3d3" &&
        elements[i] != "tri4a" /*FIXME?*/) {
      errors += current_error;
    }
    else {
      if (current_error > 0)
        std::cerr << "\t\tIGNORING " << elements[i] << " ERRORS...\n";
    }
  }
  ASSERT_TRUE(errors == 0);
}


