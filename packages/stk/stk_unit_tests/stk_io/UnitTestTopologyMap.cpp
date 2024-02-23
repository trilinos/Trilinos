// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_io/IossBridge.hpp>        // for map_ioss_topology_to_stk
#include <stk_topology/topology.hpp>    // for topology, etc
#include <gtest/gtest.h>
#include <string>                       // for operator<<, operator!=, etc
#include "gtest/gtest.h"                // for AssertHelper

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

int testElement(const std::string &name, unsigned spatialDim)
{
  int errors = 0;
  Ioss::ElementTopology *element = Ioss::ElementTopology::factory(name);
  if ( element == nullptr ) {
    // std::cerr << "\tERROR: Element type '" << name << "' could not be constructed.";
    // Must return since we have a NULL pointer and can't do further tests...
    return 1;
  }

  if (element->name() != name) {
    // This is an alias.  Don't run it through the rest of the tests
    // std::cerr << "Element '" << name << "' is an alias for element '" << element->name() << "'\n";
    return 0;
  }

  std::cerr << "Testing element '" << name << "'\n";
  // Currently not supported:
  if (element->name() == "unknown"   ||
      element->name() == "bar4"      ||
      element->name() == "edge4"     ||
      element->name() == "hex9"      ||
      element->name() == "hex16"     ||
      element->name() == "hex32"     ||
      element->name() == "hex64"     ||
      element->name() == "pyramid18" ||
      element->name() == "pyramid19" ||
      element->name() == "quad6"     ||
      element->name() == "quad12"    ||
      element->name() == "quad16"    ||
      element->name() == "tetra14"   ||
      element->name() == "tetra15"   ||
      element->name() == "tetra16"   ||
      element->name() == "tetra40"   ||
      element->name() == "tri13"     ||
      element->name() == "tri7"      ||
      element->name() == "tri9"      ||
      element->name() == "trishell7" ||
      element->name() == "wedge12"   ||
      element->name() == "wedge16"   ||
      element->name() == "wedge20"   ||
      element->name() == "wedge21"   ||
      element->name() == "wedge24"   ||
      element->name() == "wedge52") {
    std::cerr << "\tERROR (EXPECTED): No support for '" << element->name() << "'\n";
    return 0;
  }

  // Get the corresponding stk::topology:
  stk::topology cell = stk::io::map_ioss_topology_to_stk(element, spatialDim);
  if (cell == stk::topology::INVALID_TOPOLOGY) {
    std::cerr << "\tERROR: Could not find a stk::topology corresponding to the Ioss::ElementTopology element '" << name << "'.\n";
    return 1;
  }

  // See if we get the same element back when converting from
  // stk::topology to Ioss::ElementToplogy
  Ioss::ElementTopology *new_element = Ioss::ElementTopology::factory(cell.name());
  if (element->name() != new_element->name()) {
    std::cerr << "\tERROR: New name = '" << new_element->name() << "' doesn't match old name '" << element->name() << "'\n";
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

  // NOTE: stk::topology and Ioss disagree on parametric dimension.
  int add_to = element->spatial_dimension() != element->parametric_dimension() && element->is_element() ? 1 : 0;
  errors += my_assert(static_cast<int>(cell.dimension()),
                      element->parametric_dimension()+add_to,
                      "parametric dimension");
  errors += my_assert(static_cast<int>(cell.num_edges()),
                      element->number_edges(),
                      "edge count");

  // NOTE: Ioss counts edges and faces as boundaries for shell elements
  errors += my_assert(static_cast<int>(cell.num_sides()),
                      element->number_boundaries(),
                      "boundary count");

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
  unsigned spatialDim = 3;
  for (int i=0; i < element_count; i++) {
    // FIXME: Need to totally skip tetra7 for now
    if (elements[i] == "tetra7") {
      continue;
    }

    int current_error = testElement(elements[i], spatialDim);
    if (elements[i] != "node"    &&
        elements[i] != "bar2"    &&
        elements[i] != "bar3"    &&
        elements[i] != "bar4"    &&
        elements[i] != "spring2" &&
        elements[i] != "spring3") {
      errors += current_error;
    }
    else {
      if (current_error > 0)
        std::cerr << "\t\tIGNORING " << elements[i] << " ERRORS...\n";
    }
  }
  ASSERT_TRUE(errors == 0);
}
